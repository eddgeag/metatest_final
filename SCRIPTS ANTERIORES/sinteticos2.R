# =======================================================
# Paquetes
# =======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet)         # ridge/lasso/elastic-net
  library(nnet)           # multinomial
  library(MASS)           # mvrnorm + stepAIC
  library(caret)          # confusionMatrix + train + createFolds
  library(pROC)           # multiclass.roc
  library(e1071)          # svm (no usado aquí)
  library(ggplot2)
  library(ggridges)
})

set.seed(123)

# =======================================================
# Utilidades de métricas y calibración
# =======================================================
metrics_multiclass <- function(y_true, probs_mat) {
  out_na <- c(Accuracy=NA_real_, Sensibilidad=NA_real_, Especificidad=NA_real_, AUC=NA_real_)
  if (is.null(probs_mat) || nrow(probs_mat) == 0) return(out_na)
  y_use <- y_true; p_use <- probs_mat
  if (!is.null(names(y_true)) && !is.null(rownames(probs_mat))) {
    common <- intersect(names(y_true), rownames(probs_mat))
    if (length(common) >= 1) {
      y_use <- factor(y_true[common], levels=levels(y_true))
      p_use <- probs_mat[common, , drop=FALSE]
    }
  }
  if (length(y_use) != nrow(p_use)) {
    n <- min(length(y_use), nrow(p_use))
    y_use <- y_use[seq_len(n)]
    p_use <- p_use[seq_len(n), , drop=FALSE]
  }
  if (any(!is.finite(p_use))) return(out_na)
  rs <- rowSums(p_use); bad <- !is.finite(rs) | rs == 0
  if (any(bad)) p_use[bad, ] <- 1 / ncol(p_use)
  cls  <- colnames(p_use)
  pred <- factor(cls[max.col(p_use, ties.method="first")], levels=levels(y_true))
  if (length(pred) != length(y_use)) return(out_na)
  cm   <- caret::confusionMatrix(pred, y_use)
  acc  <- as.numeric(cm$overall["Accuracy"])
  byc  <- suppressWarnings(cm$byClass)
  sens <- mean(byc[,"Sensitivity"], na.rm=TRUE)
  spec <- mean(byc[,"Specificity"], na.rm=TRUE)
  auc  <- tryCatch(as.numeric(pROC::multiclass.roc(y_use, p_use)$auc), error=function(e) NA_real_)
  c(Accuracy=acc, Sensibilidad=sens, Especificidad=spec, AUC=auc)
}

platt_fit <- function(y_true, P){
  lev <- colnames(P)
  fits <- lapply(lev, function(lv){
    x <- as.numeric(P[,lv]); y <- as.integer(y_true==lv)
    o <- order(x); x <- x[o]; y <- y[o]
    iso <- isoreg(x, y)
    list(x = iso$x, y = iso$yf)
  })
  setNames(fits, lev)
}
platt_apply <- function(fits, P){
  Z <- sapply(names(fits), function(lv){
    f <- fits[[lv]]; x <- as.numeric(P[,lv])
    yhat <- approx(f$x, f$y, x, rule=2, ties="ordered")$y
    pmin(pmax(yhat, 1e-6), 1-1e-6)
  })
  Z <- sweep(Z, 1, rowSums(Z), "/")
  Z
}

# =======================================================
# Inputs
# =======================================================
dats      <- readRDS("./datos_para_modelar.rds")
tst_Data  <- readRDS("./ready_for_modeling.rds")

factors    <- dats$factores
pesos_list <- dats$pesos

# === 1) Semilla: "Real" (tu 60%) ===
Z_seed <- as.matrix(factors %>% filter(Tipo=="Real") %>% dplyr::select(where(is.numeric)))[,1:2,drop=FALSE]
y_seed <- factors %>% filter(Tipo=="Real") %>% pull(Grupo) %>% factor()
names(y_seed) <- rownames(Z_seed)  # clave para alineación
lev <- levels(y_seed)

# --- split de semilla: anchors (train) y val (holdout) ---
split_seed <- function(y, k_anchor = 2L){
  L <- levels(y); ids_all <- names(y)
  # vectorizar
  if (length(k_anchor) > 1L) {
    stopifnot(!is.null(names(k_anchor)), all(names(k_anchor) %in% L))
    k_vec <- as.integer(k_anchor[L])
  } else {
    k_vec <- rep(as.integer(k_anchor), length(L))
  }
  # muestreo por clase dejando al menos 1 para validación
  out <- unlist(mapply(function(cl, k){
    ids_cl <- ids_all[y==cl]
    if (length(ids_cl)==0) return(character(0))
    k <- min(max(0L, k), max(0L, length(ids_cl)-1L))
    if (k==0) return(character(0))
    sample(ids_cl, k, replace=FALSE)
  }, L, k_vec, SIMPLIFY=FALSE), use.names=FALSE)
  unique(out)
}



# --- anchors y holdout coherentes ---
tbl <- table(y_seed)
n_i <- as.integer(tbl)

desired   <- pmax(5L, floor(0.2 * n_i))                  # objetivo por clase
k_anchor  <- setNames(pmax(1L, pmin(desired, pmax(0L, n_i - 1L))), names(tbl))

ids_anchor <- split_seed(y_seed, k_anchor = k_anchor)     # usa la versión vectorizada
stopifnot(length(ids_anchor) > 0)

ids_val <- setdiff(names(y_seed), ids_anchor)
stopifnot(length(ids_val) > 0)

Z_fit <- Z_seed[ids_anchor, , drop = FALSE]
y_fit <- y_seed[ids_anchor]
Z_val <- Z_seed[ids_val,   , drop = FALSE]
y_val <- y_seed[ids_val]

# --- objetivos con mínimo por fold ---
k <- 5L; min_fold <- 8L
min_total  <- ceiling(min_fold * k / (k - 1))   # >=10
max_ratio  <- 3L
base_target <- 15L

raw_target <- pmin(base_target, n_i * max_ratio)
targets    <- setNames(pmax(raw_target, min_total), names(tbl))  # mantener nombres

print(k_anchor)
print(targets)

# síntesis usando SOLO anchors para estimar; n_add = target - k_anchor
gen_class <- function(Zc, n, shrink=0.4, rad_q=0.90){
  if (n <= 0) return(NULL)
  p <- ncol(Zc)
  mu <- colMeans(Zc)
  if (nrow(Zc) < 2) {
    # covarianza artificial diagonal si hay 1 punto
    S <- diag(rep(1, p))
  } else {
    S <- cov(Zc)
    S <- (1-shrink)*S + shrink*diag(diag(S))
  }
  X <- MASS::mvrnorm(n, mu=mu, Sigma=S)
  if (nrow(Zc) >= 3) {
    D2 <- mahalanobis(X, center=mu, cov=S)
    thr <- quantile(mahalanobis(Zc, mu, S), rad_q)
    idx <- which(D2 > thr)
    if (length(idx)) X[idx,] <- sweep(X[idx,], 1, sqrt(thr/D2[idx]), `*`)
  }
  X
}
project_to_latent <- function(X, W) {
  G <- crossprod(W)
  as.matrix(X) %*% W %*% MASS::ginv(G)
}
synth_Z2 <- function(Z_fit, y_fit, targets, shrink=0.2){
  lev <- levels(y_fit)
  Zs <- list(); ys <- c()
  for (lv in lev){
    Zi <- Z_fit[y_fit==lv,,drop=FALSE]
    need <- targets[[lv]] - nrow(Zi)
    if (need > 0){
      newZ <- gen_class(Zi, need, shrink=shrink)
      if (!is.null(newZ)){
        Zs[[lv]] <- newZ
        ys <- c(ys, rep(lv, nrow(newZ)))
      }
    }
  }
  Zs <- if (length(Zs)) do.call(rbind, Zs) else NULL
  list(Z=Zs, y=factor(ys, levels=lev))
}

gen <- synth_Z2(Z_fit, y_fit, targets, shrink=0.3)
str(gen)
stopifnot(!is.null(gen$Z))

# train/val finales SIN fuga
train_latent_X <- rbind(Z_fit, gen$Z)
train_y        <- factor(c(as.character(y_fit), as.character(gen$y)), levels = levels(y_seed))
val_latent_X   <- Z_val
val_y          <- y_val

# =======================================================
# Helpers de vistas + MOFA (W) + proyección/reconstrucción
# =======================================================
orient_and_name <- function(X, ids) {
  if (is.null(ids)) return(X)
  ids <- as.character(ids)
  if (nrow(X) != length(ids) && ncol(X) == length(ids)) X <- t(X)
  if (nrow(X) == length(ids)) rownames(X) <- ids
  X
}

view_map <- c(transcriptomica ="tx", proteomica="pr", metabolomica ="me", clinical = "cl")

bind_train_test <- function(vname) {
  obj <- tst_Data[[view_map[[vname]]]]
  Xt <- obj$train;  Xs <- obj$test
  Xt <- orient_and_name(Xt, obj$ids_train)
  Xs <- orient_and_name(Xs, obj$ids_test)
  common_feats <- intersect(colnames(Xt), colnames(Xs))
  if (length(common_feats) == 0) stop(sprintf("Vista %s sin intersección de features", vname))
  Xt <- Xt[, common_feats, drop=FALSE]
  Xs <- Xs[, common_feats, drop=FALSE]
  list(X_train = Xt, X_test = Xs, X_all = rbind(Xt, Xs))
}

val_X_list  <- setNames(vector("list", length(pesos_list)), names(pesos_list))
test_X_list <- setNames(vector("list", length(pesos_list)), names(pesos_list))

ids_val <- rownames(val_latent_X); stopifnot(!is.null(ids_val))
for (v in names(pesos_list)) {
  bt <- bind_train_test(v)
  if (all(ids_val %in% rownames(bt$X_all))) {
    val_X_list[[v]] <- bt$X_all[ids_val, , drop=FALSE]
  } else {
    val_X_list[[v]] <- bt$X_test
  }
  test_X_list[[v]] <- bt$X_test
}

# Renombrar transcriptómica: EntrezID -> Símbolo
if (!is.null(tst_Data$features_metadata)) {
  dict_tx <- tst_Data$features_metadata %>%
    dplyr::select(EntrezGeneID, GeneSymbol) %>% distinct() %>%
    filter(!is.na(EntrezGeneID), !is.na(GeneSymbol))
  map_entrez_to_symbol <- function(nms, dict) {
    nms_new <- dict$GeneSymbol[match(nms, dict$EntrezGeneID)]
    nms_final <- ifelse(is.na(nms_new) | nms_new=="", nms, nms_new)
    make.unique(nms_final)
  }
  if (!is.null(val_X_list$transcriptomica)) {
    colnames(val_X_list$transcriptomica) <- map_entrez_to_symbol(colnames(val_X_list$transcriptomica), dict_tx)
  }
  if (!is.null(test_X_list$transcriptomica)) {
    colnames(test_X_list$transcriptomica) <- map_entrez_to_symbol(colnames(test_X_list$transcriptomica), dict_tx)
  }
  if ("transcriptomica" %in% names(pesos_list)) {
    rownames(pesos_list$transcriptomica) <- map_entrez_to_symbol(rownames(pesos_list$transcriptomica), dict_tx)
  }
}

mk_W <- function(pesos_list) Reduce("rbind", pesos_list)

cbind_views <- function(X_list, pesos_list) {
  ord <- names(pesos_list)
  mats <- lapply(ord, function(v){
    Xi <- X_list[[v]]
    if (is.null(Xi)) stop(sprintf("Vista %s vacía", v))
    if (nrow(Xi) == nrow(pesos_list[[v]]) && ncol(Xi) != nrow(pesos_list[[v]])) Xi <- t(Xi)
    feat_order <- rownames(pesos_list[[v]])
    Xfull <- matrix(0, nrow=nrow(Xi), ncol=length(feat_order), dimnames=list(rownames(Xi), feat_order))
    common <- intersect(colnames(Xi), feat_order)
    if (length(common) > 0) Xfull[, common] <- as.matrix(Xi[, common, drop=FALSE])
    Xfull
  })
  nr <- sapply(mats, nrow); stopifnot(length(unique(nr)) == 1)
  do.call(cbind, mats)
}

reconstruct_from_latent <- function(Z, W) as.matrix(Z) %*% t(W)
project_to_latent <- function(X, W) { G <- crossprod(W); as.matrix(X) %*% W %*% solve(G) }

W <- mk_W(pesos_list); stopifnot(ncol(W)==2)

# Matrices crudo (no integradas) y latentes
X_va_raw <- cbind_views(val_X_list,  pesos_list)
X_te_raw <- cbind_views(test_X_list, pesos_list)
X_te_lat <- project_to_latent(X_te_raw, W)

X_tr_lat <- as.matrix(train_latent_X)
X_va_lat <- as.matrix(val_latent_X)

# Reconstrucciones SIN escalar (espacio de features)
X_tr_rec0 <- reconstruct_from_latent(X_tr_lat, W)
X_va_rec0 <- reconstruct_from_latent(X_va_lat,   W)
X_te_rec0 <- reconstruct_from_latent(X_te_lat,   W)

# Escalado coherente usando SOLO el TRAIN reconstruido
mu_rec  <- colMeans(X_tr_rec0)
sd_rec  <- apply(X_tr_rec0, 2, sd); sd_rec[sd_rec==0] <- 1
scale_with_train <- function(X, mu, sds){
  Xs <- sweep(X, 2, mu, "-")
  sweep(Xs, 2, sds, "/")
}
X_tr_rec <- scale_with_train(X_tr_rec0, mu_rec, sd_rec)
X_va_rec <- scale_with_train(X_va_rec0, mu_rec, sd_rec)
X_te_rec <- scale_with_train(X_te_rec0, mu_rec, sd_rec)

# =======================================================
# Visualización PCA (en latente MOFA 2D)
# =======================================================
df_train <- data.frame(X_tr_lat, Grupo=train_y, Set="Train_sint")
df_val   <- data.frame(X_va_lat, Grupo=val_y,   Set="Val_real")
df_test  <- data.frame(X_te_lat, Grupo=tst_Data$grupo_test$grupo, Set="Test_real")

df_all <- bind_rows(df_train, df_val, df_test)

pca_all <- prcomp(df_all %>% dplyr::select(where(is.numeric)), scale.=FALSE)
scores <- as.data.frame(pca_all$x[,1:2])
scores$Grupo <- df_all$Grupo
scores$Set   <- df_all$Set

# Ellipses solo si hay >=3 puntos por grupo/set
scores_ell <- scores %>% dplyr::group_by(Grupo, Set) %>% dplyr::filter(dplyr::n() >= 3) %>% dplyr::ungroup()

# Scatter global
print(
  ggplot(scores, aes(PC1, PC2, color=Grupo, shape=Set)) +
    geom_point(alpha=0.7, size=3) +
    stat_ellipse(data=scores_ell, aes(group=interaction(Grupo, Set)), linetype=2, alpha=0.3) +
    theme_minimal(base_size=14) +
    labs(title="PCA de conjuntos sintéticos vs reales (MOFA 2D)",
         subtitle="Scoreplot en PC1–PC2",
         color="Grupo", shape="Conjunto")
)

# Densidades por PC y grupo
print(
  ggplot(scores, aes(x=PC1, fill=Set, color=Set)) +
    geom_density(alpha=0.3) +
    facet_wrap(~Grupo, scales="free") +
    theme_minimal(base_size=14) +
    labs(title="Distribución de PC1 por Grupo",
         subtitle="Comparación entre Train_sint, Val_real y Test_real")
)
print(
  ggplot(scores, aes(x=PC2, fill=Set, color=Set)) +
    geom_density(alpha=0.3) +
    facet_wrap(~Grupo, scales="free") +
    theme_minimal(base_size=14) +
    labs(title="Distribución de PC2 por Grupo",
         subtitle="Comparación entre Train_sint, Val_real y Test_real")
)

# Ridges
print(
  ggplot(scores, aes(x=PC1, y=Set, fill=Set)) +
    geom_density_ridges(alpha=0.5) +
    facet_wrap(~Grupo, scales="free_x") +
    theme_minimal(base_size=14) +
    labs(title="Ridge: PC1 por grupo")
)
print(
  ggplot(scores, aes(x=PC2, y=Set, fill=Set)) +
    geom_density_ridges(alpha=0.5) +
    facet_wrap(~Grupo, scales="free_x") +
    theme_minimal(base_size=14) +
    labs(title="Ridge: PC2 por grupo")
)

# =======================================================
# Lasso multinomial entrenado con sintéticos balanceados
# =======================================================
# Reutiliza X_tr_rec/X_va_rec/X_te_rec ya escalados con el TRAIN
set.seed(123)
foldid <- caret::createFolds(train_y, k=5, list=FALSE)
cvfit_lasso <- cv.glmnet(
  x=as.matrix(X_tr_rec), y=train_y,
  family="multinomial", alpha=1,
  type.multinomial="ungrouped",
  foldid=foldid, standardize=FALSE
)

lambda_sel <- cvfit_lasso$lambda.1se
cat("Lambda seleccionado:", lambda_sel, "\n")

# Predicciones integradas
m_lasso_tr  <- predict(cvfit_lasso, newx = X_tr_rec, s=lambda_sel, type="response")[,,1]
m_lasso_va  <- predict(cvfit_lasso, newx = X_va_rec, s=lambda_sel, type="response")[,,1]
m_lasso_te  <- predict(cvfit_lasso, newx = X_te_rec, s=lambda_sel, type="response")[,,1]

# Predicciones no integradas (opcional)
X_va_raw <- cbind_views(val_X_list,  pesos_list)  # reasegurar en caso de GC
X_te_raw <- cbind_views(test_X_list, pesos_list)
m_lasso_vaN <- predict(cvfit_lasso, newx = X_va_raw, s=lambda_sel, type="response")[,,1]
m_lasso_teN <- predict(cvfit_lasso, newx = X_te_raw, s=lambda_sel, type="response")[,,1]

# Etiquetas de test (si existen)
y_test <- tst_Data$grupo_test$grupo
names(y_test) <- rownames(X_te_rec)

# priores (entrenaste balanceado; objetivo ~ val)
lev <- levels(train_y)
pi_train <- prop.table(table(train_y))[lev]
pi_tgt   <- prop.table(table(val_y))[lev]   # o tus priores esperados

adjust_priors <- function(P, pi_tr, pi_tg){
  Q <- sweep(P, 2, as.numeric(pi_tr), "/")
  Q <- sweep(Q, 2, as.numeric(pi_tg),  "*")
  sweep(Q, 1, rowSums(Q), "/")
}

# tras predecir en test integrado:
m_te_adj <- adjust_priors(m_lasso_te,  pi_train, pi_tgt)

# Métricas crudas
print(metrics_multiclass(train_y, m_lasso_tr))
print(metrics_multiclass(val_y,   m_lasso_va))
print(metrics_multiclass(val_y,   m_lasso_vaN))
print(metrics_multiclass(y_test,  m_lasso_te))
print(metrics_multiclass(y_test,  m_lasso_teN))
print(metrics_multiclass(y_test,  m_te_adj))

# Calibración Platt con VALIDACIÓN integrada -> aplicar en TEST
platt <- platt_fit(val_y, m_lasso_va)
m_lasso_te_cal  <- platt_apply(platt, m_lasso_te)
m_lasso_teN_cal <- platt_apply(platt, m_lasso_teN)

cat("\n== Métricas TEST calibrado (integrado) ==\n")
print(metrics_multiclass(y_test, m_lasso_te_cal))

library(caret)

# clases posibles
lev <- colnames(m_lasso_te_cal)

# predicciones duras
pred_test <- factor(lev[max.col(m_lasso_te_cal, ties.method="first")],
                    levels = levels(y_test))

# matriz de confusión
cm <- confusionMatrix(pred_test, y_test)
cm
cat("\n== Métricas TEST calibrado (no integrado) ==\n")
print(metrics_multiclass(y_test, m_lasso_teN_cal))
brier <- function(y, P){ mean(rowSums((model.matrix(~y-1) - P)^2)) }
ece <- function(y, P, M=10){
  y1 <- apply(model.matrix(~y-1),1,which.max); pmax <- apply(P,1,max)
  brks <- quantile(pmax, probs=seq(0,1,length.out=M+1))
  idx <- cut(pmax, brks, include.lowest=TRUE)
  m <- tapply(pmax, idx, mean); a <- tapply(as.integer(y), idx, function(v) mean(v==names(which.max(table(y)))))
  sum(abs(m - a) * table(idx))/length(pmax)
}
cat("Brier test:", brier(y_test, m_lasso_te), "\n")
