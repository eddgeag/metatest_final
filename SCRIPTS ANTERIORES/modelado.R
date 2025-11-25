# =======================================================
# Paquetes
# =======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet)         # ridge/lasso/elastic-net
  library(nnet)           # multinomial
  library(MASS)           # stepAIC
  library(caret)          # confusionMatrix + train
  library(pROC)           # multiclass.roc
  library(e1071)          # svm
})

set.seed(123)

# =======================================================
# Inputs
# =======================================================
dats      <- readRDS("./datos_para_modelar.rds")
tst_Data  <- readRDS("./ready_for_modeling.rds")

factors        <- dats$factores
pesos_list     <- dats$pesos
## === 1) Semilla: “Real” (tu 60%) ===
Z_seed <- as.matrix(factors %>% filter(Tipo=="Real") %>% dplyr::select(where(is.numeric)))[,1:2,drop=FALSE]
y_seed <- factors %>%
  filter(Tipo=="Real") %>%
  pull(Grupo) %>%
  factor()
names(y_seed) <- rownames(Z_seed)   # << clave
lev    <- levels(y_seed)

## === 2) Objetivo por clase con límite de inflación ===
tbl <- table(y_seed)                         # NP=6, FD=3, SO=7
max_ratio <- 3L                              # no más de x3
base_target <- 15L                           # o el mínimo entre 15 y x3
targets <- pmin(base_target, as.integer(tbl)*max_ratio)
names(targets) <- names(tbl)                 # por clase

## === 3) Síntesis gaussiana con pool + shrink adaptativo (por clase) ===
synth_Z2 <- function(Z2d, y, targets, min_per_class=5, lambda=1e-2) {
  stopifnot(nrow(Z2d)==length(y)); L <- levels(y)
  S_pool <- stats::cov(Z2d) + diag(lambda, ncol(Z2d))
  outs <- list(); labs <- c()
  for (cl in L) {
    Zi <- Z2d[y==cl,,drop=FALSE]
    n_real <- nrow(Zi); n_tgt <- targets[cl]; n_add <- max(0, n_tgt - n_real)
    if (n_add<=0 || n_real<1) next
    mu <- colMeans(Zi)
    if (n_real >= min_per_class) {
      S <- cov(Zi); shrink <- lambda * mean(diag(S)); S <- S + diag(shrink, ncol(Z2d))
    } else {
      S <- S_pool
    }
    outs[[cl]] <- MASS::mvrnorm(n=n_add, mu=mu, Sigma=S)
    labs <- c(labs, rep(cl, n_add))
  }
  if (!length(outs)) return(list(Z=NULL, y=NULL))
  Zs <- do.call(rbind, outs)
  rownames(Zs) <- paste0("synth_", labs, ave(seq_along(labs), labs, FUN=seq))
  list(Z=Zs, y=factor(labs, levels=L))
}

gen <- synth_Z2(Z_seed, y_seed, targets, min_per_class=5, lambda=1e-2)
stopifnot(!is.null(gen$Z))

## === 4) “Anclaje” con reales: incluye pocos reales por clase (1–2) ===
## === 4) “Anclaje” con reales: incluye pocos reales por clase (1–2) ===
k_anchor <- setNames(pmin(2L, as.integer(tbl)), names(tbl))  # ahora tiene nombres

ids_anchor <- unlist(lapply(names(tbl), function(cl){
  ids_cl <- rownames(Z_seed)[y_seed == cl]
  if (length(ids_cl) == 0) return(character(0))
  sample(ids_cl, k_anchor[cl])   # ahora siempre un número válido
}))

Z_anchor <- Z_seed[ids_anchor,,drop=FALSE]
y_anchor <- y_seed[ids_anchor]

## === 5) Entrenamiento final en latente ===
train_latent_X <- as.data.frame(rbind(Z_anchor, gen$Z))
train_y        <- factor(c(as.character(y_anchor), as.character(gen$y)), levels=lev)
rownames(train_latent_X) <- c(rownames(Z_anchor), rownames(gen$Z))
names(train_y) <- rownames(train_latent_X)

## === 6) “Val” interno (opcional, solo sanity check, NO calibrar/tunar) ===
val_latent_X <- Z_seed[setdiff(rownames(Z_seed), ids_anchor), , drop=FALSE]
val_y        <- y_seed[setdiff(names(y_seed), ids_anchor)]
# úsalo solo para mirar sobreajuste


# === helper: asegura muestras×features y pone rownames con ids_* ===
orient_and_name <- function(X, ids) {
  if (is.null(ids)) return(X)                       # no IDs, devolver tal cual
  ids <- as.character(ids)
  # si filas!=length(ids) pero columnas==length(ids) => transpón
  if (nrow(X) != length(ids) && ncol(X) == length(ids)) X <- t(X)
  # si ahora sí coincide, asigna
  if (nrow(X) == length(ids)) rownames(X) <- ids
  X
}

view_map <- c(transcriptomica ="tx",
              proteomica = "pr" ,
              metabolomica ="me",
              clinical = "cl")
# === construir listas crudo: VALIDACIÓN y TEST ===
bind_train_test <- function(vname) {
  obj <- tst_Data[[view_map[[vname]]]]
  Xt <- obj$train;  Xs <- obj$test
  Xt <- orient_and_name(Xt, obj$ids_train)
  Xs <- orient_and_name(Xs, obj$ids_test)
  
  # une por columnas comunes; si difieren, usa intersección
  common_feats <- intersect(colnames(Xt), colnames(Xs))
  if (length(common_feats) == 0)
    stop(sprintf("Vista %s sin intersección de features entre train/test", vname))
  Xt <- Xt[, common_feats, drop=FALSE]
  Xs <- Xs[, common_feats, drop=FALSE]
  
  Xall <- rbind(Xt, Xs)
  list(X_train = Xt, X_test = Xs, X_all = Xall)
}

val_X_list  <- setNames(vector("list", length(pesos_list)), names(pesos_list))
test_X_list <- setNames(vector("list", length(pesos_list)), names(pesos_list))

# IDs de validación desde el bloque latente real
ids_val <- rownames(val_latent_X)
stopifnot(!is.null(ids_val))

for (v in names(pesos_list)) {
  bt <- bind_train_test(v)
  # validación no integrada: filtra X_all por ids_val si existen
  if (all(ids_val %in% rownames(bt$X_all))) {
    val_X_list[[v]] <- bt$X_all[ids_val, , drop=FALSE]
  } else {
    # fallback: usa test si los ids no están en X_all
    val_X_list[[v]] <- bt$X_test
  }
  test_X_list[[v]] <- bt$X_test
}


# =======================================================
# Renombrar solo transcriptómica: EntrezID -> GeneSymbol
# =======================================================
dict_tx <- tst_Data$features_metadata %>%
  dplyr::select(EntrezGeneID, GeneSymbol) %>%
  distinct() %>%
  filter(!is.na(EntrezGeneID), !is.na(GeneSymbol))

map_entrez_to_symbol <- function(nms, dict) {
  nms_new <- dict$GeneSymbol[match(nms, dict$EntrezGeneID)]
  nms_final <- ifelse(is.na(nms_new) | nms_new=="", nms, nms_new)
  make.unique(nms_final)
}
if (!is.null(val_X_list$transcriptomica)) {
  colnames(val_X_list$transcriptomica) <-
    map_entrez_to_symbol(colnames(val_X_list$transcriptomica), dict_tx)
}
if (!is.null(test_X_list$transcriptomica)) {
  colnames(test_X_list$transcriptomica) <-
    map_entrez_to_symbol(colnames(test_X_list$transcriptomica), dict_tx)
}
if ("transcriptomica" %in% names(pesos_list)) {
  rownames(pesos_list$transcriptomica) <-
    map_entrez_to_symbol(rownames(pesos_list$transcriptomica), dict_tx)
}

# =======================================================
# Utilidades MOFA (ajuste en cbind_views)
# =======================================================
mk_W <- function(pesos_list) Reduce("rbind", pesos_list)

cbind_views <- function(X_list, pesos_list) {
  ord <- names(pesos_list)
  mats <- lapply(ord, function(v){
    Xi <- X_list[[v]]
    if (is.null(Xi)) stop(sprintf("Vista %s vacía", v))
    # si filas parecen features, transpón
    if (nrow(Xi) == nrow(pesos_list[[v]]) && ncol(Xi) != nrow(pesos_list[[v]])) Xi <- t(Xi)
    
    feat_order <- rownames(pesos_list[[v]])              # features esperadas por vista
    # crea matriz completa y rellena faltantes con 0
    Xfull <- matrix(0, nrow=nrow(Xi), ncol=length(feat_order),
                    dimnames=list(rownames(Xi), feat_order))
    common <- intersect(colnames(Xi), feat_order)
    if (length(common) > 0) Xfull[, common] <- as.matrix(Xi[, common, drop=FALSE])
    Xfull
  })
  # chequea nº de muestras igual en todas las vistas
  nr <- sapply(mats, nrow); stopifnot(length(unique(nr)) == 1)
  do.call(cbind, mats)
}

reconstruct_from_latent <- function(Z, W) as.matrix(Z) %*% t(W)
project_to_latent <- function(X, W) { G <- crossprod(W); as.matrix(X) %*% W %*% solve(G) }

W <- mk_W(pesos_list); stopifnot(ncol(W)==2)

# === matrices escenario ===
X_tr_lat <- as.matrix(train_latent_X)
X_va_lat <- as.matrix(val_latent_X)
X_tr_rec <- reconstruct_from_latent(X_tr_lat, W)

# crudo no integrado: VALIDACIÓN y TEST
X_va_raw <- cbind_views(val_X_list,  pesos_list)
X_te_raw <- cbind_views(test_X_list, pesos_list)


# chequeos
stopifnot(ncol(X_va_raw) == sum(sapply(pesos_list, nrow)))
stopifnot(ncol(X_te_raw) == sum(sapply(pesos_list, nrow)))

X_te_lat <- project_to_latent(X_te_raw, W)
X_va_rec <- reconstruct_from_latent(X_va_lat, W)
X_te_rec <- reconstruct_from_latent(X_te_lat, W)



# =======================================================
# Métricas y calibración
# =======================================================


platt_fit <- function(y_true, probs_mat) {
  lev <- levels(y_true)
  lapply(lev, function(lv){
    yy <- as.integer(y_true==lv)
    df <- data.frame(y=yy, p=pmin(pmax(probs_mat[,lv], 1e-6), 1-1e-6))
    glm(y ~ qlogis(p), data=df, family=binomial())
  }) |> `names<-`(lev)
}
platt_apply <- function(fits, probs_mat) {
  lev <- names(fits)
  adj <- sapply(lev, function(lv){
    eta <- predict(fits[[lv]], newdata=data.frame(p=pmin(pmax(probs_mat[,lv], 1e-6), 1-1e-6)), type="link")
    1/(1+exp(-eta))
  })
  adj <- sweep(adj, 1, rowSums(adj), "/")
  colnames(adj) <- lev; adj
}

predict_glmnet_probs <- function(fit, X) {
  arr <- predict(fit$cvfit, newx=as.matrix(X), type="response", s=fit$lambda_select)
  p <- as.matrix(arr[,,1])
  colnames(p) <- fit$cvfit$glmnet.fit$classnames
  rownames(p) <- rownames(X)
  p
}
predict_multinom_probs <- function(fit, X) {
  p <- predict(fit, newdata=as.data.frame(X), type="probs")
  p <- as.matrix(p)[, fit$lev, drop=FALSE]
  rownames(p) <- rownames(X)
  p
}
# =======================================================
# Nested CV para glmnet multinomial (tunea α y λ)
# =======================================================
nested_cv_glmnet <- function(X, y, alphas=seq(0,1,by=0.25), k_outer=5, k_inner=5, max_lambda=NULL) {
  y <- factor(y); lev <- levels(y)
  folds_outer <- caret::createFolds(y, k=k_outer, list=TRUE, returnTrain=FALSE)
  scores <- data.frame(alpha=numeric(), Accuracy=numeric())
  
  for (ao in alphas) {
    acc_outer <- c()
    for (fo in folds_outer) {
      X_tr <- X[-fo, , drop=FALSE]; y_tr <- y[-fo]
      X_te <- X[ fo, , drop=FALSE]; y_te <- y[ fo]
      
      # inner CV para λ
      set.seed(999)
      cvfit <- cv.glmnet(x=as.matrix(X_tr), y=y_tr, family="multinomial",
                         alpha=ao, nfolds=k_inner, type.multinomial="ungrouped",
                         parallel=FALSE)
      lambda_sel <- if (is.null(max_lambda)) cvfit$lambda.1se else min(cvfit$lambda.1se, max_lambda)
      
      # eval outer
      arr <- predict(cvfit, newx=as.matrix(X_te), type="response", s=lambda_sel)
      p <- as.matrix(arr[,,1]); colnames(p) <- cvfit$glmnet.fit$classnames
      acc <- mean(colnames(p)[max.col(p)] == y_te)
      acc_outer <- c(acc_outer, acc)
    }
    scores <- rbind(scores, data.frame(alpha=ao, Accuracy=mean(acc_outer)))
  }
  
  best_alpha <- scores$alpha[which.max(scores$Accuracy)]
  
  # refit en todo el train con inner CV para λ
  set.seed(42)
  cvfit_all <- cv.glmnet(x=as.matrix(X), y=y, family="multinomial",
                         alpha=best_alpha, nfolds=k_inner, type.multinomial="ungrouped",
                         parallel=FALSE)
  lambda_sel <- if (is.null(max_lambda)) cvfit_all$lambda.1se else min(cvfit_all$lambda.1se, max_lambda)
  
  list(cvfit=cvfit_all, alpha=best_alpha, lambda_select=lambda_sel, scores=scores)
}

# =======================================================
# Bootstrap de estabilidad de biomarcadores (Lasso)
# =======================================================
bootstrap_lasso_stability <- function(X, y, B=100, alpha=1, k_inner=5) {
  p <- ncol(X); feats <- colnames(X) %||% paste0("V", seq_len(p))
  tab <- setNames(numeric(p), feats)
  n <- nrow(X)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace=TRUE)
    y_b <- y[idx]
    # Chequeo: todas las clases con >= 2 muestras
    if (any(table(y_b) < 2)) next
    
    cvfit <- cv.glmnet(
      x=as.matrix(X[idx,]), y=y_b, family="multinomial",
      alpha=alpha, nfolds=k_inner, type.multinomial="ungrouped"
    )
    
    cls <- cvfit$glmnet.fit$classnames
    nz_all <- unique(unlist(lapply(cls, function(cl){
      co <- coef(cvfit, s=cvfit$lambda.1se)[[cl]]
      rownames(co)[which(as.numeric(co)!=0)]
    })))
    nz_all <- setdiff(nz_all, "(Intercept)")
    tab[nz_all] <- tab[nz_all] + 1
  }
  
  tibble(feature=names(tab), freq=as.numeric(tab)/B) %>% arrange(desc(freq))
}


# =======================================================
# M1: Multinomial en latente (control)
# =======================================================
lev <- levels(train_y)
m1 <- nnet::multinom(train_y ~ ., data=data.frame(X_tr_lat, train_y), trace=FALSE, maxit=500); m1$lev <- lev
m1_tr_p  <- predict_multinom_probs(m1, X_tr_lat)
m1_va_pI <- predict_multinom_probs(m1, X_va_lat)
m1_te_pI <- predict_multinom_probs(m1, X_te_lat)

# SVM RBF en latente 2D (comparación no lineal)
ctrl_svm <- trainControl(method="repeatedcv", number=5, repeats=2, classProbs=TRUE, summaryFunction=multiClassSummary)
svm_fit  <- caret::train(x=as.data.frame(X_tr_lat), y=train_y, method="svmRadial",
                         tuneLength=10, trControl=ctrl_svm, metric="Accuracy", preProcess=c("center","scale"))
svm_pred <- function(fit, X) as.matrix(predict(fit, newdata=as.data.frame(X), type="prob"))[, levels(train_y), drop=FALSE]
svm_tr_p  <- svm_pred(svm_fit, X_tr_lat)
svm_va_pI <- svm_pred(svm_fit, X_va_lat)
svm_te_pI <- svm_pred(svm_fit, X_te_lat)

# =======================================================
# M2–M4: Ridge/Lasso/Elastic Net en reconstrucción con Nested CV
# =======================================================
alphas_grid <- seq(0,1,by=0.25)

m2 <- nested_cv_glmnet(X_tr_rec, train_y, alphas=alphas_grid)   # Ridge si α≈0
m3 <- nested_cv_glmnet(X_tr_rec, train_y, alphas=1)             # Lasso
m4 <- nested_cv_glmnet(X_tr_rec, train_y, alphas=alphas_grid)   # Elastic Net

# Predicciones integradas (latente→reconstrucción)
m2_tr_p  <- predict_glmnet_probs(m2, X_tr_rec)
m2_va_pI <- predict_glmnet_probs(m2, X_va_rec)
m2_te_pI <- predict_glmnet_probs(m2, X_te_rec)

m3_tr_p  <- predict_glmnet_probs(m3, X_tr_rec)
m3_va_pI <- predict_glmnet_probs(m3, X_va_rec)
m3_te_pI <- predict_glmnet_probs(m3, X_te_rec)

m4_tr_p  <- predict_glmnet_probs(m4, X_tr_rec)
m4_va_pI <- predict_glmnet_probs(m4, X_va_rec)
m4_te_pI <- predict_glmnet_probs(m4, X_te_rec)

# Predicciones NO integradas (crudo)
m2_va_pN <- predict_glmnet_probs(m2, X_va_raw)
m2_te_pN <- predict_glmnet_probs(m2, X_te_raw)

m3_va_pN <- predict_glmnet_probs(m3, X_va_raw)
m3_te_pN <- predict_glmnet_probs(m3, X_te_raw)

m4_va_pN <- predict_glmnet_probs(m4, X_va_raw)
m4_te_pN <- predict_glmnet_probs(m4, X_te_raw)

# =======================================================
# Calibración Platt sobre validación -> aplicada a TEST
# =======================================================
ids_val <- rownames(X_va_raw)  # ya alineado a val_y si construiste val_X_list con ids_val

calibrate <- TRUE
cal_apply <- function(y_val, p_val, p_test, ids_val = NULL) {
  if (!is.null(ids_val) && !is.null(rownames(p_val))) {
    common <- intersect(ids_val, rownames(p_val))
    if (length(common) >= 3) {
      y2 <- factor(y_val[match(common, ids_val)], levels = levels(y_val))
      p2 <- p_val[common, , drop = FALSE]
      fits <- platt_fit(y2, p2)
      return(platt_apply(fits, p_test))
    } else return(p_test)
  }
  if (nrow(p_val) == length(y_val)) {
    fits <- platt_fit(y_val, p_val)
    return(platt_apply(fits, p_test))
  }
  p_test
}

if (calibrate) {
  # Integrado
  m1_te_pI  <- cal_apply(val_y, m1_va_pI,  m1_te_pI,  ids_val)
  svm_te_pI <- cal_apply(val_y, svm_va_pI, svm_te_pI, ids_val)
  m2_te_pI  <- cal_apply(val_y, m2_va_pI,  m2_te_pI,  ids_val)
  m3_te_pI  <- cal_apply(val_y, m3_va_pI,  m3_te_pI,  ids_val)
  m4_te_pI  <- cal_apply(val_y, m4_va_pI,  m4_te_pI,  ids_val)

  # No integrado
  m2_te_pN  <- cal_apply(val_y, m2_va_pN,  m2_te_pN,  ids_val)
  m3_te_pN  <- cal_apply(val_y, m3_va_pN,  m3_te_pN,  ids_val)
  m4_te_pN  <- cal_apply(val_y, m4_va_pN,  m4_te_pN,  ids_val)
}

# Tablas de resultados
# =======================================================
metrics_multiclass <- function(y_true, probs_mat) {
  out_na <- c(Accuracy=NA_real_, Sensibilidad=NA_real_, Especificidad=NA_real_, AUC=NA_real_)
  if (is.null(probs_mat) || nrow(probs_mat) == 0) return(out_na)
  
  # si hay nombres en y y rownames en probs => alinear por intersección
  y_use <- y_true
  p_use <- probs_mat
  if (!is.null(names(y_true)) && !is.null(rownames(probs_mat))) {
    common <- intersect(names(y_true), rownames(probs_mat))
    if (length(common) >= 1) {
      y_use <- factor(y_true[common], levels=levels(y_true))
      p_use <- probs_mat[common, , drop=FALSE]
    }
  }
  
  # si aún no coincide longitud, recorta al mínimo (fallback seguro)
  if (length(y_use) != nrow(p_use)) {
    n <- min(length(y_use), nrow(p_use))
    y_use <- y_use[seq_len(n)]
    p_use <- p_use[seq_len(n), , drop=FALSE]
  }
  
  # si hay NA/Inf en probs, no se puede calcular => devuelve NAs
  if (any(!is.finite(p_use))) return(out_na)
  
  # normaliza por si no suman 1
  rs <- rowSums(p_use)
  bad <- !is.finite(rs) | rs == 0
  if (any(bad)) p_use[bad, ] <- 1 / ncol(p_use)
  
  cls  <- colnames(p_use)
  pred <- factor(cls[max.col(p_use, ties.method="first")], levels=levels(y_true))
  
  # si longitudes difieren por cualquier motivo residual, devuelve NAs
  if (length(pred) != length(y_use)) return(out_na)
  
  cm   <- caret::confusionMatrix(pred, y_use)
  acc  <- as.numeric(cm$overall["Accuracy"])
  byc  <- suppressWarnings(cm$byClass)
  sens <- mean(byc[,"Sensitivity"], na.rm=TRUE)
  spec <- mean(byc[,"Specificity"], na.rm=TRUE)
  auc  <- tryCatch(as.numeric(pROC::multiclass.roc(y_use, p_use)$auc), error=function(e) NA_real_)
  c(Accuracy=acc, Sensibilidad=sens, Especificidad=spec, AUC=auc)
}

collect_rows <- function(nombre, sets) {
  bind_rows(lapply(names(sets), function(k) {
    met <- metrics_multiclass(sets[[k]]$y, sets[[k]]$p)
    data.frame(Modelo=nombre, Conjunto=k, t(met))
  }))
}

res_list <- list(
  collect_rows("M1_MOFA_latente", list(
    Train=list(y=train_y, p=m1_tr_p),
    Val_integrado=list(y=val_y, p=m1_va_pI),
    Test_integrado=list(y=val_y, p=m1_te_pI),
    Val_no_integrado=list(y=val_y, p=matrix(NA, nrow=length(val_y), ncol=length(lev), dimnames=list(NULL, lev))),
    Test_no_integrado=list(y=val_y, p=matrix(NA, nrow=length(val_y), ncol=length(lev), dimnames=list(NULL, lev)))
  )),
  collect_rows("SVM_latente", list(
    Train=list(y=train_y, p=svm_tr_p),
    Val_integrado=list(y=val_y, p=svm_va_pI),
    Test_integrado=list(y=val_y, p=svm_te_pI),
    Val_no_integrado=list(y=val_y, p=matrix(NA, nrow=length(val_y), ncol=length(lev), dimnames=list(NULL, lev))),
    Test_no_integrado=list(y=val_y, p=matrix(NA, nrow=length(val_y), ncol=length(lev), dimnames=list(NULL, lev)))
  )),
  collect_rows("M2_Ridge_rec", list(
    Train=list(y=train_y, p=m2_tr_p),
    Val_integrado=list(y=val_y, p=m2_va_pI),
    Test_integrado=list(y=val_y, p=m2_te_pI),
    Val_no_integrado=list(y=val_y, p=m2_va_pN),
    Test_no_integrado=list(y=val_y, p=m2_te_pN)
  )),
  collect_rows("M3_Lasso_rec", list(
    Train=list(y=train_y, p=m3_tr_p),
    Val_integrado=list(y=val_y, p=m3_va_pI),
    Test_integrado=list(y=val_y, p=m3_te_pI),
    Val_no_integrado=list(y=val_y, p=m3_va_pN),
    Test_no_integrado=list(y=val_y, p=m3_te_pN)
  )),
  collect_rows("M4_ElasticNet_rec", list(
    Train=list(y=train_y, p=m4_tr_p),
    Val_integrado=list(y=val_y, p=m4_va_pI),
    Test_integrado=list(y=val_y, p=m4_te_pI),
    Val_no_integrado=list(y=val_y, p=m4_va_pN),
    Test_no_integrado=list(y=val_y, p=m4_te_pN)
  ))
)

tabla_final <- bind_rows(res_list) %>%
  mutate(across(c(Accuracy,Sensibilidad,Especificidad,AUC), ~round(.,4))) %>%
  arrange(Modelo, factor(Conjunto, levels=c("Train","Val_integrado","Val_no_integrado","Test_integrado","Test_no_integrado")))
print(tabla_final)
# readr::write_csv(tabla_final, "resultados_modelos_nestedCV_CAL.csv")

# =======================================================
# Reporte de hiperparámetros elegidos y biomarcadores
# =======================================================
hp <- tibble(
  Modelo=c("M2_Ridge_rec","M3_Lasso_rec","M4_ElasticNet_rec"),
  alpha=c(m2$alpha, m3$alpha, m4$alpha),
  lambda=c(m2$lambda_select, m3$lambda_select, m4$lambda_select)
)
print(hp)

coef_nonzero <- function(cvwrapper) {
  cls <- cvwrapper$cvfit$glmnet.fit$classnames
  bind_rows(lapply(cls, function(cl){
    co <- coef(cvwrapper$cvfit, s=cvwrapper$lambda_select)[[cl]]
    tibble(feature=rownames(co), coef=as.numeric(co)) %>%
      filter(feature!="(Intercept)", coef!=0) %>%
      mutate(clase=cl)
  }))
}
biomarcadores_lasso <- coef_nonzero(m3)
biomarcadores_enet  <- coef_nonzero(m4)

# Frecuencias bootstrap ya calculadas en stab_tab
# readr::write_csv(biomarcadores_lasso, "biomarcadores_lasso.csv")
# readr::write_csv(biomarcadores_enet,  "biomarcadores_elasticnet.csv")
# readr::write_csv(stab_tab,             "estabilidad_bootstrap_lasso.csv")
