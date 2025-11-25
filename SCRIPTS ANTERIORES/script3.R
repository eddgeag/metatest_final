# =======================================================
# Pipeline Z=2 (MOFA) con síntesis mejorada, diagnóstico,
# visualización y modelos (multinom control + lasso)
# =======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet)
  library(nnet)
  library(MASS)
  library(caret)
  library(pROC)
  library(FNN)
  library(ggplot2)
  library(ggridges)
})

set.seed(123)

# =======================================================
# Utilidades de métricas
# =======================================================
metrics_multiclass2 <- function(y_true, probs_mat) {
  out_na <- c(Accuracy=NA, BalancedAcc=NA, MacroF1=NA, Top2Acc=NA, AUC_HT=NA, AUC_1vsAll=NA)
  if (is.null(probs_mat) || nrow(probs_mat) == 0) return(out_na)
  y_use <- factor(y_true, levels=levels(y_true))
  P <- as.matrix(probs_mat)[, levels(y_use), drop=FALSE]
  pred <- factor(colnames(P)[max.col(P, ties.method = "first")], levels = levels(y_use))
  cm   <- caret::confusionMatrix(pred, y_use)
  acc  <- as.numeric(cm$overall["Accuracy"])
  bacc <- mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE)
  f1s  <- 2*cm$byClass[,"Sensitivity"]*cm$byClass[,"Pos Pred Value"] /
    (cm$byClass[,"Sensitivity"]+cm$byClass[,"Pos Pred Value"])
  macroF1 <- mean(f1s, na.rm=TRUE)
  # top-2
  # dentro de metrics_multiclass2, reemplaza top-2 por:
  # dentro de metrics_multiclass2, reemplaza top-2 por:
  top2 <- mean(sapply(seq_len(nrow(P)), function(i){
    idx_true <- which(colnames(P) == as.character(y_use[i]))
    idx_true %in% order(P[i,], decreasing=TRUE)[1:2]
  }))
  
  # AUCs
  auc_ht <- tryCatch(as.numeric(pROC::multiclass.roc(y_use, P)$auc), error=function(e) NA)
  auc_ova <- tryCatch({
    levs <- levels(y_use)
    mean(sapply(levs, function(lv){
      yb <- factor(y_use==lv, levels=c(FALSE, TRUE))
      as.numeric(pROC::auc(pROC::roc(yb, P[,lv], quiet=TRUE)))
    }), na.rm=TRUE)
  }, error=function(e) NA)
  c(Accuracy=acc, BalancedAcc=bacc, MacroF1=macroF1, Top2Acc=top2, AUC_HT=auc_ht, AUC_1vsAll=auc_ova)
}

# =======================================================
# Calibración Platt isotónica ya provista por el usuario
# =======================================================
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
  Z <- sweep(Z, 1, rowSums(Z), "/"); Z
}

# =======================================================
# Inputs (MOFA Z=2, grupos y pesos W para proyectar TEST)
# =======================================================
dats      <- readRDS("./datos_para_modelar.rds")
tst_Data  <- readRDS("./ready_for_modeling.rds")

factors    <- dats$factores
pesos_list <- dats$pesos

Z_seed <- as.matrix(factors %>% filter(Tipo=="Real") %>% dplyr::select(where(is.numeric)))[,1:2,drop=FALSE]
y_seed <- factors %>% filter(Tipo=="Real") %>% pull(Grupo) %>% factor()
names(y_seed) <- rownames(Z_seed)
lev <- levels(y_seed)

# Split anchors/holdout
split_seed <- function(y, k_anchor = 2L){
  L <- levels(y); ids_all <- names(y)
  if (length(k_anchor) > 1L) {
    stopifnot(!is.null(names(k_anchor)), all(names(k_anchor) %in% L))
    k_vec <- as.integer(k_anchor[L])
  } else k_vec <- rep(as.integer(k_anchor), length(L))
  out <- unlist(mapply(function(cl,k){
    ids_cl <- ids_all[y==cl]
    if (!length(ids_cl)) return(character(0))
    k <- min(max(0L,k), max(0L,length(ids_cl)-1L))
    if (k==0) return(character(0))
    sample(ids_cl, k, replace=FALSE)
  }, L, k_vec, SIMPLIFY=FALSE), use.names=FALSE)
  unique(out)
}

tbl <- table(y_seed); n_i <- as.integer(tbl)
desired   <- pmax(5L, floor(0.2 * n_i))
k_anchor  <- setNames(pmax(1L, pmin(desired, pmax(0L, n_i - 1L))), names(tbl))
ids_anchor <- split_seed(y_seed, k_anchor = k_anchor)
ids_val    <- setdiff(names(y_seed), ids_anchor)

Z_fit <- Z_seed[ids_anchor, , drop=FALSE]
y_fit <- y_seed[ids_anchor]
Z_val <- Z_seed[ids_val, , drop=FALSE]
y_val <- y_seed[ids_val]

# =======================================================
# Síntesis mejorada en Z=2: Gauss truncada + SMOTE + alineación de momentos
# =======================================================
smote_Z <- function(Zc, n, k=5, noise_sd=0.05){
  if (n<=0 || nrow(Zc)<2) return(NULL)
  k <- min(k, max(1, nrow(Zc)-1))
  nn <- FNN::get.knn(Zc, k=k)$nn.index
  i <- sample(1:nrow(Zc), n, replace=TRUE)
  j <- sapply(i, function(ii) sample(nn[ii,], 1))
  t <- runif(n)
  X <- Zc[i,,drop=FALSE] + t*(Zc[j,,drop=FALSE] - Zc[i,,drop=FALSE])
  X + matrix(rnorm(n*ncol(Zc), 0, noise_sd), n, ncol(Zc))
}
# reemplaza align_moments y gen_class_mix por versiones con fallback pooled

# reemplaza align_moments y gen_class_mix por versiones con fallback pooled

align_moments <- function(Zsyn, Zreal, Zpool=NULL){
  mu_s <- colMeans(Zsyn); mu_r <- colMeans(Zreal)
  # elegir S_r robusta
  S_r <- if (nrow(Zreal) >= 3) cov(Zreal) else {
    if (!is.null(Zpool) && nrow(Zpool) >= 3) cov(Zpool) else diag(apply(Zsyn,2,sd)^2)
  }
  # S_s desde sintéticos
  S_s <- cov(Zsyn)
  Es <- eigen(S_s, symmetric=TRUE); Er <- eigen(S_r, symmetric=TRUE)
  W  <- Es$vectors %*% diag(1/sqrt(pmax(Es$values, 1e-6))) %*% t(Es$vectors)
  C  <- Er$vectors %*% diag( sqrt(pmax(Er$values, 1e-6))) %*% t(Er$vectors)
  Zc <- sweep(Zsyn, 2, mu_s, "-") %*% W %*% C
  sweep(Zc, 2, mu_r, "+")
}

gen_class_mix <- function(Zc, n_gauss, n_smote, shrink=0.3, rad_q=0.95, Zpool=NULL){
  p <- ncol(Zc); mu <- colMeans(Zc)
  # S de referencia: pooled si pocos reales
  Zref <- if (nrow(Zc) >= 3) Zc else if(!is.null(Zpool)) Zpool else Zc
  Sref <- cov(Zref); Sref <- (1-shrink)*Sref + shrink*diag(diag(Sref))
  Xg <- if(n_gauss>0) MASS::mvrnorm(n_gauss, mu, Sref) else NULL
  if (!is.null(Xg) && nrow(Zref) >= 3){
    D2  <- mahalanobis(Xg, center=colMeans(Zref), cov=Sref)
    thr <- quantile(mahalanobis(Zref, colMeans(Zref), Sref), rad_q)
    idx <- which(D2 > thr)
    if (length(idx)) Xg[idx,] <- sweep(Xg[idx,], 1, sqrt(thr / D2[idx]), `*`)
  }
  Xs <- if(n_smote>0 && nrow(Zc)>=2) smote_Z(Zc, n_smote) else NULL
  X  <- rbind(Xg, Xs)
  align_moments(X, Zc, Zpool=Zpool)
}

synth_Z_mix <- function(Z_fit, y_fit, targets, shrink=0.3, rad_q=0.95){
  lev <- levels(y_fit); Zs <- list(); ys <- c()
  for(lv in lev){
    Zc <- Z_fit[y_fit==lv,,drop=FALSE]
    need <- targets[[lv]] - nrow(Zc)
    if (need > 0){
      Xnew <- gen_class_mix(Zc, ceiling(0.5*need), floor(0.5*need),
                            shrink=shrink, rad_q=rad_q, Zpool=Z_fit)
      Zs[[lv]] <- Xnew; ys <- c(ys, rep(lv, nrow(Xnew)))
    }
  }
  Zs <- if(length(Zs)) do.call(rbind,Zs) else NULL
  list(Z=Zs, y=factor(ys, levels=lev))
}

# Targets más altos para cubrir variabilidad
k <- 5L; min_fold <- 8L
min_total  <- ceiling(min_fold * k / (k - 1))
max_ratio  <- 3L; base_target <- 15L
raw_target <- pmin(base_target, as.integer(tbl) * max_ratio)
# mínimo 40 por clase
targets    <- setNames(pmax(raw_target, 40L), names(tbl))

syn <- synth_Z_mix(Z_fit, y_fit, targets, shrink=0.3, rad_q=0.95)
stopifnot(!is.null(syn$Z))

train_latent_X <- rbind(Z_fit, syn$Z)
train_y        <- factor(c(as.character(y_fit), as.character(syn$y)), levels=lev)

# pesos: reales=1, sintéticos=0.4
w_train <- c(rep(1, nrow(Z_fit)), rep(0.4, nrow(syn$Z)))

# =======================================================
# Diagnóstico rápido de sintéticos vs reales (por clase) en Z
# =======================================================
# Distancia de Bhattacharyya entre Gaussianas por clase + diferencias de medias/covs
bhatt_class <- function(Zr, Zs){
  mu1 <- colMeans(Zr); S1 <- cov(Zr)
  mu2 <- colMeans(Zs); S2 <- cov(Zs)
  S   <- 0.5*(S1+S2)
  invS <- tryCatch(solve(S), error=function(e) MASS::ginv(S))
  term1 <- 0.125 * t(mu1-mu2) %*% invS %*% (mu1-mu2)
  term2 <- 0.5 * log( det(S) / sqrt(pmax(det(S1),1e-12)*pmax(det(S2),1e-12)) )
  as.numeric(term1 + term2)
}

diag_sinteticos <- function(Z_fit, y_fit, Z_syn, y_syn, S_pool=NULL){
  stopifnot(!is.null(S_pool))
  lev <- levels(y_fit)
  out <- lapply(lev, function(lv){
    Zr <- Z_fit[y_fit==lv,,drop=FALSE]; Zs <- Z_syn[y_syn==lv,,drop=FALSE]
    if (nrow(Zs)<2) return(NULL)
    mu_diff <- sqrt(sum((colMeans(Zr)-colMeans(Zs))^2))
    # Sigma de referencia robusta:
    S_ref <- if (nrow(Zr)>=3) cov(Zr) else S_pool
    db <- bhatt_class( # usa S_ref para Zr
      Zr = if (nrow(Zr)>=2) Zr else sweep(MASS::mvrnorm(50, colMeans(Z_syn[y_syn==lv,,drop=FALSE]), S_ref),2,0,"+"),
      Zs = Zs
    )
    data.frame(clase=lv, n_real=nrow(Zr), n_syn=nrow(Zs),
               mu_L2=round(mu_diff,4),
               cov_F=round(norm(S_ref - cov(Zs), "F"),4),
               Bhattacharyya=round(db,4))
  })
  do.call(rbind, out)
}
# úsalo así:
S_pool <- cov(Z_fit)
diag_tbl <- diag_sinteticos(Z_fit, y_fit, syn$Z, syn$y, S_pool=S_pool)


# =======================================================
# Visualización en Z: función reusable
# =======================================================
plot_z_real_vs_syn <- function(Z_fit, y_fit, Z_syn, y_syn){
  df_r <- data.frame(Z1=Z_fit[,1], Z2=Z_fit[,2], Grupo=y_fit, Set="Real_train")
  df_s <- data.frame(Z1=Z_syn[,1], Z2=Z_syn[,2], Grupo=y_syn, Set="Sintetico")
  df   <- rbind(df_r, df_s)
  p1 <- ggplot(df, aes(Z1, Z2, color=Grupo, shape=Set)) +
    geom_point(alpha=0.7, size=2.8) +
    stat_ellipse(aes(group=interaction(Grupo,Set)), linetype=2, alpha=0.25) +
    theme_minimal(base_size=14) +
    labs(title="Z=2 MOFA: Reales vs Sintéticos", subtitle="Puntos y elipses por Grupo/Set")
  p2 <- ggplot(df, aes(Z1, fill=Set, color=Set)) +
    geom_density(alpha=0.3) + facet_wrap(~Grupo, scales="free") +
    theme_minimal(base_size=14) + labs(title="Distribución de Z1 por grupo")
  p3 <- ggplot(df, aes(Z2, fill=Set, color=Set)) +
    geom_density(alpha=0.3) + facet_wrap(~Grupo, scales="free") +
    theme_minimal(base_size=14) + labs(title="Distribución de Z2 por grupo")
  print(p1); print(p2); print(p3)
}

# Llamada de visualización
plot_z_real_vs_syn(Z_fit, y_fit, syn$Z, syn$y)

# =======================================================
# Proyección y TEST en latentes (usando W y bloques de test)
# =======================================================
W <- do.call(rbind, pesos_list)
stopifnot(!is.null(rownames(W)))
feat_order <- rownames(W)

X_test_blocks <- list(
  transcriptomica = tst_Data$tx$test,
  proteomica     = tst_Data$pr$test,
  metabolomica   = tst_Data$me$test,
  clinica        = tst_Data$cl$test
)
X_test_all <- do.call(rbind, X_test_blocks)
X_test_all <- X_test_all[feat_order, , drop=FALSE]
X_test_norec <- t(X_test_all)

project_to_latent <- function(X, W){ G <- crossprod(W); X %*% W %*% MASS::ginv(G) }
X_test_lat <- project_to_latent(X_test_norec, W)

y_test <- tst_Data$grupo_test$grupo

train_y <- as.factor(as.character(train_y))
y_test <- as.factor(as.character(y_test))
stopifnot(all(levels(train_y) == levels(y_test)))

# =======================================================
# Modelos en Z: Multinom (control) y Lasso (glmnet)
# =======================================================
# Multinomial sin regularización
set.seed(123)
df_tr_lat <- data.frame(train_latent_X, y=train_y)
mn_fit <- nnet::multinom(y ~ ., data=df_tr_lat, trace=FALSE, maxit=1000)
colnames(X_test_lat) <- c("Factor1","Factor2")
probs_mn <- predict(mn_fit, newdata=as.data.frame(X_test_lat), type="probs")

# Lasso multinomial con pesos
set.seed(123)
cv_glm <- cv.glmnet(as.matrix(train_latent_X), train_y, family="multinomial",
                    alpha=1, type.measure="class", weights=w_train, nfolds=5)
probs_lasso <- predict(cv_glm, newx=as.matrix(X_test_lat), s="lambda.1se", type="response")[,,1]

# =======================================================
# Calibración con holdout real (Platt isotónica) para Lasso
# =======================================================
P_val <- predict(cv_glm, newx=as.matrix(Z_val), s="lambda.1se", type="response")[,,1]
cal    <- platt_fit(y_val, as.matrix(P_val))
P_te   <- probs_lasso
P_te_c <- platt_apply(cal, as.matrix(P_te))

# =======================================================
# Reporte de métricas y matrices de confusión
# =======================================================
mm_mn     <- metrics_multiclass2(y_test, probs_mat = probs_mn)
mm_lasso  <- metrics_multiclass2(y_test, probs_mat = probs_lasso)
mm_lassoC <- metrics_multiclass2(y_test, probs_mat = P_te_c)

print("== Métricas Multinom (Z) ==");  print(mm_mn)
print("== Métricas Lasso (Z) ==");     print(mm_lasso)
print("== Métricas Lasso Calibrado (Z) =="); print(mm_lassoC)

# Confusion matrices
pred_mn     <- factor(colnames(probs_mn)[max.col(probs_mn)], levels=levels(y_test))
pred_lasso  <- factor(colnames(probs_lasso)[max.col(probs_lasso)], levels=levels(y_test))
pred_lassoC <- factor(colnames(P_te_c)[max.col(P_te_c)], levels=levels(y_test))

cm_mn     <- caret::confusionMatrix(pred_mn,     y_test)
cm_lasso  <- caret::confusionMatrix(pred_lasso,  y_test)
cm_lassoC <- caret::confusionMatrix(pred_lassoC, y_test)

print("== CM Multinom (Z) ==");  print(cm_mn)
print("== CM Lasso (Z) ==");     print(cm_lasso)
print("== CM Lasso Calibrado (Z) =="); print(cm_lassoC)


# ===== Toggle de calibración =====
CALIBRAR <- FALSE  # <- déjalo en FALSE

# Bloque actual
P_val <- predict(cv_glm, newx=as.matrix(Z_val), s="lambda.1se", type="response")[,,1]
cal    <- platt_fit(y_val, as.matrix(P_val))
P_te   <- probs_lasso
if (CALIBRAR) {
  P_te_c <- platt_apply(cal, as.matrix(P_te))
} else {
  P_te_c <- probs_lasso  # sin calibrar
}

# Métricas y CM (dejan de “castigar” por la mala calibración)
mm_lassoC <- metrics_multiclass2(y_test, probs_mat = P_te_c)
pred_lassoC <- factor(colnames(P_te_c)[max.col(P_te_c)], levels=levels(y_test))
cm_lassoC <- caret::confusionMatrix(pred_lassoC, y_test)

###====aqui


# =======================================================
# Pipeline Z=2 (MOFA) con síntesis mejorada (mezclas),
# diagnóstico, visualización y modelos
# =======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet)
  library(nnet)
  library(MASS)
  library(caret)
  library(pROC)
  library(FNN)
  library(ggplot2)
  library(ggridges)
})

set.seed(123)

# =======================================================
# Utilidades de métricas
# =======================================================
metrics_multiclass2 <- function(y_true, probs_mat) {
  out_na <- c(Accuracy=NA, BalancedAcc=NA, MacroF1=NA, Top2Acc=NA, AUC_HT=NA, AUC_1vsAll=NA)
  if (is.null(probs_mat) || nrow(probs_mat) == 0) return(out_na)
  y_use <- factor(y_true, levels=levels(y_true))
  P <- as.matrix(probs_mat)[, levels(y_use), drop=FALSE]
  pred <- factor(colnames(P)[max.col(P, ties.method = "first")], levels = levels(y_use))
  cm   <- caret::confusionMatrix(pred, y_use)
  acc  <- as.numeric(cm$overall["Accuracy"])
  bacc <- mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE)
  f1s  <- 2*cm$byClass[,"Sensitivity"]*cm$byClass[,"Pos Pred Value"] /
    (cm$byClass[,"Sensitivity"]+cm$byClass[,"Pos Pred Value"])
  macroF1 <- mean(f1s, na.rm=TRUE)
  top2 <- mean(sapply(seq_len(nrow(P)), function(i){
    idx_true <- which(colnames(P) == as.character(y_use[i]))
    idx_true %in% order(P[i,], decreasing=TRUE)[1:2]
  }))
  auc_ht <- tryCatch(as.numeric(pROC::multiclass.roc(y_use, P)$auc), error=function(e) NA)
  auc_ova <- tryCatch({
    levs <- levels(y_use)
    mean(sapply(levs, function(lv){
      yb <- factor(y_use==lv, levels=c(FALSE, TRUE))
      as.numeric(pROC::auc(pROC::roc(yb, P[,lv], quiet=TRUE)))
    }), na.rm=TRUE)
  }, error=function(e) NA)
  c(Accuracy=acc, BalancedAcc=bacc, MacroF1=macroF1, Top2Acc=top2, AUC_HT=auc_ht, AUC_1vsAll=auc_ova)
}

# =======================================================
# Calibración Platt isotónica (toggle configurable)
# =======================================================
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
  sweep(Z, 1, rowSums(Z), "/")
}

# =======================================================
# Inputs
# =======================================================
dats      <- readRDS("./datos_para_modelar.rds")
tst_Data  <- readRDS("./ready_for_modeling.rds")
factors    <- dats$factores
pesos_list <- dats$pesos

Z_seed <- as.matrix(factors %>% filter(Tipo=="Real") %>% dplyr::select(where(is.numeric)))[,1:2,drop=FALSE]
y_seed <- factors %>% filter(Tipo=="Real") %>% pull(Grupo) %>% factor()
names(y_seed) <- rownames(Z_seed)
lev <- levels(y_seed)

# split anchors/holdout
split_seed <- function(y, k_anchor = 2L){
  L <- levels(y); ids_all <- names(y)
  if (length(k_anchor) > 1L) {
    stopifnot(!is.null(names(k_anchor)), all(names(k_anchor) %in% L))
    k_vec <- as.integer(k_anchor[L])
  } else k_vec <- rep(as.integer(k_anchor), length(L))
  out <- unlist(mapply(function(cl,k){
    ids_cl <- ids_all[y==cl]
    if (!length(ids_cl)) return(character(0))
    k <- min(max(0L,k), max(0L,length(ids_cl)-1L))
    if (k==0) return(character(0))
    sample(ids_cl, k, replace=FALSE)
  }, L, k_vec, SIMPLIFY=FALSE), use.names=FALSE)
  unique(out)
}
tbl <- table(y_seed)
desired   <- pmax(5L, floor(0.2 * as.integer(tbl)))
k_anchor  <- setNames(pmax(1L, pmin(desired, pmax(0L, as.integer(tbl) - 1L))), names(tbl))
ids_anchor <- split_seed(y_seed, k_anchor = k_anchor)
ids_val    <- setdiff(names(y_seed), ids_anchor)

Z_fit <- Z_seed[ids_anchor, , drop=FALSE]
y_fit <- y_seed[ids_anchor]
Z_val <- Z_seed[ids_val, , drop=FALSE]
y_val <- y_seed[ids_val]

# =======================================================
# Síntesis en Z=2 con mezclas por clase + control varianza
# =======================================================
smote_Z <- function(Zc, n, k=5, noise_sd=0.05){
  if (n<=0 || nrow(Zc)<2) return(NULL)
  k <- min(k, max(1, nrow(Zc)-1))
  nn <- FNN::get.knn(Zc, k=k)$nn.index
  i <- sample(1:nrow(Zc), n, replace=TRUE)
  j <- sapply(i, function(ii) sample(nn[ii,], 1))
  t <- runif(n)
  X <- Zc[i,,drop=FALSE] + t*(Zc[j,,drop=FALSE] - Zc[i,,drop=FALSE])
  X + matrix(rnorm(n*ncol(Zc), 0, noise_sd), n, ncol(Zc))
}
align_moments <- function(Zsyn, Zreal, Zpool=NULL){
  mu_s <- colMeans(Zsyn); mu_r <- colMeans(Zreal)
  S_r <- if (nrow(Zreal) >= 3) cov(Zreal) else {
    if (!is.null(Zpool) && nrow(Zpool) >= 3) cov(Zpool) else diag(apply(Zsyn,2,sd)^2)
  }
  S_s <- cov(Zsyn)
  Es <- eigen(S_s, symmetric=TRUE); Er <- eigen(S_r, symmetric=TRUE)
  W  <- Es$vectors %*% diag(1/sqrt(pmax(Es$values, 1e-6))) %*% t(Es$vectors)
  C  <- Er$vectors %*% diag( sqrt(pmax(Er$values, 1e-6))) %*% t(Er$vectors)
  Zc <- sweep(Zsyn, 2, mu_s, "-") %*% W %*% C
  sweep(Zc, 2, mu_r, "+")
}
rescale_variance_pc <- function(Zsyn, Zreal, min_sd=0.15, max_sd=Inf){
  Er <- eigen(cov(Zreal), symmetric=TRUE)
  U <- Er$vectors
  Xr <- scale(Zreal, center=colMeans(Zreal), scale=FALSE) %*% U
  Xs <- scale(Zsyn,  center=colMeans(Zsyn),  scale=FALSE) %*% U
  sd_r <- pmax(apply(Xr,2,sd), 1e-6)
  sd_t <- pmin(pmax(sd_r, min_sd), max_sd)
  Xs2  <- sweep(Xs, 2, apply(Xs,2,sd)+1e-6, "/") %*% diag(sd_t)
  sweep(Xs2 %*% t(U), 2, colMeans(Zreal), "+")
}
# util: rbind seguro con NULL
rb <- function(a,b){
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  rbind(a,b)
}

gen_class_mix <- function(Zc, n_gauss, n_smote, shrink=0.3, rad_q=0.95, Zpool=NULL){
  mu <- colMeans(Zc)
  
  # cov de referencia (pooled si pocos reales)
  Zref <- if (nrow(Zc) >= 3) Zc else if(!is.null(Zpool)) Zpool else Zc
  Sref <- cov(Zref); Sref <- (1-shrink)*Sref + shrink*diag(diag(Sref))
  
  # parte gauss
  Xg <- if (n_gauss>0) MASS::mvrnorm(n_gauss, mu, Sref) else NULL
  
  # truncado por Mahalanobis (robusto a casos borde)
  if (!is.null(Xg) && nrow(Zref) >= 3){
    D2  <- mahalanobis(Xg, center=colMeans(Zref), cov=Sref)
    Dz  <- mahalanobis(Zref, center=colMeans(Zref), cov=Sref)
    thr <- as.numeric(stats::quantile(Dz, probs=rad_q, names=FALSE))
    idx <- which(is.finite(D2) & D2 > thr)
    if (length(idx) > 0){
      s <- sqrt(thr / D2[idx])
      # re-escalado fila-a-fila, evitando drop
      Xg[idx, ] <- Xg[idx, , drop=FALSE] * matrix(s, nrow=length(idx), ncol=ncol(Xg))
    }
  }
  
  # parte SMOTE (si hay al menos 2 reales)
  Xs <- if (n_smote>0 && nrow(Zc)>=2) smote_Z(Zc, n_smote) else NULL
  
  # combinar y alinear
  X  <- rb(Xg, Xs)
  if (is.null(X)) {
    # fallback: si por alguna razón no hay X, replica leves perturbaciones alrededor de mu
    X <- matrix(rep(mu, each=max(1, n_gauss+n_smote)), ncol=ncol(Zc), byrow=FALSE) +
      matrix(rnorm(max(1, n_gauss+n_smote)*ncol(Zc), 0, 0.05), ncol=ncol(Zc))
  }
  
  X  <- align_moments(X, Zc, Zpool=Zpool)
  if (nrow(Zc) >= 3) X <- rescale_variance_pc(X, Zc, min_sd=0.12)
  X
}
synth_Z_mix_kmeans <- function(Z_fit, y_fit, targets, k_max=3, shrink=0.3, rad_q=0.95){
  lev <- levels(y_fit); Zs <- list(); ys <- c()
  for(lv in lev){
    Zc <- Z_fit[y_fit==lv,,drop=FALSE]; nr <- nrow(Zc)
    need <- targets[[lv]] - nr
    if (need <= 0) next
    if (nr >= 2) {
      k <- min(k_max, nr)
      if (k > 1 && nr > k) {
        set.seed(1)
        km <- kmeans(Zc, centers=k, nstart=20)
        alloc <- round(need * table(km$cluster)/nr)
        alloc[which.max(alloc)] <- alloc[which.max(alloc)] + (need - sum(alloc))
        for(cc in seq_len(k)){
          Zcc <- Zc[km$cluster==cc,,drop=FALSE]
          n_new <- alloc[as.character(cc)]
          if (n_new>0){
            Xnew <- gen_class_mix(Zcc, ceiling(0.5*n_new), floor(0.5*n_new),
                                  shrink=shrink, rad_q=rad_q, Zpool=Z_fit)
            Zs[[length(Zs)+1]] <- Xnew
            ys <- c(ys, rep(lv, nrow(Xnew)))
          }
        }
      } else {
        # nr == k o k==1 → sin kmeans, generar directo
        Xnew <- gen_class_mix(Zc, ceiling(0.5*need), floor(0.5*need),
                              shrink=shrink, rad_q=rad_q, Zpool=Z_fit)
        Zs[[length(Zs)+1]] <- Xnew
        ys <- c(ys, rep(lv, nrow(Xnew)))
      }
    } else {
      # nr==1 → fallback
      Xnew <- gen_class_mix(Zc, ceiling(0.5*need), floor(0.5*need),
                            shrink=shrink, rad_q=rad_q, Zpool=Z_fit)
      Zs[[length(Zs)+1]] <- Xnew
      ys <- c(ys, rep(lv, nrow(Xnew)))
    }
    
  }
  Zs <- if(length(Zs)) do.call(rbind,Zs) else NULL
  list(Z=Zs, y=factor(ys, levels=lev))
}
# targets
k <- 5L; min_fold <- 8L
min_total  <- ceiling(min_fold * k / (k - 1))
max_ratio  <- 3L; base_target <- 15L
raw_target <- pmin(base_target, as.integer(tbl) * max_ratio)
targets    <- setNames(pmax(raw_target, 40L), names(tbl))
syn <- synth_Z_mix_kmeans(Z_fit, y_fit, targets, k_max=3, shrink=0.3, rad_q=0.95)
train_latent_X <- rbind(Z_fit, syn$Z)
train_y        <- factor(c(as.character(y_fit), as.character(syn$y)), levels=lev)
w_train <- c(rep(1, nrow(Z_fit)),
             ifelse(as.character(syn$y)=="FD", 0.3, 0.4))

# =======================================================
# Visualización en Z
# =======================================================
plot_z_real_vs_syn <- function(Z_fit, y_fit, Z_syn, y_syn){
  df_r <- data.frame(Z1=Z_fit[,1], Z2=Z_fit[,2], Grupo=y_fit, Set="Real_train")
  df_s <- data.frame(Z1=Z_syn[,1], Z2=Z_syn[,2], Grupo=y_syn, Set="Sintetico")
  df   <- rbind(df_r, df_s)
  p1 <- ggplot(df, aes(Z1, Z2, color=Grupo, shape=Set)) +
    geom_point(alpha=0.7, size=2.8) +
    stat_ellipse(aes(group=interaction(Grupo,Set)), linetype=2, alpha=0.25) +
    theme_minimal(base_size=14) +
    labs(title="Z=2 MOFA: Reales vs Sintéticos", subtitle="Puntos y elipses por Grupo/Set")
  print(p1)
}
plot_z_real_vs_syn(Z_fit, y_fit, syn$Z, syn$y)

# =======================================================
# Proyección + Modelos
# =======================================================
W <- do.call(rbind, pesos_list)
feat_order <- rownames(W)
X_test_blocks <- list(
  transcriptomica = tst_Data$tx$test,
  proteomica     = tst_Data$pr$test,
  metabolomica   = tst_Data$me$test,
  clinica        = tst_Data$cl$test
)
X_test_all <- do.call(rbind, X_test_blocks)
X_test_all <- X_test_all[feat_order, , drop=FALSE]
X_test_lat <- t(X_test_all) %*% W %*% MASS::ginv(crossprod(W))

y_test <- tst_Data$grupo_test$grupo
train_y <- as.factor(as.character(train_y))
y_test  <- as.factor(as.character(y_test))

# Multinom
df_tr_lat <- data.frame(train_latent_X, y=train_y)
mn_fit <- nnet::multinom(y ~ ., data=df_tr_lat, trace=FALSE, maxit=1000)
colnames(X_test_lat) <- c("Factor1","Factor2")
probs_mn <- predict(mn_fit, newdata=as.data.frame(X_test_lat), type="probs")

# Lasso
cv_glm <- cv.glmnet(as.matrix(train_latent_X), train_y, family="multinomial",
                    alpha=1, type.measure="class", weights=w_train, nfolds=5)
probs_lasso <- predict(cv_glm, newx=as.matrix(X_test_lat), s="lambda.1se", type="response")[,,1]

# Calibración (opcional)
CALIBRAR <- FALSE
P_val <- predict(cv_glm, newx=as.matrix(Z_val), s="lambda.1se", type="response")[,,1]
cal    <- platt_fit(y_val, as.matrix(P_val))
P_te   <- probs_lasso
P_te_c <- if (CALIBRAR) platt_apply(cal, as.matrix(P_te)) else probs_lasso

# Métricas + Confusion Matrices
mm_mn     <- metrics_multiclass2(y_test, probs_mn)
mm_lasso  <- metrics_multiclass2(y_test, probs_lasso)
mm_lassoC <- metrics_multiclass2(y_test, P_te_c)
print(mm_mn); print(mm_lasso); print(mm_lassoC)

pred_mn     <- factor(colnames(probs_mn)[max.col(probs_mn)], levels=levels(y_test))
pred_lasso  <- factor(colnames(probs_lasso)[max.col(probs_lasso)], levels=levels(y_test))
pred_lassoC <- factor(colnames(P_te_c)[max.col(P_te_c)], levels=levels(y_test))
print(caret::confusionMatrix(pred_mn, y_test))
print(caret::confusionMatrix(pred_lasso, y_test))
print(caret::confusionMatrix(pred_lassoC, y_test))
