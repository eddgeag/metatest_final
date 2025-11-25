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
# util: rbind seguro con NULL
rb <- function(a,b){
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  rbind(a,b)
}

# =======================================================
# Wrapper de síntesis: controla dispersión y varianza
# =======================================================
# Distancia de Bhattacharyya entre Gaussianas
bhatt_sigma <- function(mu1, S1, mu2, S2, eps=1e-4){
  p <- length(mu1)
  S1r <- S1 + diag(eps, p); S2r <- S2 + diag(eps, p)
  S   <- 0.5*(S1r+S2r)
  invS <- tryCatch(solve(S), error=function(e) MASS::ginv(S))
  term1 <- 0.125 * t(mu1-mu2) %*% invS %*% (mu1-mu2)
  term2 <- 0.5 * log( det(S) / sqrt(pmax(det(S1r),1e-12)*pmax(det(S2r),1e-12)) )
  as.numeric(term1 + term2)
}

diag_sinteticos <- function(Z_fit, y_fit, Z_syn, y_syn, S_pool){
  lev <- levels(y_fit)
  do.call(rbind, lapply(lev, function(lv){
    Zr <- Z_fit[y_fit==lv,,drop=FALSE]; Zs <- Z_syn[y_syn==lv,,drop=FALSE]
    if (nrow(Zs)<2) return(NULL)
    mu_r <- colMeans(Zr); mu_s <- colMeans(Zs)
    S_ref <- if (nrow(Zr)>=3) cov(Zr) else S_pool
    S_s   <- cov(Zs)
    data.frame(
      clase=lv, n_real=nrow(Zr), n_syn=nrow(Zs),
      mu_L2 = round(sqrt(sum((mu_r-mu_s)^2)), 4),
      cov_F = round(norm(S_ref - S_s, "F"), 4),
      Bhattacharyya = round(bhatt_sigma(mu_r, S_ref, mu_s, S_s, eps=1e-4), 4)
    )
  }))
}

# úsalo así:
S_pool <- cov(Z_fit)

# whitening-recoloring exacto a mu/Sigma objetivo
match_moments <- function(X, mu_tgt, S_tgt){
  Sx <- cov(X)
  Ex <- eigen(Sx, symmetric=TRUE); Et <- eigen(S_tgt, symmetric=TRUE)
  Wx <- Ex$vectors %*% diag(1/sqrt(pmax(Ex$values,1e-8))) %*% t(Ex$vectors)
  Ct <- Et$vectors %*% diag(sqrt(pmax(Et$values,1e-8))) %*% t(Et$vectors)
  Zc <- sweep(X, 2, colMeans(X), "-") %*% Wx %*% Ct
  sweep(Zc, 2, mu_tgt, "+")
}

# mu+Sigma de referencia por clase, robusto en n<3
ref_mu_S <- function(Zc, Zpool=S_pool, shrink=0.2){
  mu <- colMeans(Zc)
  Sref <- if (nrow(Zc) >= 3) cov(Zc) else Zpool
  # shrink a la diagonal para estabilidad
  Sref <- (1-shrink)*Sref + shrink*diag(diag(Sref), ncol(Sref))
  list(mu=mu, S=Sref)
}

# muestreo Gaussiano truncado por radio de Mahalanobis
rgauss_trunc <- function(n, mu, S, rad_q=0.8){
  p <- length(mu)
  # umbral por cuantiles empíricos o Chi2
  thr <- qchisq(rad_q, df=p)
  X <- MASS::mvrnorm(n*2, mu, S)  # sobre-muestrea
  D2 <- mahalanobis(X, center=mu, cov=S)
  keep <- which(D2 <= thr)
  if (length(keep) < n) {
    # reescala radialmente los que exceden
    over <- setdiff(seq_len(nrow(X)), keep)
    X[over,] <- sweep(X[over,]-matrix(mu, nrow=length(over), ncol=p, byrow=TRUE),
                      1, sqrt(D2[over]/thr), `/`) + matrix(mu, nrow=length(over), ncol=p, byrow=TRUE)
    keep <- seq_len(nrow(X))
  }
  X[sample(keep, n), , drop=FALSE]
}

# --- Generador por clase con control estricto de dispersión ---
gen_class_matched <- function(Zc, n_new, shrink=0.2, rad_q=0.8, Zpool=S_pool){
  RS <- ref_mu_S(Zc, Zpool=Zpool, shrink=shrink)
  # evitar SD diminutas o excesivas: acota autovalores
  Et <- eigen(RS$S, symmetric=TRUE)
  vals <- pmin(pmax(Et$values, 0.05*median(Et$values)), 5*median(Et$values))
  S_tgt <- Et$vectors %*% diag(vals) %*% t(Et$vectors)
  X <- rgauss_trunc(n_new, RS$mu, S_tgt, rad_q=rad_q)
  match_moments(X, RS$mu, S_tgt)
}

# --- Wrapper de síntesis: usa gen_class_matched; SMOTE sólo si n>=3 ---
synth_with_controls <- function(Z_fit, y_fit, targets,
                                shrink=0.2, rad_q=0.8){
  lev <- levels(y_fit); Zs <- list(); ys <- character(0)
  for (lv in lev){
    Zc <- Z_fit[y_fit==lv,,drop=FALSE]; nr <- nrow(Zc)
    need <- targets[[lv]] - nr
    if (need <= 0) next
    if (nr < 3){
      # clase rara: sólo Gauss truncada match-moments con Sigma pool
      Xnew <- gen_class_matched(Zc, need, shrink=shrink, rad_q=rad_q, Zpool=S_pool)
    } else {
      # clase normal: Gauss trunc + luego matching exacto
      Xnew <- gen_class_matched(Zc, need, shrink=shrink, rad_q=rad_q, Zpool=S_pool)
    }
    Zs[[length(Zs)+1]] <- Xnew
    ys <- c(ys, rep(lv, nrow(Xnew)))
  }
  Zs <- if (length(Zs)) do.call(rbind,Zs) else NULL
  list(Z=Zs, y=factor(ys, levels=lev))
}

tbl <- table(y_seed)
k <- 5L; min_fold <- 8L
min_total  <- ceiling(min_fold * k / (k - 1))  # si lo usas luego
max_ratio  <- 3L; base_target <- 15L
raw_target <- pmin(base_target, as.integer(tbl) * max_ratio)
targets    <- setNames(pmax(raw_target, 40L), names(tbl))

grid <- expand.grid(
  shrink   = c(0.01, 0.02,0.03, 0.04,0.05),
  rad_q    = c(0.3,0.4,0.5,0.6, 0.7, 0.8)
)
# targets equilibrados moderados
targets <- setNames(rep(40L, length(levels(y_fit))), levels(y_fit))
score_synthesis <- function(syn, Z_fit, y_fit){
  dt <- diag_sinteticos(Z_fit, y_fit, syn$Z, syn$y, S_pool=S_pool)
  if (is.null(dt) || !nrow(dt)) return(Inf)
  mean(dt$Bhattacharyya, na.rm=TRUE) +
    0.1*mean(dt$mu_L2, na.rm=TRUE) +
    0.05*mean(dt$cov_F, na.rm=TRUE)
}
syn <- synth_with_controls(Z_fit, y_fit, targets, shrink=0.3, rad_q=0.8)

diag_tbl <- diag_sinteticos(Z_fit, y_fit, syn$Z, syn$y, S_pool=S_pool)
print(diag_tbl)

results <- purrr::pmap_dfr(grid, function(shrink, rad_q){
  syn <- try(synth_with_controls(Z_fit, y_fit, targets = targets,shrink = shrink,rad_q = rad_q), silent=TRUE)
  if (inherits(syn,"try-error") || is.null(syn$Z)) return(NULL)
  sc <- score_synthesis(syn, Z_fit, y_fit)
  data.frame(shrink, rad_q, score=sc)
})

best <- results[which.min(results$score), , drop=FALSE]
print(best)

# usa proporciones reales (o las esperadas en producción)
pri <- prop.table(table(y_seed))  # o fija manualmente
Ntot <- 120L
targets <- (Ntot * pri) |> round() |> as.integer() |> setNames(names(pri))





syn <- synth_with_controls(Z_fit, y_fit, targets, shrink=best$shrink, rad_q=best$rad_q)

train_latent_X <- rbind(Z_fit, syn$Z)
train_y        <- factor(c(as.character(y_fit), as.character(syn$y)), levels=levels(y_fit))
w_train <- rep(1, nrow(train_latent_X))  # quita sobrepesos a FD


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

diag_tbl <- diag_sinteticos(Z_fit, y_fit, syn$Z, syn$y,S_pool = S_pool)

diag_tbl
# priors reales en test vs. priors de entrenamiento (con sintéticos)
pi_train <- prop.table(table(train_y))
pi_test  <- prop.table(table(y_test))  # o estima de producción si la conoces

# corrige posterior por cociente de priors y renormaliza
probs_lasso_adj <- rebalance_post(probs_lasso, prop.table(table(train_y)), prop.table(table(y_test)))

mm_lasso_adj <- metrics_multiclass2(y_test, probs_lasso_adj)
print(mm_lasso_adj)
pred_lasso_adj <- factor(colnames(probs_lasso_adj)[max.col(probs_lasso_adj)], levels=levels(y_test))
print(caret::confusionMatrix(pred_lasso_adj, y_test))


cv_glm <- cv.glmnet(as.matrix(train_latent_X), train_y,
                    family="multinomial", alpha=1,
                    type.measure="deviance", nfolds=5)

P_val <- predict(cv_glm, newx=as.matrix(Z_val), s="lambda.1se", type="response")[,,1]
cal   <- platt_fit(y_val, as.matrix(P_val))
P_te  <- predict(cv_glm, newx=as.matrix(X_test_lat), s="lambda.1se", type="response")[,,1]
P_te  <- platt_apply(cal, as.matrix(P_te))                 # activa calibración
P_te  <- rebalance_post(P_te, prop.table(table(train_y)),  # corrige priors
                        prop.table(table(y_test)))
mm <- metrics_multiclass2(y_test, P_te)

mm


# --- decisión óptima por clase (umbral/ponderación) ---
levs <- colnames(P_val)

bal_acc_from_weights <- function(w, P, y_true){
  W <- matrix(w, nrow=nrow(P), ncol=length(w), byrow=TRUE)
  S <- P * W
  preds <- factor(levs[max.col(S, ties.method="first")], levels=levels(y_true))
  cm <- caret::confusionMatrix(preds, y_true)
  mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE)
}

# optimiza en log-espacio para w_k > 0
init <- rep(0, length(levs))
opt <- optim(
  par = init,
  fn = function(par) -bal_acc_from_weights(exp(par), P_val, y_val),
  method = "Nelder-Mead",
  control = list(maxit=2000, reltol=1e-8)
)
w_star <- exp(opt$par); names(w_star) <- levs
print(w_star)

# aplica a test
Wte <- matrix(w_star, nrow=nrow(P_te), ncol=length(w_star), byrow=TRUE)
P_opt <- P_te * Wte
P_opt <- sweep(P_opt, 1, rowSums(P_opt), "/")

mm_opt <- metrics_multiclass2(y_test, P_opt)
print(mm_opt)

pred_opt <- factor(colnames(P_opt)[max.col(P_opt)], levels=levels(y_test))
print(caret::confusionMatrix(pred_opt, y_test))


X_tr <- rbind(Z_fit, Z_val); y_tr <- factor(c(as.character(y_fit), as.character(y_val)), levels=lev)
w_tr <- 1 / as.numeric(table(y_tr))[y_tr]; w_tr <- w_tr / mean(w_tr)

cv_glm <- cv.glmnet(as.matrix(X_tr), y_tr, family="multinomial",
                    alpha=1, type.measure="deviance", weights=w_tr, nfolds=5)
P_te  <- predict(cv_glm, newx=as.matrix(X_test_lat), s="lambda.1se", type="response")[,,1]
mm    <- metrics_multiclass2(y_test, P_te); print(mm)


library(klaR)
rda_fit <- klaR::rda(x=as.data.frame(X_tr), grouping=y_tr, gamma=0, lambda=0.4)  # tunea lambda ∈ [0,1]
P_rda   <- as.matrix(predict(rda_fit, as.data.frame(X_test_lat))$posterior)
mm_rda  <- metrics_multiclass2(y_test, P_rda); print(mm_rda)



tau <- 0.55  # sube/baja según validación
is_fd <- P_te[,"FD"] >= tau

# binario NP vs SO con LDA sobre (no-FD)
idx <- which(!is_fd)
lda_fit <- MASS::lda(x=as.data.frame(X_tr[y_tr!="FD",]), grouping=droplevels(y_tr[y_tr!="FD"]))
P_bin <- predict(lda_fit, as.data.frame(X_test_lat[idx,]))$posterior

pred <- rep("FD", nrow(P_te))
pred[idx] <- ifelse(P_bin[,"NP"] >= 0.5, "NP", "SO")
pred <- factor(pred, levels=levels(y_test))
print(caret::confusionMatrix(pred, y_test))


P_val <- predict(cv_glm, newx=as.matrix(Z_val), s="lambda.1se", type="response")[,,1]
cal   <- platt_fit(y_val, as.matrix(P_val))
P_teC <- platt_apply(cal, as.matrix(P_te))
P_teC <- rebalance_post(P_teC, prop.table(table(y_tr)), prop.table(table(y_test)))
print(metrics_multiclass2(y_test, P_teC))

