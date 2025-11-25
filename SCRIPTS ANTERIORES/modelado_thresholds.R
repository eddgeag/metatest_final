# diagnóstico, visualización y modelos
# =======================================================

rm(list = ls())

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
  library(MOFA2)
  library(AnnotationDbi)
})

set.seed(123)
select_features_view <- function(mod, view, q_thr, data_list, metadata) {
  aux <- plot_top_weights(mod, view = view, factors = 1:2,
                          abs = TRUE, scale = TRUE, nfeatures = Inf)
  
  # threshold
  thr <- quantile(aux$data$value, probs = q_thr, na.rm = TRUE)
  
  # histograma con threshold
  hist(aux$data$value,
       main = paste("View", view, "- threshold =", q_thr),
       xlab = "Weight values", col = "grey", border = "white")
  abline(v = thr, col = "red", lwd = 2)
  
  # features seleccionados
  feats <- aux$data$feature[aux$data$value >= thr]
  
  # subset train y test
  sub_train <- data_list$train[as.character(feats), , drop = FALSE]
  sub_test  <- data_list$test[as.character(feats), , drop = FALSE]
  
  # mapeo Entrez->Symbolmes(metadata) & "GeneSymbol" %in% colnames(metadata)) {
    map_dic <- metadata %>% dplyr::select(EntrezGeneID, GeneSymbol) %>% distinct()
    idx <- rownames(sub_train)
    ids_map <- data.frame(EntrezGeneID = idx) %>%
      left_join(map_dic, by = "EntrezGeneID")
    new_ids <- ifelse(!is.na(ids_map$GeneSymbol), ids_map$GeneSymbol, ids_map$EntrezGeneID)
    rownames(sub_train) <- new_ids
    rownames(sub_test)  <- new_ids

  
  return(list(train = sub_train, test = sub_test))
}
project_all_views <- function(X_list, W_list) {
  W_all <- NULL
  X_all <- NULL
  
  for (v in names(W_list)) {
    Xv <- X_list[[v]]
    if (is.list(Xv) && "test" %in% names(Xv)) {
      Xv <- Xv$test
    }
    common_feats <- intersect(rownames(W_list[[v]]), rownames(Xv))
    W_sub <- as.matrix(W_list[[v]][common_feats, , drop=FALSE])
    X_sub <- as.matrix(Xv[common_feats, , drop=FALSE])
    W_all <- rbind(W_all, W_sub)
    X_all <- rbind(X_all, X_sub)
  }
  
  Z_est <- solve(t(W_all) %*% W_all) %*% t(W_all) %*% X_all
  Z_est <- t(Z_est)
  return(Z_est)
}

refactor_levels <- function(x) as.factor(as.character(x))



build_sets_from_grid <- function(thr_grid){
  out_train <- vector("list", nrow(thr_grid))
  out_test  <- vector("list", nrow(thr_grid))
  X_list <- list(
    transcriptomica = tst_Data$tx$test,
    proteomica      = tst_Data$pr$test,
    metabolomica    = tst_Data$me$test,
    clinical        = tst_Data$cl$test
  )
  factores_Test <- project_all_views(X_list, pesos_list)
  
  for(i in seq_len(nrow(thr_grid))){
    q_tx <- thr_grid$q_tx[i]; q_pr <- thr_grid$q_pr[i]
    q_me <- thr_grid$q_me[i]; q_cl <- thr_grid$q_cl[i]
    
    # aplicar selección de features
    tx <- select_features_view(mod, view = 1, q_thr = q_tx,
                               data_list = tst_Data$tx,
                               metadata  = tst_Data$features_metadata)
    pr <- select_features_view(mod, view = 2, q_thr = q_pr,
                               data_list = tst_Data$pr,
                               metadata  = tst_Data$features_metadata)
    me <- select_features_view(mod, view = 3, q_thr = q_me,
                               data_list = tst_Data$me,
                               metadata  = tst_Data$features_metadata)
    cl <- select_features_view(mod, view = 4, q_thr = q_cl,
                               data_list = tst_Data$cl,
                               metadata  = tst_Data$features_metadata)
    
    # concatenar
    train_feats   <- rbind(tx$train, pr$train, me$train, cl$train)
    train_feats_t <- as.data.frame(t(train_feats))
    final_train   <- cbind(Z_seed, train_feats_t)
    
    test_feats   <- rbind(tx$test, pr$test, me$test, cl$test)
    test_feats_t <- as.data.frame(t(test_feats))
    final_test   <- cbind(factores_Test, test_feats_t)
    
    out_train[[i]] <- final_train
    out_test[[i]]  <- final_test
  }
  
  list(
    train   = out_train,
    test    = out_test,
    y_train = y_seed,
    y_test  = tst_Data$grupo_test$grupo
  )
}

# =======================================================
# Uso
# =======================================================
# =======================================================
# Inputs
# =======================================================
dats       <- readRDS("./datos_para_modelar.rds")
tst_Data   <- readRDS("./ready_for_modeling.rds")
factors    <- dats$factores
pesos_list <- dats$pesos

Z_seed <- as.matrix(factors %>% filter(Tipo=="Real") %>% dplyr::select(where(is.numeric)))[,1:2,drop=FALSE]
y_seed <- factors %>% filter(Tipo=="Real") %>% pull(Grupo) %>% factor()
names(y_seed) <- rownames(Z_seed)
lev <- levels(y_seed)

mod <- readRDS("./modelo.rds")

v <- c(0.8, 0.85, 0.9, 0.95)
thr_grid <- expand.grid(
  q_tx = v,
  q_pr = v,
  q_me = v,
  q_cl = v
)


return.list <- build_sets_from_grid(thr_grid)

detach("package:MOFA2", unload=TRUE)


X_train <- return.list$train[[1]]
X_test <- return.list$test[[1]]
y_train <- refactor_levels(return.list$y_train)
y_test <- refactor_levels(return.list$y_test)
# 1) Referencia explícita
# referencia
y_train <- relevel(refactor_levels(return.list$y_train), ref = "NP")
y_test  <- relevel(refactor_levels(return.list$y_test),  ref = "NP")

df_train <- data.frame(y = y_train, X_train)
df_test  <- data.frame(y = y_test,  X_test)

# priors solo para no-referencia
make_multinom_priors <- function(y_levels, ref = "NP", sd_b = 0.5, sd_int = 5){
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(
      set_prior(paste0("normal(0,", sd_b, ")"), class = "b",         dpar = dpar),
      set_prior(paste0("normal(0,", sd_int,")"), class = "Intercept", dpar = dpar)
    )
  }))
}
priors <- make_multinom_priors(levels(df_train$y), ref = "NP", sd_b = 0.5, sd_int = 5)



suppressPackageStartupMessages({
  library(brms); library(nnet); library(caret); library(pROC); library(dplyr)
})

# ---------- utilidades ----------
metrics_from_probs <- function(prob_mat, y_true){
  lev <- levels(y_true)
  stopifnot(setequal(colnames(prob_mat), lev))
  prob_mat <- prob_mat[, lev, drop = FALSE]
  
  pred_labels <- apply(prob_mat, 1, function(p) lev[which.max(p)])
  pred_class  <- factor(pred_labels, levels = lev)
  
  cm <- caret::confusionMatrix(pred_class, y_true)
  
  auc_obj <- tryCatch(pROC::multiclass.roc(y_true, prob_mat), error = function(e) NULL)
  auc_num <- if (is.null(auc_obj)) NA_real_ else as.numeric(auc_obj$auc)
  
  Y   <- model.matrix(~ y_true - 1); colnames(Y) <- lev
  eps <- 1e-15
  logloss <- -mean(rowSums(Y * log(pmax(pmin(prob_mat, 1 - eps), eps))))
  
  tibble(
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    BalAcc   = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
    Kappa    = as.numeric(cm$overall["Kappa"]),
    AUC      = auc_num,
    LogLoss  = logloss
  )
}

brms_probs <- function(fit, newdata){
  pp <- posterior_epred(fit, newdata = newdata)         # draws x N x K
  M  <- apply(pp, c(2,3), mean)
  cls <- dimnames(pp)[[3]]; if (is.null(cls)) cls <- fit$family$names
  stopifnot(ncol(M) == length(cls))
  colnames(M) <- cls
  M
}

# ---------- BRMS MCMC ----------
probs_mcmc <- brms_probs(fit, df_test)
m_mcmc <- metrics_from_probs(probs_mcmc, df_test$y)

# ---------- BRMS VB ----------
probs_vb <- brms_probs(fit_vb, df_test)
m_vb <- metrics_from_probs(probs_vb, df_test$y)

# ---------- NNET multinomial ----------
fit_nnet <- nnet::multinom(y ~ ., data = df_train, trace = FALSE, MaxNWts = 20000)
probs_nnet <- predict(fit_nnet, newdata = df_test, type = "probs")
if (is.null(colnames(probs_nnet))) colnames(probs_nnet) <- levels(df_test$y)
m_nnet <- metrics_from_probs(probs_nnet, df_test$y)

# ---------- Comparación y “mejor” ----------
res_comp <- bind_rows(
  MCMC = m_mcmc,
  VB   = m_vb,
  NNET = m_nnet,
  .id = "Model"
)

# ranking: mayor BalAcc, menor LogLoss, mayor AUC
ranked <- res_comp %>%
  mutate(
    r1 = rank(-BalAcc,  ties.method = "min"),
    r2 = rank( LogLoss, ties.method = "min"),
    r3 = rank(-AUC,     ties.method = "min"),
    RankSum = r1 + r2 + r3
  ) %>%
  arrange(RankSum, desc(BalAcc), AUC, LogLoss)

print(res_comp)
cat("\nMejor según RankSum:\n")
print(ranked[1, ])



