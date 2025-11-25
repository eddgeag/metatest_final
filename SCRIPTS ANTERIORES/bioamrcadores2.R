# =============================================================
# PIPELINE DOBLE: (A) PREDICCIÓN   (B) DESCUBRIMIENTO DE BIOMARCADORES
# =============================================================
# - Datos esperados: df_train_final.rds, df_test_final.rds (con columna 'y')
# - Clases: NP / FD / SO (referencia NP)
# - Este script separa claramente la parte predictiva de la parte de biomarcadores.
# - "Mandos" para relajar / endurecer selección al inicio del script.
# =============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet)
  library(brms)
  library(caret)
  library(pROC)
  library(randomForest)
})

set.seed(123)

# -----------------------------
# 0) PARÁMETROS (MANDOS)
# -----------------------------
REF_CLASS <- "NP"
# Para preselección por LASSO (solo B)
RELAX_MIN_K <- 5     # nº mínimo de features para parar (más bajo = más laxo)
RELAX_MAX_K <- 100   # nº máximo permitido
N_FOLDS_GLMENT <- 3  # folds para cv.glmnet con N pequeño
# Umbral de estabilidad (bootstrap glmnet)
STAB_B <- 50
STAB_THRESHOLD <- 0.2  # proporción (0.2 = 20%)
# Priors BRMS (más grandes = más laxos)
BRMS_SD_B <- 1
BRMS_SD_INT <- 5
# ¿BRMS usa todas las variables (TRUE) o solo las seleccionadas por LASSO (FALSE)?
BRMS_USE_ALL <- TRUE   # A: para predicción mantenemos todo (recupera performance)

# -----------------------------
# 1) UTILIDADES COMUNES
# -----------------------------
refactor_levels <- function(x) factor(as.character(x))

sanitize_names <- function(v){
  v2 <- gsub("[^A-Za-z0-9]+", "_", v)
  v2 <- gsub("__+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- ifelse(grepl("^[0-9]", v2), paste0("X", v2), v2)
  v2 <- make.names(v2, unique = FALSE)
  v2 <- gsub("\\.+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  make.unique(v2, sep = "_")
}

metrics_from_probs <- function(prob_mat, y_true){
  lev <- levels(y_true); prob_mat <- prob_mat[, lev, drop=FALSE]
  pred_labels <- apply(prob_mat, 1, function(p) lev[which.max(p)])
  pred_class  <- factor(pred_labels, levels=lev)
  cm <- caret::confusionMatrix(pred_class, y_true)
  auc_obj <- tryCatch(pROC::multiclass.roc(y_true, prob_mat), error=function(e) NULL)
  auc_num <- if(is.null(auc_obj)) NA_real_ else as.numeric(auc_obj$auc)
  Y <- model.matrix(~ y_true - 1); colnames(Y) <- lev
  eps <- 1e-15
  logloss <- -mean(rowSums(Y * log(pmax(pmin(prob_mat, 1-eps), eps))))
  tibble(Accuracy=as.numeric(cm$overall["Accuracy"]),
         BalAcc  = mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE),
         Kappa   = as.numeric(cm$overall["Kappa"]),
         AUC=auc_num, LogLoss=logloss)
}

brms_probs <- function(fit, newdata){
  pp <- posterior_epred(fit, newdata=newdata)
  M  <- apply(pp, c(2,3), mean)
  cls <- dimnames(pp)[[3]]; if(is.null(cls)) cls <- fit$family$names
  colnames(M) <- cls; M
}

make_multinom_priors <- function(y_levels, ref=REF_CLASS, sd_b=BRMS_SD_B, sd_int=BRMS_SD_INT){
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(set_prior(paste0("student_t(3,0,", sd_b, ")"), class="b", dpar=dpar),
      set_prior(paste0("student_t(3,0,", sd_int,")"), class="Intercept", dpar=dpar))
  }))
}

# -----------------------------
# 2) CARGA Y ARMONIZACIÓN DE DATOS
# -----------------------------

df_train <- readRDS("./df_train_final.rds")
df_test  <- readRDS("./df_test_final.rds")

df_train$y <- stats::relevel(refactor_levels(df_train$y), ref = REF_CLASS)
df_test$y  <- factor(df_test$y, levels = levels(df_train$y))

common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
new_names <- make.unique(sanitize_names(common_cols), sep="_")
colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
colnames(df_test)[match(common_cols,  colnames(df_test))]  <- new_names

df_train <- df_train[, c("y", new_names), drop=FALSE]
df_test  <- df_test[,  c("y", new_names), drop=FALSE]

# =============================================================
# A) MODELO PREDICTIVO (SIN PRESELECCIÓN AGRESIVA)
#    - Mantiene TODAS las variables para maximizar performance con N pequeño
#    - Reporta: BRMS + varios clasificadores clásicos
# =============================================================

cat("\n================ (A) PREDICCIÓN ================\n")

# --- BRMS usando todas las variables ---
priors_all <- make_multinom_priors(levels(df_train$y))
ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

fit_brm_all <- brms::brm(
  y ~ ., data = df_train, family = brms::categorical(), prior = priors_all,
  chains = 4, iter = 8000, warmup = 4000, refresh = 0, control = ctrl
)

cat("\nResumen BRMS (todas las features)\n"); print(summary(fit_brm_all))
cat("\nLOO/WAIC\n"); print(loo::loo(fit_brm_all, moment_match=TRUE)); print(loo::waic(fit_brm_all))

P_all <- brms_probs(fit_brm_all, df_test)[, levels(df_test$y), drop=FALSE]
cat("\nMétricas TEST (BRMS-todo)\n"); print(metrics_from_probs(P_all, df_test$y))
cat("\nMatriz de confusión (BRMS-todo)\n"); print(caret::confusionMatrix(
  factor(colnames(P_all)[max.col(P_all, ties.method="first")], levels=levels(df_test$y)),
  df_test$y
))

# --- Modelos clásicos con CV (todas las variables) ---
# Control CV
multiSummary <- function (data, lev = NULL, model = NULL) {
  acc <- caret::postResample(data$pred, data$obs)[["Accuracy"]]
  kap <- caret::postResample(data$pred, data$obs)[["Kappa"]]
  cm  <- caret::confusionMatrix(data$pred, data$obs)
  bal <- mean(cm$byClass[, "Balanced Accuracy"], na.rm=TRUE)
  prob <- as.matrix(data[, lev])
  auc <- tryCatch(pROC::multiclass.roc(data$obs, prob)$auc, error=function(e) NA)
  Y <- model.matrix(~ data$obs - 1); colnames(Y) <- lev
  eps <- 1e-15
  logloss <- -mean(rowSums(Y * log(pmax(pmin(prob, 1-eps), eps))))
  c(Accuracy=acc, BalAcc=bal, Kappa=kap, AUC=auc, LogLoss=logloss)
}

ctrl_cv <- trainControl(method="cv", number=3, classProbs=TRUE,
                        summaryFunction=multiSummary, savePredictions="final")

set.seed(123)
fit_ridge_all <- train(y ~ ., data=df_train, method="glmnet", trControl=ctrl_cv,
                       tuneGrid=expand.grid(alpha=0, lambda=10^seq(-3,1,length=20)))
fit_enet_all  <- train(y ~ ., data=df_train, method="glmnet", trControl=ctrl_cv,
                       tuneLength=10)
fit_svm_all   <- train(y ~ ., data=df_train, method="svmRadial", trControl=ctrl_cv, tuneLength=5)
fit_nb_all    <- train(y ~ ., data=df_train, method="nb",         trControl=ctrl_cv, tuneLength=5)
fit_rf_all    <- train(y ~ ., data=df_train, method="rf",         trControl=ctrl_cv, tuneLength=5)

models_pred <- list(Ridge=fit_ridge_all, ENet=fit_enet_all, SVM=fit_svm_all,
                    NB=fit_nb_all, RF=fit_rf_all)

pred_results <- bind_rows(lapply(models_pred, function(fit){
  pr <- predict(fit, newdata=df_test, type="prob")
  metrics_from_probs(as.matrix(pr), df_test$y)
}), .id="Model")

cat("\nResultados TEST (clásicos, todas las vars)\n"); print(pred_results)

# =============================================================
# B) DESCUBRIMIENTO DE BIOMARCADORES (SELECCIÓN + ESTABILIDAD + IMPORTANCIA)
#    - Mantiene la parte interpretativa separada de la predictiva
# =============================================================

cat("\n================ (B) BIOMARCADORES ================\n")

# --- 1) LASSO multinomial RELAJADO para preselección ---
class_counts <- table(df_train$y)
df_train_bal <- if(any(class_counts < max(class_counts))) {
  caret::upSample(df_train[,-1], df_train$y, yname="y")
} else df_train

x_tr <- as.matrix(df_train_bal[, setdiff(names(df_train_bal), "y")])
y_tr <- factor(df_train_bal$y, levels=levels(df_train$y))

cvfit <- cv.glmnet(x_tr, y_tr, family="multinomial", type.multinomial="grouped",
                   alpha=1, standardize=TRUE, nfolds=N_FOLDS_GLMENT)

get_nonzero <- function(fit, lambda){
  coefs <- coef(fit, s=lambda)
  feats <- unique(unlist(lapply(coefs, function(m){
    nm <- rownames(m); if(is.null(nm) || nrow(m)<=1) return(character(0))
    nz <- as.vector(m[-1, , drop=FALSE]) != 0
    nm <- nm[-1]; nm[nz]
  })))
  feats[feats %in% colnames(x_tr)]
}

relaxed_select <- function(fit, min_k=RELAX_MIN_K, max_k=RELAX_MAX_K){
  lambdas <- sort(fit$lambda, decreasing=FALSE)
  for(l in lambdas){
    nz <- get_nonzero(fit, l)
    if(length(nz) >= min_k) return(head(nz, max_k))
  }
  nz_min <- get_nonzero(fit, min(lambdas))
  if(length(nz_min)) return(head(nz_min, max_k))
  nz_1se <- get_nonzero(fit, cvfit$lambda.1se)
  if(length(nz_1se)) return(head(nz_1se, max_k))
  nz_min2 <- get_nonzero(fit, cvfit$lambda.min)
  if(length(nz_min2)) return(head(nz_min2, max_k))
  character(0)
}

sel <- relaxed_select(cvfit$glmnet.fit)
cat("\n# Features seleccionados (relajado): ", length(sel), "\n"); print(sel)

# Datasets reducidos SOLO para bloque de biomarcadores
keep_bio <- unique(c("y", sel))
df_train_sel <- df_train[, keep_bio, drop=FALSE]
df_test_sel  <- df_test[,  keep_bio, drop=FALSE]

# --- 2) Estabilidad de selección (bootstrap + rebalanceo) ---
bootstrap_lasso <- function(df, B=STAB_B){
  out <- vector("list", B)
  levs <- levels(df$y)
  for(b in seq_len(B)){
    idx <- sample(seq_len(nrow(df)), size=nrow(df), replace=TRUE)
    dfb <- df[idx, , drop=FALSE]
    cc <- table(dfb$y)
    if(any(cc < max(cc))) dfb <- caret::upSample(dfb[,-1], dfb$y, yname="y")
    x <- as.matrix(dfb[, setdiff(names(dfb), "y")])
    y <- factor(dfb$y, levels=levs)
    nfolds <- max(3, min(5, length(y)))
    cvb <- try(cv.glmnet(x, y, family="multinomial", type.multinomial="grouped",
                         alpha=1, standardize=TRUE, nfolds=nfolds), silent=TRUE)
    if(inherits(cvb, "try-error")) { out[[b]] <- tibble(feature=character(0)); next }
    selb <- get_nonzero(cvb$glmnet.fit, cvb$lambda.1se)
    if(length(selb) == 0) selb <- get_nonzero(cvb$glmnet.fit, cvb$lambda.min)
    out[[b]] <- tibble(feature=selb)
  }
  bind_rows(out, .id="boot") %>% count(feature, name="freq") %>% arrange(desc(freq))
}

stab_tbl <- if(length(sel) > 0) bootstrap_lasso(df_train_sel, B=STAB_B) else tibble(feature=character(0), freq=integer(0))
cat("\nTop estabilidad (bootstrap glmnet)\n"); print(head(stab_tbl, 20))

stable_candidates <- stab_tbl %>% filter(freq >= STAB_THRESHOLD*STAB_B) %>% pull(feature)
cat("\n# Candidatos estables (>", STAB_THRESHOLD*100, "%) : ", length(stable_candidates), "\n", sep=""); print(stable_candidates)

# --- 3) Importancia de variables (Ridge/ENet/RF) sobre df_train_sel ---
coef_importance <- function(fit){
  if(!inherits(fit$finalModel, "glmnet") && !inherits(fit$finalModel, "cv.glmnet")) return(NULL)
  lambda <- fit$bestTune$lambda
  coefs <- try(coef(fit$finalModel, s=lambda), silent=TRUE)
  if(inherits(coefs, "try-error")) coefs <- coef(fit$finalModel)
  if(is.list(coefs)){
    v <- Reduce(`+`, lapply(coefs, function(m){
      a <- abs(as.numeric(m[-1, , drop=FALSE])); names(a) <- rownames(m)[-1]; a
    }))
  } else {
    v <- abs(as.numeric(coefs[-1])); names(v) <- rownames(coefs)[-1]
  }
  sort(v, decreasing=TRUE)
}

# Entrenar modelos de importancia SOLO con features seleccionadas
if(length(sel) >= 2){
  set.seed(123)
  fit_ridge_bio <- train(y ~ ., data=df_train_sel, method="glmnet", trControl=ctrl_cv,
                         tuneGrid=expand.grid(alpha=0, lambda=10^seq(-3,1,length=20)))
  fit_enet_bio  <- train(y ~ ., data=df_train_sel, method="glmnet", trControl=ctrl_cv, tuneLength=10)
  fit_rf_bio    <- train(y ~ ., data=df_train_sel, method="rf",     trControl=ctrl_cv, tuneLength=5)
  
  imp_ridge <- coef_importance(fit_ridge_bio)
  imp_enet  <- coef_importance(fit_enet_bio)
  imp_rf    <- varImp(fit_rf_bio, scale=TRUE)
  
  cat("\nImportancia Ridge (abs coef)\n"); print(head(imp_ridge, 20))
  cat("\nImportancia ENet (abs coef)\n"); print(head(imp_enet, 20))
  cat("\nImportancia RF\n"); print(imp_rf)
  
  # Ranking consolidado (maneja NULLs)
  names_ridge <- if(length(imp_ridge)) names(imp_ridge) else character(0)
  names_enet  <- if(length(imp_enet))  names(imp_enet)  else character(0)
  names_rf    <- rownames(imp_rf$importance)
  all_features <- unique(c(names_ridge, names_enet, names_rf))
  
  importance_tbl <- data.frame(
    Feature = all_features,
    Ridge   = if(length(imp_ridge)) imp_ridge[all_features] else rep(NA_real_, length(all_features)),
    ENet    = if(length(imp_enet))  imp_enet[all_features]  else rep(NA_real_, length(all_features)),
    RF      = imp_rf$importance[all_features, 1]
  ) %>%
    mutate(across(-Feature, ~replace_na(., 0))) %>%
    mutate(AverageRank = (rank(-Ridge, ties.method="average") +
                            rank(-ENet,  ties.method="average") +
                            rank(-RF,    ties.method="average"))/3) %>%
    arrange(AverageRank)
  
  cat("\nRanking consolidado de importancia\n"); print(importance_tbl)
  
} else {
  cat("\n[AVISO] Muy pocas features seleccionadas para importancia consolidada.\n")
}

# --- 4) BRMS orientado a biomarcadores (opcional, sobre subset estable) ---
if(length(stable_candidates) >= 1){
  keep2 <- c("y", intersect(stable_candidates, colnames(df_train_sel)))
  df_train_stab <- df_train_sel[, keep2, drop=FALSE]
  df_test_stab  <- df_test_sel[,  keep2, drop=FALSE]
  pri2 <- make_multinom_priors(levels(df_train_stab$y))
  fit_brm_stab <- brms::brm(
    y ~ ., data = df_train_stab, family = brms::categorical(), prior = pri2,
    chains = 4, iter = 6000, warmup = 3000, refresh = 0, control = ctrl
  )
  P2 <- brms_probs(fit_brm_stab, df_test_stab)[, levels(df_test_stab$y), drop=FALSE]
  cat("\nMétricas TEST (BRMS con features estables)\n"); print(metrics_from_probs(P2, df_test_stab$y))
}

cat("\n================ FIN ====================\n")
