# ==========================
# MODELO FINAL: PRESELECCIÓN LASSO + BRMS + ESTABILIDAD (ROBUSTO EN N PEQUEÑO)
# ==========================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet)
  library(brms)
  library(caret)
  library(pROC)
})

set.seed(123)

# --- Utilidades ---
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

make_multinom_priors <- function(y_levels, ref="NP", sd_b=0.3, sd_int=3){
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(set_prior(paste0("student_t(3,0,", sd_b, ")"), class="b", dpar=dpar),
      set_prior(paste0("student_t(3,0,", sd_int,")"), class="Intercept", dpar=dpar))
  }))
}

extract_biomarkers_brms <- function(fit, top_k=50){
  fx <- as.data.frame(brms::fixef(fit, robust=TRUE)); fx$term <- rownames(fx)
  fx <- fx[!grepl("Intercept", fx$term), ]
  fx$feature <- sub("^b_", "", sub("^.*?b_", "b_", fx$term))
  fx$class   <- sub("^mu([^_]+).*", "\\1", fx$term)
  fx$nonzero <- (fx$Q2.5 > 0) | (fx$Q97.5 < 0)
  agg <- fx %>% group_by(feature) %>% summarise(any_nonzero=any(nonzero),
                                                max_abs_est=max(abs(Estimate), na.rm=TRUE),
                                                n_classes=sum(nonzero), .groups="drop") %>%
    arrange(desc(any_nonzero), desc(n_classes), desc(max_abs_est))
  list(per_class = fx[fx$nonzero, c("class","feature","Estimate","Q2.5","Q97.5")],
       ranking   = head(agg, top_k))
}

# ==========================
# (1) Cargar datos y armonizar columnas
# ==========================

df_train <- readRDS("./df_train_final.rds")
df_test  <- readRDS("./df_test_final.rds")

ref <- "NP"
df_train$y <- stats::relevel(refactor_levels(df_train$y), ref = ref)
df_test$y  <- factor(df_test$y, levels = levels(df_train$y))

common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
new_names <- make.unique(sanitize_names(common_cols), sep="_")
colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
colnames(df_test)[match(common_cols,  colnames(df_test))]  <- new_names

df_train <- df_train[, c("y", new_names), drop=FALSE]
df_test  <- df_test[,  c("y", new_names), drop=FALSE]

# ==========================
# (2) LASSO multinomial para preselección (con balanceo interno)
# ==========================
# ==========================
# (2) LASSO multinomial para preselección (relajado)
# ==========================
set.seed(123)

# Balanceo solo para glmnet
class_counts <- table(df_train$y)
df_train_bal <- if(any(class_counts < max(class_counts))) {
  caret::upSample(df_train[,-1], df_train$y, yname="y")
} else df_train

x_tr <- as.matrix(df_train_bal[, setdiff(names(df_train_bal), "y")])
y_tr <- factor(df_train_bal$y, levels=levels(df_train$y))

nfolds <- max(3, min(5, length(y_tr)))
cvfit <- cv.glmnet(x_tr, y_tr, family="multinomial",
                   type.multinomial="grouped",
                   alpha=1, standardize=TRUE, nfolds=nfolds)

get_nonzero <- function(fit, lambda){
  coefs <- coef(fit, s=lambda)  # lista por clase
  feats <- unique(unlist(lapply(coefs, function(m){
    nm <- rownames(m); if(is.null(nm) || nrow(m)<=1) return(character(0))
    nz <- as.vector(m[-1, , drop=FALSE]) != 0
    nm <- nm[-1]; nm[nz]
  })))
  feats[feats %in% colnames(x_tr)]
}

# --- Selección relajada: buscamos el primer lambda (más permisivo) que alcance min_k
relaxed_select <- function(fit, min_k=5, max_k=50){
  lambdas <- sort(fit$lambda, decreasing=FALSE)  # de fuerte a débil penalización
  for(l in lambdas){
    nz <- get_nonzero(fit, l)
    if(length(nz) >= min_k) return(head(nz, max_k))
  }
  # Fallbacks
  nz_min <- get_nonzero(fit, min(lambdas))
  if(length(nz_min)) return(head(nz_min, max_k))
  nz_1se <- get_nonzero(fit, cvfit$lambda.1se)
  if(length(nz_1se)) return(head(nz_1se, max_k))
  nz_min2 <- get_nonzero(fit, cvfit$lambda.min)
  if(length(nz_min2)) return(head(nz_min2, max_k))
  character(0)
}

sel <- relaxed_select(cvfit$glmnet.fit, min_k=3, max_k=100)
message("# Features seleccionados (relajado): ", length(sel))
print(sel)

# Reducimos datasets con lo seleccionado (automático)
keep <- unique(c("y", sel))
# df_train_sel <- df_train[, keep, drop=FALSE]
# df_test_sel  <- df_test[,  keep, drop=FALSE]

df_train_sel <- df_train
df_test_sel  <- df_test

# ==========================
# (3) BRMS final sobre el set reducido
# ==========================

priors <- make_multinom_priors(levels(df_train_sel$y), ref = ref, sd_b = 1, sd_int = 5)
ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

fit_brm <- brms::brm(
  y ~ ., data = df_train_sel, family = brms::categorical(), prior = priors,
  chains = 4, iter = 8000, warmup = 4000, refresh = 0, control = ctrl
)

cat("\n=== Resumen BRMS (set reducido) ===\n")
print(summary(fit_brm))
cat("\n=== LOO (moment_match=TRUE) ===\n"); print(loo::loo(fit_brm, moment_match=TRUE))
cat("\n=== WAIC ===\n"); print(loo::waic(fit_brm))

# ==========================
# (4) Métricas predictivas en test
# ==========================

P <- brms_probs(fit_brm, df_test_sel)[, levels(df_test_sel$y), drop=FALSE]
cat("\n=== Métricas TEST ===\n")
print(metrics_from_probs(P, df_test_sel$y))
cat("\n=== Matriz de confusión ===\n")
print(caret::confusionMatrix(
  factor(colnames(P)[max.col(P, ties.method="first")], levels=levels(df_test_sel$y)),
  df_test_sel$y
))

# ==========================
# (5) Biomarcadores robustos (BRMS)
# ==========================

biomarkers <- extract_biomarkers_brms(fit_brm, top_k=50)
cat("\n=== Coeficientes con IC que excluye 0 (por clase) ===\n"); print(biomarkers$per_class)
cat("\n=== Ranking agregado (prioriza IC≠0, nº clases y |beta|) ===\n"); print(biomarkers$ranking)

# ==========================
# (6) Estabilidad de selección vía bootstrap (GLMNET)
#     (solo en TRAIN para no tocar TEST) — con re-balanceo en cada resample
# ==========================

bootstrap_lasso <- function(df, B=50){
  out <- vector("list", B)
  levs <- levels(df$y)
  for(b in seq_len(B)){
    idx <- sample(seq_len(nrow(df)), size=nrow(df), replace=TRUE)
    dfb <- df[idx, , drop=FALSE]
    # rebalancear dentro del bootstrap si hace falta
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

cat("\n=== Bootstrap de estabilidad (GLMNET) ===\n")
stab_tbl <- bootstrap_lasso(df_train_sel, B=50)
print(head(stab_tbl, 20))

# Candidatos estables (>=40% de las corridas)
stable_candidates <- stab_tbl %>% filter(freq >= 0.2*50) %>% pull(feature)
cat("\n# Candidatos estables (>=40%): ", length(stable_candidates), "\n"); print(stable_candidates)

# ==========================
# (7) (Opcional) Reajuste BRMS con sólo candidatos estables
# ==========================

if(length(stable_candidates) >= 1){
  keep2 <- c("y", intersect(stable_candidates, colnames(df_train_sel)))
  df_train_stab <- df_train_sel[, keep2, drop=FALSE]
  df_test_stab  <- df_test_sel[,  keep2, drop=FALSE]
  pri2 <- make_multinom_priors(levels(df_train_stab$y), ref = ref, sd_b = 1, sd_int = 5)
  fit_brm_stab <- brms::brm(
    y ~ ., data = df_train_stab, family = brms::categorical(), prior = pri2,
    chains = 4, iter = 6000, warmup = 3000, refresh = 0, control = ctrl
  )
  cat("\n=== Métricas TEST (BRMS con features estables) ===\n")
  P2 <- brms_probs(fit_brm_stab, df_test_stab)[, levels(df_test_stab$y), drop=FALSE]
  print(metrics_from_probs(P2, df_test_stab$y))
  print(caret::confusionMatrix(
    factor(colnames(P2)[max.col(P2, ties.method="first")], levels=levels(df_test_stab$y)),
    df_test_stab$y
  ))
}



##===aqui====

# ==========================
# MODELOS CLÁSICOS CON FEATURES SELECCIONADAS
# ==========================

suppressPackageStartupMessages({
  library(dplyr); library(caret); library(pROC)
  library(glmnet); library(klaR); library(e1071); library(randomForest)
})

set.seed(123)

# Dataset reducido a features candidatos
# Dataset reducido: usa las seleccionadas por LASSO (sel) automáticamente
keep <- unique(c("y", sel))
df_train_sel <- df_train[, keep, drop=FALSE]
df_test_sel  <- df_test[,  keep, drop=FALSE]


# Métricas personalizadas (incluye logLoss y AUC multiclass)
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
  
  out <- c(Accuracy=acc, BalAcc=bal, Kappa=kap, AUC=auc, LogLoss=logloss)
  out
}

ctrl <- trainControl(
  method="cv", number=3,  # CV interna
  classProbs=TRUE,
  summaryFunction=multiSummary,
  savePredictions="final"
)

# --------------------------
# Modelos
# --------------------------
set.seed(123)

# Ridge (alpha=0)
fit_ridge <- train(y ~ ., data=df_train_sel,
                   method="glmnet",
                   trControl=ctrl,
                   tuneGrid=expand.grid(alpha=0, lambda=10^seq(-3,1,length=20)))

# Lasso (alpha=1)
fit_lasso <- train(y ~ ., data=df_train_sel,
                   method="glmnet",
                   trControl=ctrl,
                   tuneGrid=expand.grid(alpha=1, lambda=10^seq(-3,1,length=20)))

# Elastic Net (alpha entre 0 y 1)
fit_enet <- train(y ~ ., data=df_train_sel,
                  method="glmnet",
                  trControl=ctrl,
                  tuneLength=10)  # caret explora alpha y lambda

# SVM radial
fit_svm <- train(y ~ ., data=df_train_sel,
                 method="svmRadial",
                 trControl=ctrl,
                 tuneLength=5)

# Naive Bayes
fit_nb <- train(y ~ ., data=df_train_sel,
                method="nb",
                trControl=ctrl,
                tuneLength=5)

# Random Forest
fit_rf <- train(y ~ ., data=df_train_sel,
                method="rf",
                trControl=ctrl,
                tuneLength=5)

# --------------------------
# Evaluación en TEST
# --------------------------
model_list <- list(Ridge=fit_ridge, Lasso=fit_lasso, ENet=fit_enet,
                   SVM=fit_svm, NB=fit_nb, RF=fit_rf)

test_results <- lapply(model_list, function(fit){
  probs <- predict(fit, newdata=df_test_sel, type="prob")
  metrics_from_probs(as.matrix(probs), df_test_sel$y)
})

test_results <- bind_rows(test_results, .id="Model")
print(test_results)



# ==========================
# IMPORTANCIA DE VARIABLES (Ridge, ENet, RF)
# ==========================

# --- Ridge / Elastic Net ---
coef_importance <- function(fit){
  if(!inherits(fit$finalModel, "glmnet") && !inherits(fit$finalModel, "cv.glmnet"))
    return(NULL)
  
  # Usar el mejor lambda que eligió caret
  lambda <- fit$bestTune$lambda
  
  coefs <- try(coef(fit$finalModel, s=lambda), silent=TRUE)
  if(inherits(coefs, "try-error")) {
    # fallback: tomar coeficientes sin s
    coefs <- coef(fit$finalModel)
  }
  
  if(is.list(coefs)){
    # multinomial: lista por clase
    v <- Reduce(`+`, lapply(coefs, function(m){
      a <- abs(as.numeric(m[-1, , drop=FALSE]))
      names(a) <- rownames(m)[-1]
      a
    }))
  } else {
    v <- abs(as.numeric(coefs[-1]))
    names(v) <- rownames(coefs)[-1]
  }
  
  v <- sort(v, decreasing=TRUE)
  return(v)
}


imp_ridge <- coef_importance(fit_ridge)
imp_enet  <- coef_importance(fit_enet)
imp_rf    <- varImp(fit_rf, scale=TRUE)

cat("\n=== Importancia Ridge (abs coef) ===\n")
print(head(imp_ridge, 20))

cat("\n=== Importancia ENet (abs coef) ===\n")
print(head(imp_enet, 20))

cat("\n=== Importancia RF ===\n")
print(imp_rf)

print(varImp(fit_rf, scale=TRUE))
# --- Consolidar ranking ---
# --- Consolidar ranking ---
all_features <- unique(c(names(imp_ridge), names(imp_enet), rownames(imp_rf$importance)))
importance_tbl <- data.frame(
  Feature = all_features,
  Ridge   = imp_ridge[all_features],
  ENet    = imp_enet[all_features],
  RF      = imp_rf$importance[all_features, 1]
)

importance_tbl <- importance_tbl %>%
  mutate(across(-Feature, ~replace_na(., 0))) %>%
  mutate(AverageRank = (rank(-Ridge, ties.method="average") +
                          rank(-ENet, ties.method="average") +
                          rank(-RF, ties.method="average"))/3) %>%
  arrange(AverageRank)

cat("\n=== Ranking consolidado de importancia ===\n")
print(importance_tbl)


