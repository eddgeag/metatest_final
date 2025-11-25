

rm(list = ls())

##===librerias =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(glmnet); library(nnet); library(MASS)
  library(caret); library(pROC); library(FNN)
  library(ggplot2); library(ggridges)
  library(MOFA2); library(AnnotationDbi);library(brms)
})

set.seed(123)

##===Funciones =====

refactor_levels <- function(x) factor(as.character(x))

sanitize_names <- function(v){
  v2 <- gsub("[^A-Za-z0-9]+", "_", v)
  v2 <- gsub("__+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- ifelse(grepl("^[0-9]", v2), paste0("X", v2), v2)
  v2 <- gsub("_+$", "", v2)
  v2 <- make.names(v2, unique = FALSE)
  v2 <- gsub("\\.+", "_", v2)
  v2 <- gsub("^_+|_+$", "", v2)
  v2 <- gsub("__+", "_", v2)
  make.unique(v2, sep = "_")
}

std_and_prune <- function(df){
  y <- df$y
  X <- df[, setdiff(colnames(df),"y"), drop=FALSE]
  nzv <- caret::nearZeroVar(X)
  if(length(nzv)) X <- X[, -nzv, drop=FALSE]
  if(ncol(X) > 2){
    cc <- stats::cor(X, use="pairwise.complete.obs")
    bad <- suppressWarnings(caret::findCorrelation(cc, cutoff=0.95, names=FALSE))
    if(length(bad)) X <- X[, -bad, drop=FALSE]
  }
  mu <- vapply(X, mean, 0.0, na.rm=TRUE)
  sd <- vapply(X,   sd, 0.0, na.rm=TRUE)
  sd[sd == 0] <- 1
  Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
  data.frame(y=y, Xs, check.names=FALSE)
}


select_features_view <- function(mod, view, q_thr, data_list, metadata = NULL){
  aux <- plot_top_weights(mod, view=view, factors=1:2, abs=TRUE, scale=TRUE, nfeatures=Inf)
  thr <- stats::quantile(aux$data$value, probs=q_thr, na.rm=TRUE)
  feats <- aux$data$feature[aux$data$value >= thr]
  sub_train <- data_list$train[as.character(feats), , drop=FALSE]
  sub_test  <- data_list$test[as.character(feats),  , drop=FALSE]
  if(!is.null(metadata) && all(c("EntrezGeneID","GeneSymbol") %in% colnames(metadata))){
    map_dic <- metadata %>% dplyr::select(EntrezGeneID, GeneSymbol) %>% distinct()
    idx <- rownames(sub_train)
    ids_map <- data.frame(EntrezGeneID = idx) %>% dplyr::left_join(map_dic, by = "EntrezGeneID")
    new_ids <- ifelse(!is.na(ids_map$GeneSymbol), ids_map$GeneSymbol, ids_map$EntrezGeneID)
    rownames(sub_train) <- new_ids; rownames(sub_test) <- new_ids
  }
  list(train=sub_train, test=sub_test)
}

project_all_views <- function(X_list, W_list){
  W_all <- NULL; X_all <- NULL
  for(v in names(W_list)){
    Xv <- X_list[[v]]; if(is.list(Xv) && "test" %in% names(Xv)) Xv <- Xv$test
    common_feats <- intersect(rownames(W_list[[v]]), rownames(Xv))
    if(length(common_feats) == 0) next
    W_sub <- as.matrix(W_list[[v]][common_feats, , drop=FALSE])
    X_sub <- as.matrix(Xv[common_feats, , drop=FALSE])
    W_all <- rbind(W_all, W_sub); X_all <- rbind(X_all, X_sub)
  }
  stopifnot(nrow(W_all) == nrow(X_all))
  Z_est <- solve(t(W_all) %*% W_all) %*% t(W_all) %*% X_all
  t(Z_est)
}

build_sets_from_grid <- function(thr_grid, tst_Data, pesos_list, Z_seed, y_seed, mod){
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
    tx <- select_features_view(mod, 1, thr_grid$q_tx[i], tst_Data$tx, tst_Data$features_metadata)
    pr <- select_features_view(mod, 2, thr_grid$q_pr[i], tst_Data$pr, tst_Data$features_metadata)
    me <- select_features_view(mod, 3, thr_grid$q_me[i], tst_Data$me, tst_Data$features_metadata)
    cl <- select_features_view(mod, 4, thr_grid$q_cl[i], tst_Data$cl, tst_Data$features_metadata)
    train_feats <- rbind(tx$train, pr$train, me$train, cl$train)
    test_feats  <- rbind(tx$test,  pr$test,  me$test,  cl$test)
    out_train[[i]] <- cbind(Z_seed, t(train_feats)) %>% as.data.frame()
    out_test[[i]]  <- cbind(factores_Test, t(test_feats)) %>% as.data.frame()
  }
  list(train=out_train, test=out_test, y_train=y_seed, y_test=tst_Data$grupo_test$grupo)
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
         BalAcc  =mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE),
         Kappa   =as.numeric(cm$overall["Kappa"]),
         AUC=auc_num, LogLoss=logloss)
}

make_multinom_priors <- function(y_levels, ref="NP", sd_b=0.3, sd_int=3){
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(set_prior(paste0("student_t(3,0,", sd_b, ")"), class="b", dpar=dpar),
      set_prior(paste0("student_t(3,0,", sd_int,")"), class="Intercept", dpar=dpar))
  }))
}

brms_probs <- function(fit, newdata){
  pp <- posterior_epred(fit, newdata=newdata)
  M  <- apply(pp, c(2,3), mean)
  cls <- dimnames(pp)[[3]]; if(is.null(cls)) cls <- fit$family$names
  colnames(M) <- cls; M
}

count_coefs_nnet <- function(fit){
  cf <- coef(fit); if(is.list(cf)) cf <- do.call(rbind, cf)
  cf <- as.matrix(cf)
  if(ncol(cf) > 0) cf_noint <- cf[, -1, drop=FALSE] else cf_noint <- cf
  sum(abs(cf_noint) > 0)
}

count_coefs_brms <- function(fit){
  fe <- tryCatch(as.data.frame(brms::fixef(fit)), error=function(e) NULL)
  if(is.null(fe)) return(NA_integer_)
  keep <- !grepl("Intercept", rownames(fe))
  sum(abs(fe$Estimate[keep]) > 0)
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


evaluate_all <- function(return.list, thr_grid, ref="NP",
                         do_brms=FALSE, brms_iter=2000, brms_chains=2,
                         brms_algorithm=c("mcmc","vb")){
  options(contrasts = c("contr.treatment","contr.poly"))
  brms_algorithm <- match.arg(brms_algorithm)
  res_all <- list()
  
  detect_backend <- function(){
    if (requireNamespace("cmdstanr", quietly=TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent=TRUE)
      if (!inherits(v, "try-error")) return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  rank_one <- function(df){
    df %>% mutate(r1 = rank(-BalAcc, ties.method="min"),
                  r2 = rank(LogLoss,  ties.method="min"),
                  r3 = rank(-AUC,    ties.method="min"),
                  RankSum = r1 + r2 + r3) %>%
      arrange(RankSum, desc(BalAcc), AUC, LogLoss)
  }
  
  # construye control según backend
  make_ctrl <- function(backend){
    if(identical(backend, "cmdstanr")){
      list(adapt_delta=0.95, max_treedepth=12)                 # sin init_r
    } else {
      list(adapt_delta=0.95, max_treedepth=12, init_r=0.1)     # rstan permite init_r
    }
  }
  
  fit_brms_with_fallback <- function(df_train, ref_lbl, algo=brms_algorithm,
                                     iter=brms_iter, chains=brms_chains, backend=backend_choice){
    priors <- make_multinom_priors(levels(df_train$y), ref=ref_lbl, sd_b=0.3, sd_int=3)
    family_cat <- brms::categorical()
    vb_err <- NULL; mcmc_err <- NULL
    
    if(algo %in% c("vb","meanfield","fullrank")){
      fit_vb <- tryCatch(
        brm(y ~ ., data=df_train, family=family_cat, prior=priors,
            algorithm="fullrank", iter=max(4000, iter), refresh=0, backend=backend,
            control=list(tol_rel_obj=1e-4, eval_elbo=200, adapt_engaged=TRUE)),
        error=function(e){ vb_err <<- conditionMessage(e); NULL }
      )
      if(!is.null(fit_vb)){
        kk <- tryCatch({ ll <- brms::loo(fit_vb); max(ll$diagnostics$pareto_k) }, error=function(e) Inf)
        if(is.finite(kk) && kk <= 0.7) return(list(fit=fit_vb, label="BRMS_VB", err=NULL))
        vb_err <- sprintf("VB Pareto-k=%.2f", kk)
      }
      if(!is.null(vb_err)) message("[VB] ", vb_err)
    }
    
    fit_mcmc <- tryCatch(
      brm(y ~ ., data=df_train, family=family_cat, prior=priors,
          chains=max(2, chains), iter=max(2000, iter), warmup=floor(max(2000, iter)/2),
          refresh=0, backend=backend, control=make_ctrl(backend)),
      error=function(e){ mcmc_err <<- conditionMessage(e); NULL }
    )
    if(!is.null(fit_mcmc)) return(list(fit=fit_mcmc, label="BRMS_MCMC", err=NULL))
    
    errs <- c(vb_err, mcmc_err); errs <- errs[!vapply(errs, is.null, logical(1))]
    if(length(errs)==0) errs <- NA_character_
    list(fit=NULL, label=NA_character_, err=paste(unlist(errs), collapse=" | "))
  }
  
  for(i in seq_along(return.list$train)){
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref=ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test),  ref=ref)
    
    df_train <- data.frame(y=y_train, X_train, check.names=FALSE)
    df_test  <- data.frame(y=y_test,  X_test,  check.names=FALSE)
    
    df_train$y <- droplevels(df_train$y)
    df_test$y  <- factor(df_test$y, levels=levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- sanitize_names(common_cols); new_names <- make.unique(new_names, sep="_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols,  colnames(df_test))]  <- new_names
    
    df_train <- df_train[, c("y", new_names), drop=FALSE]
    df_test  <- df_test[,  c("y", new_names), drop=FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    fit_nnet <- nnet::multinom(y ~ ., data=df_train, trace=FALSE, MaxNWts=200000)
    probs_nnet <- stats::predict(fit_nnet, newdata=df_test, type="probs")
    if(is.null(colnames(probs_nnet))) colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% mutate(Model="NNET", GridID=i, p=p, n_coefs=count_coefs_nnet(fit_nnet))
    
    if(do_brms){
      fr <- fit_brms_with_fallback(df_train, ref_lbl=ref, algo=brms_algorithm,
                                   iter=brms_iter, chains=brms_chains, backend=backend_choice)
      if(!is.null(fr$fit)){
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop=FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% mutate(Model=fr$label, GridID=i, p=p, n_coefs=count_coefs_brms(fr$fit))
        attr(rows[[fr$label]], "biomarkers") <- extract_biomarkers_brms(fr$fit, top_k=100)
      } else if(!is.null(fr$err)){
        brms_errs[[length(brms_errs)+1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy,BalAcc,Kappa,AUC,LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% mutate(GridID = dplyr::row_number()), by="GridID")
  
  rank_one <- rank_one
  best_by_grid <- res_tbl %>% group_by(GridID) %>% group_modify(~ rank_one(.x) %>% slice_head(n=1)) %>% ungroup()
  best_overall <- rank_one(res_tbl) %>% slice_head(n=1)
  
  list(results=res_tbl, best_by_grid=best_by_grid, best_overall=best_overall,
       brms_errors=unique(unlist(brms_errs)), backend=backend_choice)
}

##===ejecucion=====

dats       <- readRDS("./datos_para_modelar.rds")
tst_Data   <- readRDS("./ready_for_modeling.rds")
factors    <- dats$factores
pesos_list <- dats$pesos
Z_seed <- as.matrix(factors %>% dplyr::filter(Tipo=="Real") %>% dplyr::select(where(is.numeric)))[,1:2, drop=FALSE]
y_seed <- factors %>% dplyr::filter(Tipo=="Real") %>% dplyr::pull(Grupo) %>% factor(); names(y_seed) <- rownames(Z_seed)
mod <- readRDS("./modelo.rds")

##===grilla de thresholds=====
thr_grid <- expand.grid(
  q_tx = c(0.90, 0.95,0.97),
  q_pr = c(0.70, 0.85),
  q_me = c(0.70, 0.85,0.95,0.97),
  q_cl = c(0.65,0.70)
)

return.list <- build_sets_from_grid(thr_grid, tst_Data, pesos_list, Z_seed, y_seed, mod)
detach("package:MOFA2", unload=TRUE)
EVAL <- evaluate_all(return.list, thr_grid, ref="NP", do_brms=TRUE,
                     brms_algorithm="mcmc", brms_iter=4000, brms_chains=4)



evaluacion <- EVAL$best_by_grid
evaluacion <- evaluacion[evaluacion$Accuracy<1,]
grid_id <- evaluacion$GridID[which.max(evaluacion$BalAcc)]
thr_grid



# ==== BLOQUE: evaluación detallada Grid 8 (métricas, ROC multiclase, matriz de confusión) ====
grid_id <- 33L ## fue la mejor por log loss
ref <- "NP"
backend <- EVAL$backend

# Datos consistentes con evaluate_all
X_train <- return.list$train[[grid_id]]
X_test  <- return.list$test[[grid_id]]
y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
y_test  <- stats::relevel(refactor_levels(return.list$y_test),  ref = ref)

df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
df_test  <- data.frame(y = y_test,  X_test,  check.names = FALSE)
df_train$y <- droplevels(df_train$y)
df_test$y  <- factor(df_test$y, levels = levels(df_train$y))

common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
new_names <- make.unique(sanitize_names(common_cols), sep = "_")
colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
colnames(df_test)[match(common_cols,  colnames(df_test))]  <- new_names
df_train <- df_train[, c("y", new_names), drop = FALSE]
df_test  <- df_test[,  c("y", new_names), drop = FALSE]

# BRMS_MCMC (control según backend)
priors <- make_multinom_priors(levels(df_train$y), ref = ref, sd_b = 0.3, sd_int = 3)
ctrl <- if (identical(backend, "cmdstanr")){ list(adapt_delta = 0.95, max_treedepth = 12)
}else {list(adapt_delta = 0.95, max_treedepth = 12, init_r = 0.1)}

fit_brm <- brms::brm(
  y ~ ., data = df_train, family = brms::categorical(), prior = priors,
  chains = 4, iter = 8000, warmup = 4000, refresh = 0, backend = backend, control = ctrl
)






saveRDS(df_train,"./df_train_final.rds")
saveRDS(df_test,"./df_test_final.rds")

# Probabilidades y métricas
P <- brms_probs(fit_brm, df_test)[, levels(df_test$y), drop = FALSE]
metrics_grid8 <- metrics_from_probs(P, df_test$y)

# Matriz de confusión
pred_lab <- factor(colnames(P)[max.col(P, ties.method = "first")], levels = levels(df_test$y))
cm_grid8 <- caret::confusionMatrix(pred_lab, df_test$y)


# === Diagnóstico numérico del ajuste BRMS====

# Resumen de parámetros (coeficientes, errores, Rhat, ESS)
summary(fit_brm)

# Medidas de convergencia
diag <- list(
  rhat   = max(summary(fit_brm)$fixed[,"Rhat"], na.rm=TRUE),
  ess_bulk = min(summary(fit_brm)$fixed[,"Bulk_ESS"], na.rm=TRUE),
  ess_tail = min(summary(fit_brm)$fixed[,"Tail_ESS"], na.rm=TRUE)
)
print(diag)

# Log-likelihood y WAIC/LOO
loo_brm  <- loo::loo(fit_brm, moment_match = TRUE)
waic_brm <- loo::waic(fit_brm)

print(loo_brm)
print(waic_brm)

# Medidas de ajuste predictivo
P <- brms_probs(fit_brm, df_test)[, levels(df_test$y), drop = FALSE]
metrics_grid8 <- metrics_from_probs(P, df_test$y)
print(metrics_grid8)




# ROC multiclase (one-vs-all + macro AUC)
lev <- levels(df_test$y)
roc_list <- lapply(lev, function(lv) pROC::roc(as.integer(df_test$y == lv), as.numeric(P[, lv]), quiet = TRUE))
names(roc_list) <- lev
mroc <- pROC::multiclass.roc(df_test$y, P)

# Curvas ROC por clase con ggplot
roc_df <- do.call(rbind, lapply(names(roc_list), function(lv){
  r <- roc_list[[lv]]
  data.frame(specificity = rev(r$specificities),
             sensitivity = rev(r$sensitivities),
             class = lv)
}))
gg_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = class)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  labs(title = sprintf("ROC one-vs-all (macro AUC = %.3f)", as.numeric(mroc$auc)),
       x = "1 - Especificidad", y = "Sensibilidad") +
  theme_minimal()

# Salidas
print(metrics_grid8)
print(cm_grid8)
print(gg_roc)


























