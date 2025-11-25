# =======================================================
# Pipeline MOFA -> selección de features -> modelos NNET/BRMS
# Limpio, con orden correcto de funciones y fix de 'vb_err not found'
# =======================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
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

# -------------------------------------------------------
# 1) UTILIDADES MOFA: selección de features por vista
# -------------------------------------------------------
select_features_view <- function(mod, view, q_thr, data_list, metadata = NULL) {
  aux <- plot_top_weights(
    mod,
    view = view,
    factors = 1:2,
    abs = TRUE,
    scale = TRUE,
    nfeatures = Inf
  )
  thr <- stats::quantile(aux$data$value, probs = q_thr, na.rm = TRUE)
  feats <- aux$data$feature[aux$data$value >= thr]
  
  sub_train <- data_list$train[as.character(feats), , drop = FALSE]
  sub_test  <- data_list$test[as.character(feats), , drop = FALSE]
  
  if (!is.null(metadata) && all(c("EntrezGeneID", "GeneSymbol") %in% colnames(metadata))) {
    map_dic <- metadata %>% dplyr::select(EntrezGeneID, GeneSymbol) %>% distinct()
    idx <- rownames(sub_train)
    ids_map <- data.frame(EntrezGeneID = idx) %>% dplyr::left_join(map_dic, by = "EntrezGeneID")
    new_ids <- ifelse(!is.na(ids_map$GeneSymbol), ids_map$GeneSymbol, ids_map$EntrezGeneID)
    rownames(sub_train) <- new_ids
    rownames(sub_test)  <- new_ids
  }
  list(train = sub_train, test = sub_test)
}

# -------------------------------------------------------
# 2) Proyección Z(test) usando pesos W por vista
# -------------------------------------------------------
project_all_views <- function(X_list, W_list) {
  W_all <- NULL; X_all <- NULL
  for (v in names(W_list)) {
    Xv <- X_list[[v]]
    if (is.list(Xv) && "test" %in% names(Xv)) Xv <- Xv$test
    common_feats <- intersect(rownames(W_list[[v]]), rownames(Xv))
    if (!length(common_feats)) next
    W_sub <- as.matrix(W_list[[v]][common_feats, , drop = FALSE])
    X_sub <- as.matrix(Xv[common_feats, , drop = FALSE])
    W_all <- rbind(W_all, W_sub)
    X_all <- rbind(X_all, X_sub)
  }
  stopifnot(nrow(W_all) == nrow(X_all))
  Z_est <- solve(t(W_all) %*% W_all) %*% t(W_all) %*% X_all
  t(Z_est)
}

refactor_levels <- function(x) as.factor(as.character(x))

# -------------------------------------------------------
# 3) Construcción de sets por grid de umbrales
# -------------------------------------------------------
build_sets_from_grid <- function(thr_grid, tst_Data, pesos_list, Z_seed, y_seed, mod) {
  out_train <- vector("list", nrow(thr_grid))
  out_test  <- vector("list", nrow(thr_grid))
  
  X_list <- list(
    transcriptomica = tst_Data$tx$test,
    proteomica      = tst_Data$pr$test,
    metabolomica    = tst_Data$me$test,
    clinical        = tst_Data$cl$test
  )
  factores_Test <- project_all_views(X_list, pesos_list)
  
  for (i in seq_len(nrow(thr_grid))) {
    tx <- select_features_view(mod, 1, thr_grid$q_tx[i], tst_Data$tx, tst_Data$features_metadata)
    pr <- select_features_view(mod, 2, thr_grid$q_pr[i], tst_Data$pr, tst_Data$features_metadata)
    me <- select_features_view(mod, 3, thr_grid$q_me[i], tst_Data$me, tst_Data$features_metadata)
    cl <- select_features_view(mod, 4, thr_grid$q_cl[i], tst_Data$cl, tst_Data$features_metadata)
    
    train_feats   <- rbind(tx$train, pr$train, me$train, cl$train)
    test_feats    <- rbind(tx$test,  pr$test,  me$test,  cl$test)
    
    out_train[[i]] <- cbind(Z_seed, as.data.frame(t(train_feats)))
    out_test[[i]]  <- cbind(factores_Test, as.data.frame(t(test_feats)))
  }
  
  list(train = out_train, test = out_test, y_train = y_seed, y_test = tst_Data$grupo_test$grupo)
}

# -------------------------------------------------------
# 4) Métricas, helpers y limpieza
# -------------------------------------------------------
metrics_from_probs <- function(prob_mat, y_true) {
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
  tibble(Accuracy = as.numeric(cm$overall["Accuracy"]),
         BalAcc   = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
         Kappa    = as.numeric(cm$overall["Kappa"]),
         AUC      = auc_num,
         LogLoss  = logloss)
}

make_multinom_priors <- function(y_levels, ref = "NP", sd_b = 0.3, sd_int = 3) {
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(
      brms::set_prior(paste0("student_t(3,0,", sd_b, ")"), class = "b", dpar = dpar),
      brms::set_prior(paste0("student_t(3,0,", sd_int,")"), class = "Intercept", dpar = dpar)
    )
  }))
}

brms_probs <- function(fit, newdata) {
  pp <- brms::posterior_epred(fit, newdata = newdata)
  M  <- apply(pp, c(2, 3), mean)
  cls <- dimnames(pp)[[3]]; if (is.null(cls)) cls <- fit$family$names
  colnames(M) <- cls
  M
}

count_coefs_nnet <- function(fit) {
  cf <- coef(fit); if (is.list(cf)) cf <- do.call(rbind, cf)
  cf <- as.matrix(cf)
  if (ncol(cf) > 0) cf <- cf[, -1, drop = FALSE]
  sum(abs(cf) > 0)
}

count_coefs_brms <- function(fit) {
  fe <- tryCatch(as.data.frame(brms::fixef(fit)), error = function(e) NULL)
  if (is.null(fe)) return(NA_integer_)
  keep <- !grepl("Intercept", rownames(fe))
  sum(abs(fe$Estimate[keep]) > 0)
}

std_and_prune <- function(df){
  y <- df$y
  X <- df[, setdiff(colnames(df),"y"), drop=FALSE]
  nzv <- caret::nearZeroVar(X)
  if(length(nzv)) X <- X[, -nzv, drop=FALSE]
  if(ncol(X) > 2){
    cc <- stats::cor(X, use="pairwise.complete.obs")
    bad <- caret::findCorrelation(cc, cutoff=0.95, names=FALSE)
    if(length(bad)) X <- X[, -bad, drop=FALSE]
  }
  mu <- vapply(X, mean, 0.0, na.rm=TRUE)
  sd <- vapply(X,   sd, 0.0, na.rm=TRUE); sd[sd == 0] <- 1
  Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
  data.frame(y=y, Xs, check.names=FALSE)
}

extract_biomarkers_brms <- function(fit, top_k=50){
  fx <- as.data.frame(brms::fixef(fit, robust=TRUE)); fx$term <- rownames(fx)
  fx <- fx[!grepl("Intercept", fx$term), ]
  fx$feature <- sub("^b_", "", sub("^.*?b_", "b_", fx$term))
  fx$class   <- sub("^mu([^_]+).*", "\\1", fx$term)
  fx$nonzero <- (fx$Q2.5 > 0) | (fx$Q97.5 < 0)
  agg <- fx %>% dplyr::group_by(feature) %>% dplyr::summarise(
    any_nonzero = any(nonzero),
    max_abs_est = max(abs(Estimate), na.rm=TRUE),
    n_classes   = sum(nonzero), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(any_nonzero), dplyr::desc(n_classes), dplyr::desc(max_abs_est))
  list(per_class = fx[fx$nonzero, c("class","feature","Estimate","Q2.5","Q97.5")],
       ranking   = head(agg, top_k))
}

# -------------------------------------------------------
# 5) Núcleo: evaluate_all (única definición)
#    Fix: inicializar vb_err/mcmc_err para evitar 'object not found'
# -------------------------------------------------------
evaluate_all <- function(return.list,
                         thr_grid,
                         ref = "NP",
                         do_brms = FALSE,
                         brms_iter = 2000,
                         brms_chains = 2,
                         brms_algorithm = c("mcmc", "vb")) {
  options(contrasts = c("contr.treatment", "contr.poly"))
  brms_algorithm <- match.arg(brms_algorithm)
  res_all <- list()
  
  detect_backend <- function() {
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
      if (!inherits(v, "try-error")) return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  sanitize_names <- function(v) {
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
  
  rank_one <- function(df) {
    df %>% dplyr::mutate(
      r1 = rank(-BalAcc, ties.method = "min"),
      r2 = rank(LogLoss, ties.method = "min"),
      r3 = rank(-AUC,    ties.method = "min"),
      RankSum = r1 + r2 + r3
    ) %>% dplyr::arrange(RankSum, dplyr::desc(BalAcc), AUC, LogLoss)
  }
  
  fit_brms_with_fallback <- function(df_train, ref_lbl,
                                     algo = brms_algorithm,
                                     iter = brms_iter,
                                     chains = brms_chains,
                                     backend = backend_choice) {
    priors <- make_multinom_priors(levels(df_train$y), ref = ref_lbl, sd_b = 0.3, sd_int = 3)
    family_cat <- brms::categorical()
    
    vb_err   <- NA_character_
    mcmc_err <- NA_character_
    
    if (algo %in% c("vb", "meanfield", "fullrank")) {
      fit_vb <- tryCatch(
        brms::brm(
          y ~ ., data = df_train, family = family_cat, prior = priors,
          algorithm = "fullrank", iter = max(4000, iter), refresh = 0,
          backend = backend,
          control = list(tol_rel_obj = 1e-4, eval_elbo = 200, adapt_engaged = TRUE)
        ),
        error = function(e) { vb_err <<- conditionMessage(e); NULL }
      )
      if (!is.null(fit_vb)) {
        kk <- tryCatch({ ll <- brms::loo(fit_vb); max(ll$diagnostics$pareto_k) }, error = function(e) Inf)
        if (is.finite(kk) && kk <= 0.7)
          return(list(fit = fit_vb, label = "BRMS_VB", err = NULL))
        vb_err <- sprintf("VB Pareto-k=%.2f", kk)
      }
      if (!is.na(vb_err)) message("[VB] ", vb_err)
    }
    
    fit_mcmc <- tryCatch(
      brms::brm(
        y ~ ., data = df_train, family = family_cat, prior = priors,
        chains = max(2, chains), iter = max(2000, iter), warmup = floor(max(2000, iter)/2),
        refresh = 0, backend = backend,
        control = list(adapt_delta = 0.95, max_treedepth = 12, init_r = 0.1)
      ),
      error = function(e) { mcmc_err <<- conditionMessage(e); NULL }
    )
    if (!is.null(fit_mcmc)) return(list(fit = fit_mcmc, label = "BRMS_MCMC", err = NULL))
    
    list(fit = NULL, label = NA_character_, err = paste(stats::na.omit(c(vb_err, mcmc_err)), collapse = " | "))
  }
  
  for (i in seq_along(return.list$train)) {
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test),  ref = ref)
    
    df_train <- data.frame(y = y_train, X_train, check.names = FALSE) |> std_and_prune()
    df_test  <- data.frame(y = y_test,  X_test,  check.names = FALSE)
    
    df_test  <- df_test[, intersect(colnames(df_test), colnames(df_train)), drop = FALSE]
    df_train$y <- droplevels(df_train$y)
    df_test$y  <- factor(df_test$y, levels = levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- make.unique(sanitize_names(common_cols), sep = "_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols, colnames(df_test))]   <- new_names
    df_train <- df_train[, c("y", new_names), drop = FALSE]
    df_test  <- df_test[,  c("y", new_names), drop = FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    # NNET
    fit_nnet <- nnet::multinom(y ~ ., data = df_train, trace = FALSE, MaxNWts = 200000)
    probs_nnet <- stats::predict(fit_nnet, newdata = df_test, type = "probs")
    if (is.null(colnames(probs_nnet))) colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% dplyr::mutate(Model = "NNET", GridID = i, p = p, n_coefs = count_coefs_nnet(fit_nnet))
    
    # BRMS
    if (do_brms) {
      fr <- fit_brms_with_fallback(df_train, ref_lbl = ref, algo = brms_algorithm,
                                   iter = brms_iter, chains = brms_chains, backend = backend_choice)
      if (!is.null(fr$fit)) {
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop = FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% dplyr::mutate(Model = fr$label, GridID = i, p = p, n_coefs = count_coefs_brms(fr$fit))
        attr(rows[[fr$label]], "biomarkers") <- extract_biomarkers_brms(fr$fit, top_k = 100)
      } else if (!is.null(fr$err)) {
        brms_errs[[length(brms_errs) + 1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy, BalAcc, Kappa, AUC, LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()), by = "GridID")
  
  best_by_grid <- res_tbl %>% dplyr::group_by(GridID) %>% dplyr::group_modify(~ rank_one(.x) %>% dplyr::slice_head(n = 1)) %>% dplyr::ungroup()
  best_overall <- rank_one(res_tbl) %>% dplyr::slice_head(n = 1)
  
  list(results = res_tbl, best_by_grid = best_by_grid, best_overall = best_overall,
       brms_errors = unique(unlist(brms_errs)), backend = backend_choice)
}

# -------------------------------------------------------
# 6) CARGA DE INSUMOS + GRID + RETURN.LIST
# -------------------------------------------------------
dats       <- readRDS("./datos_para_modelar.rds")
tst_Data   <- readRDS("./ready_for_modeling.rds")
factors    <- dats$factores
pesos_list <- dats$pesos
mod        <- readRDS("./modelo.rds")

Z_seed <- as.matrix(factors %>% dplyr::filter(Tipo == "Real") %>% dplyr::select(where(is.numeric)))[, 1:2, drop = FALSE]
y_seed <- factors %>% dplyr::filter(Tipo == "Real") %>% dplyr::pull(Grupo) %>% factor()
names(y_seed) <- rownames(Z_seed)

v <- c(0.7, 0.98)
thr_grid <- expand.grid(q_tx = v, q_pr = v, q_me = v, q_cl = v)

return.list <- build_sets_from_grid(thr_grid, tst_Data, pesos_list, Z_seed, y_seed, mod)

detach("package:MOFA2", unload = TRUE)

# -------------------------------------------------------
# 7) EJECUCIÓN
#    Por defecto usa MCMC para evitar problemas de VB.
# -------------------------------------------------------
suppressPackageStartupMessages({ library(brms); library(nnet); library(caret); library(pROC); library(dplyr) })

EVAL <- evaluate_all(
  return.list, thr_grid,
  ref = "NP",
  do_brms = TRUE,
  brms_algorithm = "mcmc",  # seguro; usar "vb" si se desea probar VB
  brms_iter = 4000,
  brms_chains = 4
)
