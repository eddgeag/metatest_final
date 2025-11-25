# =======================================================
# Loop por cada combinación de thr_grid con salida por lista
# Métricas por modelo, número de coeficientes y ranking de "mejor" modelo
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
# Utilidad: selección de features por vista desde MOFA
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
  # threshold
  thr <- stats::quantile(aux$data$value, probs = q_thr, na.rm = TRUE)
  
  # features seleccionados
  feats <- aux$data$feature[aux$data$value >= thr]
  
  # subset train y test
  sub_train <- data_list$train[as.character(feats), , drop = FALSE]
  sub_test  <- data_list$test[as.character(feats), , drop = FALSE]
  
  # mapeo Entrez -> Symbol (opcional si metadata disponible)
  if (!is.null(metadata) &&
      all(c("EntrezGeneID", "GeneSymbol") %in% colnames(metadata))) {
    map_dic <- metadata %>% dplyr::select(EntrezGeneID, GeneSymbol) %>% distinct()
    idx <- rownames(sub_train)
    ids_map <- data.frame(EntrezGeneID = idx) %>% dplyr::left_join(map_dic, by = "EntrezGeneID")
    new_ids <- ifelse(!is.na(ids_map$GeneSymbol),
                      ids_map$GeneSymbol,
                      ids_map$EntrezGeneID)
    rownames(sub_train) <- new_ids
    rownames(sub_test)  <- new_ids
  }
  
  list(train = sub_train, test = sub_test)
}

# -------------------------------------------------------
# Proyección Z para test a partir de pesos W de MOFA
# -------------------------------------------------------
project_all_views <- function(X_list, W_list) {
  W_all <- NULL
  X_all <- NULL
  for (v in names(W_list)) {
    Xv <- X_list[[v]]
    if (is.list(Xv) && "test" %in% names(Xv))
      Xv <- Xv$test
    common_feats <- intersect(rownames(W_list[[v]]), rownames(Xv))
    if (length(common_feats) == 0)
      next
    W_sub <- as.matrix(W_list[[v]][common_feats, , drop = FALSE])
    X_sub <- as.matrix(Xv[common_feats, , drop = FALSE])
    W_all <- rbind(W_all, W_sub)
    X_all <- rbind(X_all, X_sub)
  }
  stopifnot(nrow(W_all) == nrow(X_all))
  Z_est <- solve(t(W_all) %*% W_all) %*% t(W_all) %*% X_all
  t(Z_est)
}

refactor_levels <- function(x)
  as.factor(as.character(x))

# -------------------------------------------------------
# Construcción de sets por grid de umbrales
# -------------------------------------------------------
build_sets_from_grid <- function(thr_grid,
                                 tst_Data,
                                 pesos_list,
                                 Z_seed,
                                 y_seed,
                                 mod) {
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
    q_tx <- thr_grid$q_tx[i]
    q_pr <- thr_grid$q_pr[i]
    q_me <- thr_grid$q_me[i]
    q_cl <- thr_grid$q_cl[i]
    
    tx <- select_features_view(
      mod,
      view = 1,
      q_thr = q_tx,
      data_list = tst_Data$tx,
      metadata = tst_Data$features_metadata
    )
    pr <- select_features_view(
      mod,
      view = 2,
      q_thr = q_pr,
      data_list = tst_Data$pr,
      metadata = tst_Data$features_metadata
    )
    me <- select_features_view(
      mod,
      view = 3,
      q_thr = q_me,
      data_list = tst_Data$me,
      metadata = tst_Data$features_metadata
    )
    cl <- select_features_view(
      mod,
      view = 4,
      q_thr = q_cl,
      data_list = tst_Data$cl,
      metadata = tst_Data$features_metadata
    )
    
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
# Carga de insumos
# =======================================================
dats       <- readRDS("./datos_para_modelar.rds")
tst_Data   <- readRDS("./ready_for_modeling.rds")
factors    <- dats$factores
pesos_list <- dats$pesos

Z_seed <- as.matrix(factors %>% dplyr::filter(Tipo == "Real") %>% dplyr::select(where(is.numeric)))[, 1:2, drop = FALSE]
y_seed <- factors %>% dplyr::filter(Tipo == "Real") %>% dplyr::pull(Grupo) %>% factor()
names(y_seed) <- rownames(Z_seed)
lev <- levels(y_seed)

mod <- readRDS("./modelo.rds")

v <- c(0.7, 0.98)
thr_grid <- expand.grid(
  q_tx = v,
  q_pr = v,
  q_me = v,
  q_cl = v
)

return.list <- build_sets_from_grid(thr_grid, tst_Data, pesos_list, Z_seed, y_seed, mod)

detach("package:MOFA2", unload = TRUE)

# =======================================================
# Métricas y modelos por cada elemento de return.list$train/test
# =======================================================
suppressPackageStartupMessages({
  library(brms)
  library(nnet)
  library(caret)
  library(pROC)
  library(dplyr)
})

# ---------- utilidades métricas ----------
metrics_from_probs <- function(prob_mat, y_true) {
  lev <- levels(y_true)
  stopifnot(setequal(colnames(prob_mat), lev))
  prob_mat <- prob_mat[, lev, drop = FALSE]
  pred_labels <- apply(prob_mat, 1, function(p)
    lev[which.max(p)])
  pred_class  <- factor(pred_labels, levels = lev)
  cm <- caret::confusionMatrix(pred_class, y_true)
  auc_obj <- tryCatch(
    pROC::multiclass.roc(y_true, prob_mat),
    error = function(e)
      NULL
  )
  auc_num <- if (is.null(auc_obj))
    NA_real_
  else
    as.numeric(auc_obj$auc)
  Y   <- model.matrix( ~ y_true - 1)
  colnames(Y) <- lev
  eps <- 1e-15
  logloss <- -mean(rowSums(Y * log(pmax(
    pmin(prob_mat, 1 - eps), eps
  ))))
  tibble(
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    BalAcc   = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
    Kappa    = as.numeric(cm$overall["Kappa"]),
    AUC      = auc_num,
    LogLoss  = logloss
  )
}

# ---------- utilidades brms ----------
# priors más fuertes (shrinkage tipo ridge)
make_multinom_priors <- function(y_levels, ref="NP", sd_b=0.2, sd_int=2){
  levs <- setdiff(y_levels, ref)
  do.call(c, lapply(levs, function(lev){
    dpar <- paste0("mu", lev)
    c(
      set_prior(paste0("student_t(3,0,", sd_b, ")"), class="b", dpar=dpar),
      set_prior(paste0("student_t(3,0,", sd_int,")"), class="Intercept", dpar=dpar)
    )
  }))
}
brms_probs <- function(fit, newdata) {
  pp <- posterior_epred(fit, newdata = newdata)         # draws x N x K
  M  <- apply(pp, c(2, 3), mean)
  cls <- dimnames(pp)[[3]]
  if (is.null(cls))
    cls <- fit$family$names
  colnames(M) <- cls
  M
}

count_coefs_nnet <- function(fit) {
  cf <- coef(fit) # matrix: classes x (p + 1)
  if (is.list(cf))
    cf <- do.call(rbind, cf)
  cf <- as.matrix(cf)
  # quitar intercepto si está presente en columna llamada "(Intercept)" o última
  # nnet::multinom pone intercepto como primer columna
  if (ncol(cf) > 0)
    cf_noint <- cf[, -1, drop = FALSE]
  else
    cf_noint <- cf
  sum(abs(cf_noint) > 0)
}

count_coefs_brms <- function(fit) {
  fe <- tryCatch(
    as.data.frame(brms::fixef(fit)),
    error = function(e)
      NULL
  )
  if (is.null(fe))
    return(NA_integer_)
  # filtrar betas (no interceptos)
  trms <- rownames(fe)
  keep <- !grepl("Intercept", trms)
  sum(abs(fe$Estimate[keep]) > 0)
}

# ---------- pipeline por grid ----------


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
  
  # --- backend y registro de errores BRMS ---
  detect_backend <- function() {
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
      if (!inherits(v, "try-error"))
        return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  # --- saneo de nombres para Stan ---
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
  
  # --- ranking helper ---
  rank_one <- function(df) {
    df %>% dplyr::mutate(
      r1 = rank(-BalAcc, ties.method = "min"),
      r2 = rank(LogLoss, ties.method = "min"),
      r3 = rank(-AUC, ties.method = "min"),
      RankSum = r1 + r2 + r3
    ) %>% dplyr::arrange(RankSum, dplyr::desc(BalAcc), AUC, LogLoss)
  }
  
  # --- BRMS con VB fullrank y fallback a MCMC ---
  fit_brms_with_fallback <- function(df_train,
                                     ref_lbl,
                                     algo = brms_algorithm,
                                     iter = brms_iter,
                                     chains = brms_chains,
                                     backend = backend_choice) {
    priors <- make_multinom_priors(
      levels(df_train$y),
      ref = ref_lbl,
      sd_b = 0.3,
      sd_int = 3
    )
    family_cat <- brms::categorical()
    
    # 1) VB fullrank
    if (algo %in% c("vb", "meanfield", "fullrank")) {
      vb_err <- NULL
      fit_vb <- tryCatch(
        brm(
          y ~ .,
          data = df_train,
          family = family_cat,
          prior = priors,
          algorithm = "fullrank",
          iter = max(4000, iter),
          refresh = 0,
          backend = backend,
          control = list(
            tol_rel_obj = 1e-4,
            eval_elbo = 200,
            adapt_engaged = TRUE
          )
        ),
        error = function(e) {
          vb_err <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(fit_vb)) {
        kk <- tryCatch({
          ll <- brms::loo(fit_vb)
          max(ll$diagnostics$pareto_k)
        }, error = function(e)
          Inf)
        if (is.finite(kk) &&
            kk <= 0.7)
          return(list(
            fit = fit_vb,
            label = "BRMS_VB",
            err = NULL
          ))
        vb_err <- sprintf("VB Pareto-k=%.2f", kk)
      }
      if (!is.null(vb_err))
        message("[VB] ", vb_err)
    }
    
    # 2) MCMC
    mcmc_err <- NULL
    fit_mcmc <- tryCatch(
      brm(
        y ~ .,
        data = df_train,
        family = family_cat,
        prior = priors,
        chains = max(2, chains),
        iter = max(2000, iter),
        warmup = floor(max(2000, iter) / 2),
        refresh = 0,
        backend = backend,
        control = list(
          adapt_delta = 0.95,
          max_treedepth = 12,
          init_r = 0.1
        )
      ),
      error = function(e) {
        mcmc_err <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fit_mcmc))
      return(list(
        fit = fit_mcmc,
        label = "BRMS_MCMC",
        err = NULL
      ))
    
    list(
      fit = NULL,
      label = NA_character_,
      err = paste(na.omit(c(
        vb_err, mcmc_err
      )), collapse = " | ")
    )
  }
  
  # --- loop por grids ---
  for (i in seq_along(return.list$train)) {
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test), ref = ref)
    
    df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
    df_test  <- data.frame(y = y_test, X_test, check.names = FALSE)
    
    # niveles y columnas
    df_train$y <- droplevels(df_train$y)
    df_test$y  <- factor(df_test$y, levels = levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- sanitize_names(common_cols)
    new_names <- make.unique(new_names, sep = "_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols, colnames(df_test))] <- new_names
    df_train <- df_train[, c("y", new_names), drop = FALSE]
    df_test  <- df_test[, c("y", new_names), drop = FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    # NNET
    fit_nnet <- nnet::multinom(y ~ .,
                               data = df_train,
                               trace = FALSE,
                               MaxNWts = 200000)
    probs_nnet <- stats::predict(fit_nnet, newdata = df_test, type = "probs")
    if (is.null(colnames(probs_nnet)))
      colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% dplyr::mutate(
      Model = "NNET",
      GridID = i,
      p = p,
      n_coefs = count_coefs_nnet(fit_nnet)
    )
    
    # BRMS
    if (do_brms) {
      fr <- fit_brms_with_fallback(
        df_train,
        ref_lbl = ref,
        algo = brms_algorithm,
        iter = brms_iter,
        chains = brms_chains,
        backend = backend_choice
      )
      if (!is.null(fr$fit)) {
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop = FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% dplyr::mutate(
          Model = fr$label,
          GridID = i,
          p = p,
          n_coefs = count_coefs_brms(fr$fit)
        )
      } else if (!is.null(fr$err)) {
        brms_errs[[length(brms_errs) + 1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  # --- consolidación y ranking ---
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy, BalAcc, Kappa, AUC, LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()), by = "GridID")
  
  best_by_grid <- res_tbl %>%
    dplyr::group_by(GridID) %>%
    dplyr::group_modify( ~ rank_one(.x) %>% dplyr::slice_head(n = 1)) %>%
    dplyr::ungroup()
  
  best_overall <- rank_one(res_tbl) %>% dplyr::slice_head(n = 1)
  
  list(
    results = res_tbl,
    best_by_grid = best_by_grid,
    best_overall = best_overall,
    brms_errors = unique(unlist(brms_errs)),
    backend = backend_choice
  )
}

# Ejemplo de ejecución
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
  
  # --- backend y registro de errores BRMS ---
  detect_backend <- function() {
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
      if (!inherits(v, "try-error"))
        return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  # --- saneo de nombres para Stan ---
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
  
  # --- ranking helper ---
  rank_one <- function(df) {
    df %>% dplyr::mutate(
      r1 = rank(-BalAcc, ties.method = "min"),
      r2 = rank(LogLoss, ties.method = "min"),
      r3 = rank(-AUC, ties.method = "min"),
      RankSum = r1 + r2 + r3
    ) %>% dplyr::arrange(RankSum, dplyr::desc(BalAcc), AUC, LogLoss)
  }
  
  # --- BRMS con VB fullrank y fallback a MCMC ---
  fit_brms_with_fallback <- function(df_train,
                                     ref_lbl,
                                     algo = brms_algorithm,
                                     iter = brms_iter,
                                     chains = brms_chains,
                                     backend = backend_choice) {
    priors <- make_multinom_priors(
      levels(df_train$y),
      ref = ref_lbl,
      sd_b = 0.3,
      sd_int = 3
    )
    family_cat <- brms::categorical()
    
    # 1) VB fullrank
    if (algo %in% c("vb", "meanfield", "fullrank")) {
      vb_err <- NULL
      fit_vb <- tryCatch(
        brm(
          y ~ .,
          data = df_train,
          family = family_cat,
          prior = priors,
          algorithm = "fullrank",
          iter = max(4000, iter),
          refresh = 0,
          backend = backend,
          control = list(
            tol_rel_obj = 1e-4,
            eval_elbo = 200,
            adapt_engaged = TRUE
          )
        ),
        error = function(e) {
          vb_err <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(fit_vb)) {
        kk <- tryCatch({
          ll <- brms::loo(fit_vb)
          max(ll$diagnostics$pareto_k)
        }, error = function(e)
          Inf)
        if (is.finite(kk) &&
            kk <= 0.7)
          return(list(
            fit = fit_vb,
            label = "BRMS_VB",
            err = NULL
          ))
        vb_err <- sprintf("VB Pareto-k=%.2f", kk)
      }
      if (!is.null(vb_err))
        message("[VB] ", vb_err)
    }
    
    # 2) MCMC
    mcmc_err <- NULL
    fit_mcmc <- tryCatch(
      brm(
        y ~ .,
        data = df_train,
        family = family_cat,
        prior = priors,
        chains = max(2, chains),
        iter = max(2000, iter),
        warmup = floor(max(2000, iter) / 2),
        refresh = 0,
        backend = backend,
        control = list(
          adapt_delta = 0.95,
          max_treedepth = 12,
          init_r = 0.1
        )
      ),
      error = function(e) {
        mcmc_err <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fit_mcmc))
      return(list(
        fit = fit_mcmc,
        label = "BRMS_MCMC",
        err = NULL
      ))
    
    list(
      fit = NULL,
      label = NA_character_,
      err = paste(na.omit(c(
        vb_err, mcmc_err
      )), collapse = " | ")
    )
  }
  
  # --- loop por grids ---
  for (i in seq_along(return.list$train)) {
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test), ref = ref)
    
    df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
    df_test  <- data.frame(y = y_test, X_test, check.names = FALSE)
    
    # niveles y columnas
    df_train$y <- droplevels(df_train$y)
    df_test$y  <- factor(df_test$y, levels = levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- sanitize_names(common_cols)
    new_names <- make.unique(new_names, sep = "_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols, colnames(df_test))] <- new_names
    df_train <- df_train[, c("y", new_names), drop = FALSE]
    df_test  <- df_test[, c("y", new_names), drop = FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    # NNET
    fit_nnet <- nnet::multinom(y ~ .,
                               data = df_train,
                               trace = FALSE,
                               MaxNWts = 200000)
    probs_nnet <- stats::predict(fit_nnet, newdata = df_test, type = "probs")
    if (is.null(colnames(probs_nnet)))
      colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% dplyr::mutate(
      Model = "NNET",
      GridID = i,
      p = p,
      n_coefs = count_coefs_nnet(fit_nnet)
    )
    
    # BRMS
    if (do_brms) {
      fr <- fit_brms_with_fallback(
        df_train,
        ref_lbl = ref,
        algo = brms_algorithm,
        iter = brms_iter,
        chains = brms_chains,
        backend = backend_choice
      )
      if (!is.null(fr$fit)) {
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop = FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% dplyr::mutate(
          Model = fr$label,
          GridID = i,
          p = p,
          n_coefs = count_coefs_brms(fr$fit)
        )
      } else if (!is.null(fr$err)) {
        brms_errs[[length(brms_errs) + 1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  # --- consolidación y ranking ---
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy, BalAcc, Kappa, AUC, LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()), by = "GridID")
  
  best_by_grid <- res_tbl %>%
    dplyr::group_by(GridID) %>%
    dplyr::group_modify( ~ rank_one(.x) %>% dplyr::slice_head(n = 1)) %>%
    dplyr::ungroup()
  
  best_overall <- rank_one(res_tbl) %>% dplyr::slice_head(n = 1)
  
  list(
    results = res_tbl,
    best_by_grid = best_by_grid,
    best_overall = best_overall,
    brms_errors = unique(unlist(brms_errs)),
    backend = backend_choice
  )
}
# --- antes del loop (útil en todo el script) ---
std_and_prune <- function(df){
  y <- df$y
  X <- df[, setdiff(colnames(df),"y"), drop=FALSE]
  # quitar NZV
  nzv <- caret::nearZeroVar(X)
  if(length(nzv)) X <- X[, -nzv, drop=FALSE]
  # quitar colineales (threshold conservador)
  if(ncol(X) > 2){
    cc <- stats::cor(X, use="pairwise.complete.obs")
    bad <- caret::findCorrelation(cc, cutoff=0.95, names=FALSE)
    if(length(bad)) X <- X[, -bad, drop=FALSE]
  }
  # estandarizar
  mu <- vapply(X, mean, 0.0, na.rm=TRUE)
  sd <- vapply(X,   sd, 0.0, na.rm=TRUE)
  sd[sd == 0] <- 1
  Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
  data.frame(y=y, Xs, check.names=FALSE)
}

extract_biomarkers_brms <- function(fit, top_k=50){
  fx <- as.data.frame(brms::fixef(fit, robust=TRUE))  # Estimate, Q2.5, Q97.5
  fx$term <- rownames(fx)
  fx <- fx[!grepl("Intercept", fx$term), ]
  # term viene como "b_Feature" por cada dpar mu<Clase>; mantenemos rasgos
  fx$feature <- sub("^b_", "", sub("^.*?b_", "b_", fx$term))
  fx$class   <- sub("^mu([^_]+).*", "\\1", fx$term)
  fx$nonzero <- (fx$Q2.5 > 0) | (fx$Q97.5 < 0)
  agg <- fx %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      any_nonzero = any(nonzero),
      max_abs_est = max(abs(Estimate), na.rm=TRUE),
      n_classes   = sum(nonzero),
      .groups="drop"
    ) %>%
    dplyr::arrange(dplyr::desc(any_nonzero), dplyr::desc(n_classes), dplyr::desc(max_abs_est))
  list(
    per_class = fx[fx$nonzero, c("class","feature","Estimate","Q2.5","Q97.5")],
    ranking   = head(agg, top_k)
  )
}

# Ejemplo de ejecución
EVAL <- evaluate_all <- function(return.list,
                                 thr_grid,
                                 ref = "NP",
                                 do_brms = FALSE,
                                 brms_iter = 2000,
                                 brms_chains = 2,
                                 brms_algorithm = c("mcmc", "vb")) {
  options(contrasts = c("contr.treatment", "contr.poly"))
  brms_algorithm <- match.arg(brms_algorithm)
  res_all <- list()
  
  # --- backend y registro de errores BRMS ---
  detect_backend <- function() {
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
      if (!inherits(v, "try-error"))
        return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  # --- saneo de nombres para Stan ---
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
  
  # --- ranking helper ---
  rank_one <- function(df) {
    df %>% dplyr::mutate(
      r1 = rank(-BalAcc, ties.method = "min"),
      r2 = rank(LogLoss, ties.method = "min"),
      r3 = rank(-AUC, ties.method = "min"),
      RankSum = r1 + r2 + r3
    ) %>% dplyr::arrange(RankSum, dplyr::desc(BalAcc), AUC, LogLoss)
  }
  
  # --- BRMS con VB fullrank y fallback a MCMC ---
  fit_brms_with_fallback <- function(df_train, ref_lbl, iter, chains, backend){
    priors <- make_multinom_priors(levels(df_train$y), ref=ref_lbl, sd_b=0.2, sd_int=2)
    family_cat <- brms::categorical()
    fit_mcmc <- brm(
      y ~ .,
      data = df_train,
      family = family_cat,
      prior = priors,
      chains = max(4, chains),
      iter   = max(4000, iter),
      warmup = floor(max(4000, iter)/2),
      refresh= 0,
      backend= backend,
      control = list(adapt_delta = 0.995, max_treedepth = 15)
    )
    list(fit=fit_mcmc, label="BRMS_MCMC", err=NULL)
  }
  
  # --- loop por grids ---
  for (i in seq_along(return.list$train)) {
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test), ref = ref)
    
    df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
    df_train <- std_and_prune(df_train)
    
    df_test  <- data.frame(y = y_test, X_test, check.names = FALSE)
    df_test  <- df_test[, colnames(df_train), drop=FALSE]
    
    # niveles y columnas
    df_train$y <- refactor_levels(df_train$y)
    df_test$y  <- factor(df_test$y, levels = levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- sanitize_names(common_cols)
    new_names <- make.unique(new_names, sep = "_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols, colnames(df_test))] <- new_names
    df_train <- df_train[, c("y", new_names), drop = FALSE]
    df_test  <- df_test[, c("y", new_names), drop = FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    # NNET
    fit_nnet <- nnet::multinom(y ~ .,
                               data = df_train,
                               trace = FALSE,
                               MaxNWts = 200000)
    probs_nnet <- stats::predict(fit_nnet, newdata = df_test, type = "probs")
    if (is.null(colnames(probs_nnet)))
      colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% dplyr::mutate(
      Model = "NNET",
      GridID = i,
      p = p,
      n_coefs = count_coefs_nnet(fit_nnet)
    )
    
    # BRMS
    if (do_brms) {
      fr <- fit_brms_with_fallback(
        df_train,
        ref_lbl = ref,
        algo = brms_algorithm,
        iter = brms_iter,
        chains = brms_chains,
        backend = backend_choice
      )
      if (!is.null(fr$fit)) {
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop=FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% dplyr::mutate(Model=fr$label, GridID=i, p=p, n_coefs=count_coefs_brms(fr$fit))
        # biomarcadores
        biom_i <- extract_biomarkers_brms(fr$fit, top_k=100)
        attr(rows[[fr$label]], "biomarkers") <- biom_i   # queda colgado al resultado de este grid
      } else if (!is.null(fr$err)) {
        brms_errs[[length(brms_errs) + 1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  # --- consolidación y ranking ---
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy, BalAcc, Kappa, AUC, LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()), by = "GridID")
  
  best_by_grid <- res_tbl %>%
    dplyr::group_by(GridID) %>%
    dplyr::group_modify( ~ rank_one(.x) %>% dplyr::slice_head(n = 1)) %>%
    dplyr::ungroup()
  
  best_overall <- rank_one(res_tbl) %>% dplyr::slice_head(n = 1)
  
  list(
    results = res_tbl,
    best_by_grid = best_by_grid,
    best_overall = best_overall,
    brms_errors = unique(unlist(brms_errs)),
    backend = backend_choice
  )
}

# Ejemplo de ejecución
EVAL <-evaluate_all <- function(return.list, thr_grid, ref = "NP", do_brms = FALSE,
brms_iter = 2000, brms_chains = 2, brms_algorithm = c("mcmc", "vb")) {
  options(contrasts = c("contr.treatment", "contr.poly"))
  brms_algorithm <- match.arg(brms_algorithm)
  res_all <- list()
  
  # --- backend y registro de errores BRMS ---
  detect_backend <- function() {
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
      if (!inherits(v, "try-error"))
        return("cmdstanr")
    }
    "rstan"
  }
  backend_choice <- detect_backend()
  brms_errs <- list()
  
  # --- saneo de nombres para Stan ---
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
  
  # --- ranking helper ---
  rank_one <- function(df) {
    df %>% dplyr::mutate(
      r1 = rank(-BalAcc, ties.method = "min"),
      r2 = rank(LogLoss, ties.method = "min"),
      r3 = rank(-AUC, ties.method = "min"),
      RankSum = r1 + r2 + r3
    ) %>% dplyr::arrange(RankSum, dplyr::desc(BalAcc), AUC, LogLoss)
  }
  
  # --- BRMS con VB fullrank y fallback a MCMC ---
  fit_brms_with_fallback <- function(df_train,
                                     ref_lbl,
                                     algo = brms_algorithm,
                                     iter = brms_iter,
                                     chains = brms_chains,
                                     backend = backend_choice) {
    priors <- make_multinom_priors(
      levels(df_train$y),
      ref = ref_lbl,
      sd_b = 0.3,
      sd_int = 3
    )
    family_cat <- brms::categorical()
    
    # 1) VB fullrank
    if (algo %in% c("vb", "meanfield", "fullrank")) {
      vb_err <- NULL
      fit_vb <- tryCatch(
        brm(
          y ~ .,
          data = df_train,
          family = family_cat,
          prior = priors,
          algorithm = "fullrank",
          iter = max(4000, iter),
          refresh = 0,
          backend = backend,
          control = list(
            tol_rel_obj = 1e-4,
            eval_elbo = 200,
            adapt_engaged = TRUE
          )
        ),
        error = function(e) {
          vb_err <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(fit_vb)) {
        kk <- tryCatch({
          ll <- brms::loo(fit_vb)
          max(ll$diagnostics$pareto_k)
        }, error = function(e)
          Inf)
        if (is.finite(kk) &&
            kk <= 0.7)
          return(list(
            fit = fit_vb,
            label = "BRMS_VB",
            err = NULL
          ))
        vb_err <- sprintf("VB Pareto-k=%.2f", kk)
      }
      if (!is.null(vb_err))
        message("[VB] ", vb_err)
    }
    
    # 2) MCMC
    mcmc_err <- NULL
    fit_mcmc <- tryCatch(
      brm(
        y ~ .,
        data = df_train,
        family = family_cat,
        prior = priors,
        chains = max(2, chains),
        iter = max(2000, iter),
        warmup = floor(max(2000, iter) / 2),
        refresh = 0,
        backend = backend,
        control = list(
          adapt_delta = 0.95,
          max_treedepth = 12,
          init_r = 0.1
        )
      ),
      error = function(e) {
        mcmc_err <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fit_mcmc))
      return(list(
        fit = fit_mcmc,
        label = "BRMS_MCMC",
        err = NULL
      ))
    
    list(
      fit = NULL,
      label = NA_character_,
      err = paste(na.omit(c(
        vb_err, mcmc_err
      )), collapse = " | ")
    )
  }
  
  # --- loop por grids ---
  for (i in seq_along(return.list$train)) {
    X_train <- return.list$train[[i]]
    X_test  <- return.list$test[[i]]
    y_train <- stats::relevel(refactor_levels(return.list$y_train), ref = ref)
    y_test  <- stats::relevel(refactor_levels(return.list$y_test), ref = ref)
    
    df_train <- data.frame(y = y_train, X_train, check.names = FALSE)
    df_test  <- data.frame(y = y_test, X_test, check.names = FALSE)
    
    # niveles y columnas
    df_train$y <- droplevels(df_train$y)
    df_test$y  <- factor(df_test$y, levels = levels(df_train$y))
    
    common_cols <- setdiff(intersect(colnames(df_train), colnames(df_test)), "y")
    new_names <- sanitize_names(common_cols)
    new_names <- make.unique(new_names, sep = "_")
    colnames(df_train)[match(common_cols, colnames(df_train))] <- new_names
    colnames(df_test)[match(common_cols, colnames(df_test))] <- new_names
    df_train <- df_train[, c("y", new_names), drop = FALSE]
    df_test  <- df_test[, c("y", new_names), drop = FALSE]
    
    p <- ncol(df_train) - 1L
    rows <- list()
    
    # NNET
    fit_nnet <- nnet::multinom(y ~ .,
                               data = df_train,
                               trace = FALSE,
                               MaxNWts = 200000)
    probs_nnet <- stats::predict(fit_nnet, newdata = df_test, type = "probs")
    if (is.null(colnames(probs_nnet)))
      colnames(probs_nnet) <- levels(df_test$y)
    m_nnet <- metrics_from_probs(probs_nnet, df_test$y)
    rows[["NNET"]] <- m_nnet %>% dplyr::mutate(
      Model = "NNET",
      GridID = i,
      p = p,
      n_coefs = count_coefs_nnet(fit_nnet)
    )
    
    # BRMS
    if (do_brms) {
      fr <- fit_brms_with_fallback(
        df_train,
        ref_lbl = ref,
        algo = brms_algorithm,
        iter = brms_iter,
        chains = brms_chains,
        backend = backend_choice
      )
      if (!is.null(fr$fit)) {
        P <- brms_probs(fr$fit, df_test)[, levels(df_test$y), drop = FALSE]
        m <- metrics_from_probs(P, df_test$y)
        rows[[fr$label]] <- m %>% dplyr::mutate(
          Model = fr$label,
          GridID = i,
          p = p,
          n_coefs = count_coefs_brms(fr$fit)
        )
      } else if (!is.null(fr$err)) {
        brms_errs[[length(brms_errs) + 1]] <- sprintf("Grid %d: %s", i, fr$err)
      }
    }
    
    res_all[[i]] <- dplyr::bind_rows(rows)
  }
  
  # --- consolidación y ranking ---
  res_tbl <- dplyr::bind_rows(res_all) %>%
    dplyr::mutate(dplyr::across(c(Accuracy, BalAcc, Kappa, AUC, LogLoss), as.numeric)) %>%
    dplyr::left_join(thr_grid %>% dplyr::mutate(GridID = dplyr::row_number()), by = "GridID")
  
  best_by_grid <- res_tbl %>%
    dplyr::group_by(GridID) %>%
    dplyr::group_modify( ~ rank_one(.x) %>% dplyr::slice_head(n = 1)) %>%
    dplyr::ungroup()
  
  best_overall <- rank_one(res_tbl) %>% dplyr::slice_head(n = 1)
  
  list(
    results = res_tbl,
    best_by_grid = best_by_grid,
    best_overall = best_overall,
    brms_errors = unique(unlist(brms_errs)),
    backend = backend_choice
  )
}

# Ejemplo de ejecución
EVAL <- evaluate_all(
  return.list, thr_grid,
  ref="NP",
  do_brms=TRUE,
  brms_algorithm="mcmc",   # ya no VB
  brms_iter=4000,
  brms_chains=4
)

  