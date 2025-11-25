





##==librerias====

library(mixOmics)
library(dplyr)
library(ggplot2)
library(caret)
library(pROC)
library(boot)
library(glmnet)
library(doParallel)
library(boot)
library(tibble)
library(purrr)
library(nnet)

#==funciones====
run_splsda_once <- function(seed,
                            df_train,
                            df_test,
                            ncomp = 2,
                            folds = 3,
                            nrepeat = 10,
                            keepX_grid = seq(5, 70, 5)) {
  set.seed(seed)
  
  X <- df_train[, -1]
  Y <- df_train$y
  
  # Tuneo
  tuneado <- tune.splsda(
    X,
    Y,
    ncomp = ncomp,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    test.keepX = keepX_grid,
    progressBar = FALSE
  )
  
  mdl_final <- mixOmics::splsda(
    X = X,
    Y = Y,
    ncomp = ncomp,
    keepX = tuneado$choice.keepX,
    scale = TRUE
  )
  
  # Predicciones en test
  pred <- predict(mdl_final, newdata = df_test[, -c(1)])
  y_pred <- pred$MajorityVote$max.dist[, paste0("comp", ncomp)]
  
  cm <- confusionMatrix(factor(y_pred, levels = levels(df_test$y)),
                        factor(df_test$y, levels = levels(df_test$y)))
  
  acc <- cm$overall["Accuracy"]
  bal_acc <- mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE)
  BER <- 1 - bal_acc
  
  # AUROC por componente
  aumdl <- mixOmics::auroc(
    mdl_final,
    newdata = df_test[, -c(1)],
    outcome.test = df_test$y,
    roc.comp = 2
  )
  
  # Cargas
  cargas <- plotLoadings(
    mdl_final,
    contrib = "max",
    method = "median",
    plot = FALSE,
    ndisplay = 200
  )$X
  
  # Features globales y por clase
  feats_all <- rownames(cargas)
  feats_classes <- split(rownames(cargas), cargas$GroupContrib)
  
  list(
    metrics = data.frame(
      seed = seed,
      Accuracy = acc,
      BalAccuracy = bal_acc,
      BER = BER
    ),
    feats_all = feats_all,
    feats_classes = feats_classes
  )
}


run_splsda_seeds <- function(seeds = 1:50, df_train, df_test) {
  res_list <- lapply(seeds, function(s)
    run_splsda_once(s, df_train, df_test))
  
  # Métricas
  metrics_all <- bind_rows(lapply(res_list, function(x)
    x$metrics))
  metrics_summary <- metrics_all %>% summarise(across(-seed, mean, na.rm = TRUE))
  
  # Frecuencia de features (global)
  feats_all <- unlist(lapply(res_list, function(x)
    x$feats_all))
  feats_all_tab <- as.data.frame(table(feats_all)) %>% arrange(desc(Freq))
  
  # Frecuencia de features por clase
  feats_classes_all <- lapply(res_list, function(x)
    x$feats_classes)
  feats_class_tab <- bind_rows(lapply(names(feats_classes_all[[1]]), function(cl) {
    feats <- unlist(lapply(feats_classes_all, function(fc)
      fc[[cl]]))
    as.data.frame(table(feats)) %>% mutate(Class = cl)
  })) %>% arrange(Class, desc(Freq))
  
  list(
    metrics_all = metrics_all,
    metrics_summary = metrics_summary,
    feats_all_tab = feats_all_tab,
    feats_class_tab = feats_class_tab
  )
}
boot_metrics_plsda <- function(data, indices, X_test, y_test, mdl) {
  y_true <- y_test[indices]
  Xb     <- X_test[indices, , drop = FALSE]
  
  pred <- predict(mdl, newdata = Xb)
  y_pred <- pred$MajorityVote$max.dist[, paste0("comp", 2)]
  y_pred <- factor(y_pred, levels = levels(y_test))
  
  acc <- mean(y_pred == y_true)
  
  acc_bal <- mean(sapply(levels(y_true), function(cl) {
    if (sum(y_true == cl) == 0)
      return(NA)
    mean(y_pred[y_true == cl] == cl)
  }), na.rm = TRUE)
  
  BER <- 1 - acc_bal
  c(acc, acc_bal, BER)
}


boot_metrics <- function(data,
                         indices,
                         X_test,
                         y_test,
                         fit,
                         best_lambda) {
  # indices ya es el bootstrap sample
  y_true <- y_test[indices]
  Xb     <- X_test[indices, , drop = FALSE]
  
  preds_class <- predict(fit,
                         newx = Xb,
                         s = best_lambda,
                         type = "class")
  preds_class <- factor(preds_class, levels = levels(y_test))
  
  acc <- mean(preds_class == y_true)
  
  acc_bal <- mean(sapply(levels(y_test), function(cl) {
    if (sum(y_true == cl) == 0)
      return(NA)  # evita divisiones vacías
    mean(preds_class[y_true == cl] == cl)
  }), na.rm = TRUE)
  
  BER <- 1 - acc_bal
  return(c(acc, acc_bal, BER))
}



boot_metrics_multinom <- function(data, indices, df_test, fit, feats_all) {
  y_true <- df_test$y[indices]
  Xb     <- df_test[indices, feats_all, drop = FALSE]
  
  # predicciones
  preds_class <- predict(fit, newdata = Xb)
  
  acc <- mean(preds_class == y_true)
  
  acc_bal <- mean(sapply(levels(y_true), function(cl) {
    if (sum(y_true == cl) == 0)
      return(NA)
    mean(preds_class[y_true == cl] == cl)
  }), na.rm = TRUE)
  
  BER <- 1 - acc_bal
  c(acc, acc_bal, BER)
}


# === Helper: p-valor empírico (un/bi-caudal) ===
p_emp <- function(x,
                  null,
                  side = c("greater", "less", "two.sided")) {
  side <- match.arg(side)
  if (side == "less")
    return(mean(x >= null, na.rm = TRUE))
  if (side == "greater")
    return(mean(x <= null, na.rm = TRUE))
  p <- 2 * min(mean(x >= null, na.rm = TRUE), mean(x <= null, na.rm = TRUE))
  pmin(p, 1)
}


run_full_pipeline <- function(df_train,
                              df_test,
                              freq_thr = 0.4,
                              seeds = 50) {
  # df_train <- readRDS("./df_train_final.rds")[, -c(2:3)]
  # df_test  <- readRDS("./df_test_final.rds")[, -c(2:3)]
  #
  #   freq_thr = 0.4
  #   seeds = 3
  # 1. SPLS-DA FEATURES
  results_splsda <- run_splsda_seeds(seeds = 1:seeds,
                                     df_train = df_train,
                                     df_test = df_test)
  
  results_splsda$feats_class_tab$Freqs <- results_splsda$feats_class_tab$Freq /
    max(results_splsda$feats_class_tab$Freq)
  
  sel_feats <- results_splsda$feats_class_tab %>%
    filter(Freqs >= freq_thr * max(Freqs)) %>%
    pull(feats)
  
  
  
  sel_feats <- as.character(sel_feats)
  X_train <- as.matrix(df_train[, sel_feats])
  y_train <- df_train$y
  X_test  <- as.matrix(df_test[, sel_feats])
  y_test  <- df_test$y
  
  
  # --- Bootstrap metrics for sPLS-DA ---
  X_test_splsda <- df_test[, sel_feats]
  y_test_splsda <- df_test$y
  mdl_splsda <- mixOmics::splsda(
    X = df_train[, sel_feats],
    Y = df_train$y,
    ncomp = 2,
    scale = TRUE
  )
  
  set.seed(123)
  boot_res_splsda <- boot(
    data = 1:nrow(X_test_splsda),
    statistic = function(data, indices)
      boot_metrics_plsda(data, indices, X_test_splsda, y_test_splsda, mdl_splsda),
    R = 2000
  )
  
  boot_ci_acc <- boot.ci(boot_res_splsda, index = 1, type = "perc")
  boot_ci_bal <- boot.ci(boot_res_splsda, index = 2, type = "perc")
  boot_ci_BER <- boot.ci(boot_res_splsda, index = 3, type = "perc")
  
  extract_ci <- function(ci)
    c(ci$percent[4], ci$percent[5])
  
  ci_tab_splsda <- tibble(
    Metric   = c("Accuracy", "BalAccuracy", "BER"),
    Estimate = boot_res_splsda$t0,
    CI_low   = c(
      extract_ci(boot_ci_acc)[1],
      extract_ci(boot_ci_bal)[1],
      extract_ci(boot_ci_BER)[1]
    ),
    CI_high  = c(
      extract_ci(boot_ci_acc)[2],
      extract_ci(boot_ci_bal)[2],
      extract_ci(boot_ci_BER)[2]
    )
  )
  
  # === AUCs + p-values para sPLS-DA ===
  auc_splsda <- mixOmics::auroc(
    mdl_splsda,
    newdata      = X_test_splsda,
    outcome.test = y_test_splsda,
    roc.comp     = 2,
    # Comp1 y Comp2
    plot         = FALSE
  )
  
  auc_tab_splsda <- purrr::imap_dfr(auc_splsda, function(comp_res, comp_name) {
    tibble::tibble(
      Component = comp_name,
      Class     = rownames(comp_res),
      AUC       = as.numeric(comp_res[, "AUC"]),
      p_value   = as.numeric(comp_res[, "p-value"])
    )
  })
  
  # 2. elastic (GLMNET)
  weights <- 1 / table(y_train)[y_train]
  number_folds <- min(table(df_train$y))
  n_cores <- parallel::detectCores() - 2
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = number_folds,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = TRUE,
    allowParallel = TRUE
  )
  
  set.seed(123)
  seed_ <- sample(1:10000, seeds)
  results <- list()
  
  for (s in seed_) {
    set.seed(s)
    train_df <- data.frame(y = y_train, X_train)
    cv_model <- caret::train(
      y ~ .,
      data = train_df,
      method = "glmnet",
      trControl = ctrl,
      tuneGrid = expand.grid(
        alpha = seq(0.001, 1, length = 10),
        # 0 = ridge, 1 = elastic, intermedios = elastic net
        lambda = seq(0.001, 1, length = 10)
      ),
      family = "multinomial",
      weights = weights
    )
    
    lambda_opt <- cv_model$bestTune$lambda
    coefs <- coef(cv_model$finalModel, s = lambda_opt)
    results[[as.character(s)]] <- list(
      seed   = s,
      lambda = lambda_opt,
      coefs  = coefs,
      perf   = getTrainPerf(cv_model),
      fit    = cv_model$finalModel
    )
  }
  stopCluster(cl)
  lambda_vec <- sapply(results, function(x)
    x$lambda)
  perf_df    <- bind_rows(lapply(results, function(x)
    x$perf), .id = "seed")
  best_idx   <- which.max(rowMeans(perf_df[, c("TrainAccuracy", "TrainKappa")]))
  final_entry <- results[[best_idx]]
  
  coef_freq_list <- list()
  for (s in names(results)) {
    coefs_s <- results[[s]]$coefs
    for (cls in names(coefs_s)) {
      cf <- as.matrix(coefs_s[[cls]])
      nz <- rownames(cf)[which(cf != 0)]
      nz <- nz[nz != "(Intercept)"]
      if (length(nz) > 0) {
        coef_freq_list[[length(coef_freq_list) + 1]] <-
          data.frame(seed = s,
                     class = cls,
                     feature = nz)
      }
    }
  }
  coef_freq_df <- bind_rows(coef_freq_list)
  coef_freq_tab <- coef_freq_df %>%
    group_by(class, feature) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    mutate(RelFreq = Freq / length(results))
  
  coef_freq_tab <- coef_freq_tab[coef_freq_tab$RelFreq > freq_thr, ]
  
  
  set.seed(123)
  boot_res <- boot(
    data = 1:length(y_test),
    statistic = function(data, indices)
      boot_metrics(
        data,
        indices,
        X_test,
        y_test,
        fit = final_entry$fit,
        best_lambda = final_entry$lambda
      ),
    R = 2000
  )
  
  boot_ci_acc <- boot.ci(boot_res, index = 1, type = "perc")
  boot_ci_bal <- boot.ci(boot_res, index = 2, type = "perc")
  boot_ci_BER <- boot.ci(boot_res, index = 3, type = "perc")
  
  
  
  probs_test <- predict(
    final_entry$fit,
    newx = X_test,
    s = final_entry$lambda,
    type = "response"
  )
  
  # probs_test: array muestras x clases x 1
  probs_mat <- as.matrix(probs_test[, , 1])
  
  aucs <- lapply(levels(y_test), function(cl) {
    roc_resp <- ifelse(y_test == cl, 1, 0)
    roc_obj <- roc(roc_resp, probs_mat[, cl], quiet = TRUE)
    list(class = cl,
         auc = auc(roc_obj),
         ci  = ci.auc(roc_obj))
  })
  
  
  
  
  # Bootstrap CIs + estimadores
  boot_ci_list <- list(Accuracy    = boot_ci_acc,
                       BalAccuracy = boot_ci_bal,
                       BER         = boot_ci_BER)
  
  extract_ci <- function(x)
    c(x$percent[4], x$percent[5])
  
  boot_ci_tab <- imap_dfr(seq_along(boot_ci_list), function(i, name) {
    ci_vals <- extract_ci(boot_ci_list[[i]])
    tibble(
      Metric   = names(boot_ci_list)[i],
      Estimate = boot_res$t0[i],
      CI_low   = ci_vals[1],
      CI_high  = ci_vals[2]
    )
  })
  
  # ---- AUC por clase con IC y p-valor empírico
  auc_tab <- map_dfr(levels(y_test), function(cl) {
    # Etiquetas binarias one-vs-rest
    roc_resp <- ifelse(y_test == cl, 1, 0)
    roc_obj  <- roc(roc_resp, probs_mat[, cl], quiet = TRUE)
    
    # IC DeLong
    ci_vals <- ci.auc(roc_obj)
    
    # Bootstrap para p-valor empírico vs 0.5
    boot_auc <- boot(
      data = 1:length(y_test),
      statistic = function(data, indices) {
        y_true <- y_test[indices]
        probs  <- probs_mat[indices, cl]
        roc_resp <- ifelse(y_true == cl, 1, 0)
        # chequear que haya positivos y negativos
        if (length(unique(roc_resp)) < 2)
          return(NA)
        roc_obj  <- roc(roc_resp, probs, quiet = TRUE)
        as.numeric(auc(roc_obj))
      },
      R = 2000
    )
    
    # quitar NA de bootstraps inválidos
    boot_vals <- boot_auc$t[!is.na(boot_auc$t)]
    
    pval <- if (length(boot_vals) > 0) {
      p_emp(boot_vals, 0.5, side = "greater")
    } else {
      NA
    }
    
    tibble(
      Metric   = paste0("AUC_", cl),
      Estimate = as.numeric(auc(roc_obj)),
      CI_low   = as.numeric(ci_vals[1]),
      CI_high  = as.numeric(ci_vals[3]),
      p_value  = pval
    )
  })
  
  
  # ---- Resultado final: una sola tabla
  final_tab <- bind_rows(boot_ci_tab, auc_tab)
  
  # 3. PLS-DA (mixOmics)
  feats_splsda <- results_splsda$feats_class_tab
  merged <- dplyr::inner_join(
    feats_splsda %>% dplyr::rename(Class_sPLSDA = Class, Freq_sPLSDA = Freq),
    coef_freq_tab %>% dplyr::rename(
      Class_elastic = class,
      Freq_elastic = Freq,
      RelFreq_elastic = RelFreq
    ),
    by = c("feats" = "feature")
  )
  
  feats_splsda.feats <- unique(merged$feats)
  mdl <- mixOmics::plsda(
    X = df_train[, feats_splsda.feats],
    Y = df_train$y,
    ncomp = 2,
    scale = TRUE
  )
  
  cargas1 <- plotLoadings(
    mdl,
    contrib = "max",
    method = "median",
    plot = FALSE,
    comp = 1
  )$X %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::mutate(Comp = 1)
  
  cargas2 <- plotLoadings(
    mdl,
    contrib = "max",
    method = "median",
    plot = FALSE,
    comp = 2
  )$X %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::mutate(Comp = 2)
  
  cargas <- bind_rows(cargas1, cargas2) %>%
    dplyr::select(Feature, Comp, NP, FD, SO, GroupContrib, importance) %>%
    dplyr::mutate(
      Sign_NP = sign(NP),
      Sign_FD = sign(FD),
      Sign_SO = sign(SO)
    )
  
  
  feats_all <- unique(cargas$Feature)
  
  # bootstrap metrics plsda
  X_test_plsda <- df_test[, feats_splsda.feats]
  y_test_plsda <- df_test$y
  set.seed(123)
  boot_res_plsda <- boot(
    data = 1:nrow(X_test_plsda),
    statistic = function(data, indices)
      boot_metrics_plsda(data, indices, X_test_plsda, y_test_plsda, mdl),
    R = 2000
  )
  
  boot_ci_acc <- boot.ci(boot_res_plsda, index = 1, type = "perc")
  boot_ci_bal <- boot.ci(boot_res_plsda, index = 2, type = "perc")
  boot_ci_BER <- boot.ci(boot_res_plsda, index = 3, type = "perc")
  extract_ci <- function(ci)
    c(ci$percent[4], ci$percent[5])
  ci_tab_plsda <- tibble(
    Metric   = c("Accuracy", "BalAccuracy", "BER"),
    Estimate = boot_res_plsda$t0,
    CI_low   = c(
      extract_ci(boot_ci_acc)[1],
      extract_ci(boot_ci_bal)[1],
      extract_ci(boot_ci_BER)[1]
    ),
    CI_high  = c(
      extract_ci(boot_ci_acc)[2],
      extract_ci(boot_ci_bal)[2],
      extract_ci(boot_ci_BER)[2]
    )
  )
  
  # Calcular AUCs
  auc_plsda <- mixOmics::auroc(
    mdl,
    newdata     = X_test_plsda,
    outcome.test = y_test_plsda,
    roc.comp    = 2,
    # incluye Comp1 y Comp2
    plot        = FALSE
  )
  
  # Reorganizar resultados en tabla larga
  auc_tab_plsda <- purrr::imap_dfr(auc_plsda, function(comp_res, comp_name) {
    tibble::tibble(
      Component = comp_name,
      Class     = rownames(comp_res),
      AUC       = comp_res[, "AUC"],
      p_value   = comp_res[, "p-value"]
    )
  })
  
  
  
  
  # 4. MULTINOMIAL (nnet)
  
  fit_final <- multinom(y ~ ., data = df_train[, c("y", feats_all)], trace = FALSE)
  coef_mat <- coef(fit_final)
  coef_tab <- as.data.frame(coef_mat) %>%
    tibble::rownames_to_column("Class") %>%
    tidyr::pivot_longer(-Class, names_to = "Feature", values_to = "Coef")
  ref_class <- levels(df_train$y)[1]
  coef_ref <- tibble(Class = ref_class,
                     Feature = colnames(coef_mat),
                     Coef = 0)
  coef_tab <- bind_rows(coef_tab, coef_ref) %>%
    arrange(Class, Feature) %>%
    mutate(Sign = sign(Coef))
  
  set.seed(123)
  boot_res <- boot(
    data = 1:nrow(df_test),
    statistic = function(data, indices)
      boot_metrics_multinom(data, indices, df_test, fit_final, feats_all),
    R = 2000
  )
  boot_ci_acc <- boot.ci(boot_res, index = 1, type = "perc")
  boot_ci_bal <- boot.ci(boot_res, index = 2, type = "perc")
  boot_ci_BER <- boot.ci(boot_res, index = 3, type = "perc")
  ci_tab_multinom <- tibble(
    Metric   = c("Accuracy", "BalAccuracy", "BER"),
    Estimate = boot_res$t0,
    CI_low   = c(
      extract_ci(boot_ci_acc)[1],
      extract_ci(boot_ci_bal)[1],
      extract_ci(boot_ci_BER)[1]
    ),
    CI_high  = c(
      extract_ci(boot_ci_acc)[2],
      extract_ci(boot_ci_bal)[2],
      extract_ci(boot_ci_BER)[2]
    )
  )
  
  
  # === Probabilidades multinomial ===
  probs_multinom <- predict(fit_final, newdata = df_test[, feats_all, drop = FALSE], type = "probs")
  
  # ---- AUC por clase con IC y p-valor empírico
  auc_tab_multinom <- map_dfr(levels(df_test$y), function(cl) {
    # Respuesta binaria one-vs-rest
    roc_resp <- ifelse(df_test$y == cl, 1, 0)
    roc_obj  <- roc(roc_resp, probs_multinom[, cl], quiet = TRUE)
    
    # IC DeLong
    ci_vals <- ci.auc(roc_obj)
    
    # Bootstrap para p-valor empírico vs 0.5
    boot_auc <- boot(
      data = 1:nrow(df_test),
      statistic = function(data, indices) {
        y_true <- df_test$y[indices]
        probs  <- probs_multinom[indices, cl]
        roc_resp <- ifelse(y_true == cl, 1, 0)
        if (length(unique(roc_resp)) < 2)
          return(NA_real_)
        roc_obj <- roc(roc_resp, probs, quiet = TRUE)
        as.numeric(auc(roc_obj))
      },
      R = 2000
    )
    
    boot_vals <- boot_auc$t[!is.na(boot_auc$t)]
    
    pval <- if (length(boot_vals) > 0) {
      p_emp(boot_vals, 0.5, side = "greater")  # H1: AUC > 0.5
    } else {
      NA_real_
    }
    
    tibble(
      Metric   = paste0("AUC_", cl),
      Estimate = as.numeric(auc(roc_obj)),
      CI_low   = as.numeric(ci_vals[1]),
      CI_high  = as.numeric(ci_vals[3]),
      p_value  = pval
    )
  })
  
  # Consolidar con las métricas bootstrap de Accuracy, BalAccuracy, BER
  final_tab_multinom <- bind_rows(ci_tab_multinom, auc_tab_multinom)
  
  final_tab_multinom
  
  # 5. CONSOLIDADO FINAL
  K <- nlevels(y_test_splsda)
  p_acc_splsda <- p_emp(boot_res_splsda$t[, 1], 1 / K, "greater")
  p_bal_splsda <- p_emp(boot_res_splsda$t[, 2], 1 / K, "greater")
  p_ber_splsda <- p_emp(boot_res_splsda$t[, 3], 1 - 1 / K, "less")
  
  
  ci_tab_splsda <- ci_tab_splsda %>%
    mutate(p_value = c(p_acc_splsda, p_bal_splsda, p_ber_splsda))
  
  boot_res_elastic <- boot_res           # renombra el objeto existente
  K <- nlevels(y_test)
  p_acc_en  <- p_emp(boot_res_elastic$t[, 1], 1 / K, "greater")
  p_bal_en  <- p_emp(boot_res_elastic$t[, 2], 1 / K, "greater")
  p_ber_en  <- p_emp(boot_res_elastic$t[, 3], 1 - 1 / K, "less")
  
  
  K <- nlevels(y_test_plsda)
  p_acc_pls  <- p_emp(boot_res_plsda$t[, 1], 1 / K, "greater")
  p_bal_pls  <- p_emp(boot_res_plsda$t[, 2], 1 / K, "greater")
  p_ber_pls  <- p_emp(boot_res_plsda$t[, 3], 1 - 1 / K, "less")
  
  ci_tab_plsda <- ci_tab_plsda %>%
    mutate(p_value = c(p_acc_pls, p_bal_pls, p_ber_pls))
  
  
  
  boot_ci_tab <- boot_ci_tab %>%
    mutate(p_value = c(p_acc_en, p_bal_en, p_ber_en))
  
  boot_res_multinom <- boot_res         # renombra el objeto existente
  K <- nlevels(df_test$y)
  p_acc_mn  <- p_emp(boot_res_multinom$t[, 1], 1 / K, "greater")
  p_bal_mn  <- p_emp(boot_res_multinom$t[, 2], 1 / K, "greater")
  p_ber_mn  <- p_emp(boot_res_multinom$t[, 3], 1 - 1 / K, "less")
  
  ci_tab_multinom <- ci_tab_multinom %>%
    mutate(p_value = c(p_acc_mn, p_bal_mn, p_ber_mn))
  
  
  # --- Consolidación de métricas ---
  
  # Añadir columna Modelo a cada tabla
  tab_splsda    <- ci_tab_splsda    %>% mutate(Model = "sPLSDA")
  tab_elastic   <- boot_ci_tab      %>% mutate(Model = "ElasticNet")
  tab_plsda     <- ci_tab_plsda     %>% mutate(Model = "PLSDA")
  tab_multinom  <- ci_tab_multinom  %>% mutate(Model = "Multinom")
  
  # Añadir las AUCs (ya incluyen p_value)
  auc_splsda_tidy   <- auc_tab_splsda   %>% mutate(Model = "sPLSDA")
  auc_elastic_tidy  <- auc_tab          %>% mutate(Model = "ElasticNet")
  auc_multinom_tidy <- auc_tab_multinom %>% mutate(Model = "Multinom")
  
  # Nota: auc_plsda ya está separado como auc_tab_plsda
  auc_plsda_tidy    <- auc_tab_plsda    %>% mutate(Model = "PLSDA")
  
  # Unificar métricas de accuracy + aucs
  metrics_all_long <- bind_rows(
    tab_splsda,
    tab_elastic,
    tab_plsda,
    tab_multinom,
    auc_splsda_tidy,
    auc_elastic_tidy,
    auc_multinom_tidy,
    auc_plsda_tidy
  )
  
  # Reordenar columnas
  metrics_all_long <- metrics_all_long %>%
    select(Model,
           Metric,
           Estimate,
           CI_low,
           CI_high,
           p_value,
           everything())
  
  # Vista final
  print(metrics_all_long)
  
  freq_splsda <- results_splsda$feats_class_tab
  freq_elastic  <- coef_freq_tab
  plsda_loadings <- cargas %>%
    dplyr::select(Feature, NP, FD, SO, GroupContrib, importance, Comp) %>%
    dplyr::mutate(
      Sign_NP = sign(NP),
      Sign_FD = sign(FD),
      Sign_SO = sign(SO)
    )
  
  coef_final_tab <- coef_tab
  
  # --- consolidación de todas las fuentes
  
  # 1. Selecciona columnas clave de cada método
  splsda_tab <- results_splsda$feats_class_tab %>%
    group_by(feats) %>%
    summarise(
      Freq_sPLSDA  = max(Freq, na.rm = TRUE),
      # no sumar
      Class_sPLSDA = paste0(unique(Class), collapse = ";"),
      # mantener trazabilidad
      Freqs        = max(Freqs, na.rm = TRUE),
      # conserva escala [0,1]
      .groups = "drop"
    ) %>%
    rename(Feature = feats)
  
  elastic_tab <- coef_freq_tab %>%
    group_by(feature) %>%
    summarise(
      Class_elastic   = paste0(unique(class), collapse = ";"),
      Freq_elastic    = max(Freq, na.rm = TRUE),
      # no sumar
      RelFreq_elastic = max(RelFreq, na.rm = TRUE),
      # se mantiene en [0,1]
      .groups = "drop"
    ) %>%
    rename(Feature = feature)
  
  plsda_tab <- plsda_loadings %>%
    dplyr::select(Feature,
                  Comp,
                  GroupContrib,
                  NP,
                  FD,
                  SO,
                  Sign_NP,
                  Sign_FD,
                  Sign_SO,
                  importance)
  
  # plsda_tab <- plsda_tab[abs(plsda_tab$importance) > quantile(abs(plsda_tab$importance),0.6),]
  
  multinom_tab <- coef_tab %>%
    dplyr::select(Feature, Class_Multinom = Class, Coef, Sign)
  
  # 2. Resumir multinomial → 1 fila por feature, concatenando info
  multinom_summary <- multinom_tab %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      Multinom_Class = paste0(Class_Multinom, collapse = ";"),
      Multinom_Coef  = paste0(round(Coef, 3), collapse = ";"),
      Multinom_Sign  = paste0(Sign, collapse = ";"),
      .groups = "drop"
    )
  
  # 3. Resumir PLSDA → 1 fila por feature, agrupando comp
  # 3. Resumir PLSDA → 1 fila por feature, agrupando comp
  plsda_summary <- plsda_tab %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      PLSDA_Comp     = paste0(Comp, collapse = ";"),
      PLSDA_Class    = paste0(GroupContrib, collapse = ";"),
      PLSDA_Sign_NP  = paste0(Sign_NP, collapse = ";"),
      PLSDA_Sign_FD  = paste0(Sign_FD, collapse = ";"),
      PLSDA_Sign_SO  = paste0(Sign_SO, collapse = ";"),
      PLSDA_importance_mean = mean(abs(importance), na.rm = TRUE),
      PLSDA_importance_sum  = sum(abs(importance), na.rm = TRUE),
      PLSDA_importance_max  = max(abs(importance), na.rm = TRUE),
      # signo dominante según el comp con mayor |importance|
      PLSDA_sign_dom = sign(importance[which.max(abs(importance))]),
      .groups = "drop"
    )
  
  
  # 4. Unir todo
  # 4. Unir todo
  features_final <- splsda_tab %>%
    dplyr::left_join(elastic_tab, by = "Feature") %>%
    dplyr::left_join(plsda_summary, by = "Feature") %>%
    dplyr::left_join(multinom_summary, by = "Feature")
  
  # --- Definir criterios separados por bloque ---
  features_final <- features_final %>%
    mutate(
      # Criterios basados en frecuencias
      crit_splsda   = ifelse(!is.na(Freqs) &
                               Freqs >= freq_thr, 1, 0),
      crit_elastic  = ifelse(!is.na(RelFreq_elastic) &
                               RelFreq_elastic >= freq_thr, 1, 0),
      
      # Criterio basado en importancias (usa la suma, pero puedes cambiar a mean o max)
      crit_plsda    = ifelse(
        !is.na(PLSDA_importance_sum) &
          PLSDA_importance_sum >= quantile(
            features_final$PLSDA_importance_sum,
            probs = freq_thr,
            na.rm = TRUE
          ),
        1,
        0
      ),
      
      # Criterio basado en coeficientes multinomiales (≠ 0)
      crit_multinom = ifelse(!is.na(Multinom_Coef) &
                               max(abs(
                                 as.numeric(strsplit(Multinom_Coef, ";")[[1]])
                               ), na.rm = TRUE) >= 0.2, 1, 0)      ,
      # Score global
      Score     = crit_splsda + crit_elastic + crit_plsda + crit_multinom,
      Candidate = ifelse(Score >= 3, "Yes", "No")
    )
  
  
  return(list(metrics = metrics_all_long, features = features_final))
}

#=====Ejecución=====


df_train <- readRDS("./df_train_final.rds")
df_test  <- readRDS("./df_test_final.rds")



t0 <- Sys.time()

sin_factores_7 <- run_full_pipeline(
  df_train = df_train[, -c(2:3)],
  df_test = df_test[, -c(2:3)],
  freq_thr = 0.7,
  seeds = 100
)

sin_factores_65 <- run_full_pipeline(
  df_train = df_train[, -c(2:3)],
  df_test = df_test[, -c(2:3)],
  freq_thr = 0.65,
  seeds = 100
)

sin_factores_6 <- run_full_pipeline(
  df_train = df_train[, -c(2:3)],
  df_test = df_test[, -c(2:3)],
  freq_thr = 0.6,
  seeds = 100
)
t1 <- Sys.time()
print(t1 - t0)
# #
#
# t2 <- Sys.time()
# 
# con_factores <- run_full_pipeline(
#   df_train = df_train,
#   df_test = df_test,
#   freq_thr = 0.7,
#   seeds = 50
# )
# 
# 
# t3 <- Sys.time()
# 
# print(t3 - t2)


# res <- list(sin_factores = sin_factores, con_factores = con_factores)


# metrics <- sin_factores$metrics

# metrics


saveRDS(sin_factores_6, "./biomarcadores_6.rds")
saveRDS(sin_factores_65, "./biomarcadores_65.rds")
saveRDS(sin_factores_7, "./biomarcadores_7.rds")
