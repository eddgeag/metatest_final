# ============================================================
# VALIDACIÓN INTERNA DEL PANEL FINAL
# Reconstrucción de variables originales sin integrar
# + z-transformation
# + modelo multivariado bayesiano por condición
# + contraparte univariante rápida con glm()
#
# No guarda archivos
# ============================================================

library(data.table)
library(brms)
library(posterior)

# ============================================================
# PARTE I. RECONSTRUIR VARIABLES ORIGINALES SIN INTEGRAR
# ============================================================

model <- readRDS("./ready_for_modeling.rds")

# ------------------------------------------------------------
# 1. Función para limpiar nombres
# ------------------------------------------------------------

clean_names <- function(x) {
  x <- gsub("-", "_", x)
  x <- gsub("/", "_", x)
  x <- gsub("\\+", "_", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("%", "X", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ------------------------------------------------------------
# 2. Metadata transcriptómica
# ------------------------------------------------------------

tx_meta <- as.data.table(model$features_metadata)

tx_meta_small <- unique(
  tx_meta[, .(
    EntrezGeneID,
    GeneSymbol,
    GeneName,
    RefSeqAccession,
    EnsemblID,
    GenomicCoordinates,
    Cytoband,
    Description
  )]
)

tx_meta_small[
  is.na(GeneSymbol) | GeneSymbol == "",
  GeneSymbol := EntrezGeneID
]

tx_meta_small <- tx_meta_small[!duplicated(EntrezGeneID)]

# ------------------------------------------------------------
# 3. Función para reconstruir valores originales por bloque
# ------------------------------------------------------------

unscale_block_to_wide <- function(block, block_name, tx_meta_small = NULL) {
  
  train_z <- block$train
  test_z  <- block$test
  
  center <- attr(train_z, "scaled:center")
  scale  <- attr(train_z, "scaled:scale")
  
  if (is.null(center) || is.null(scale)) {
    stop(paste0("El bloque ", block_name, " no tiene scaled:center o scaled:scale."))
  }
  
  center <- center[rownames(train_z)]
  scale  <- scale[rownames(train_z)]
  
  train_original <- sweep(train_z, 1, scale, "*")
  train_original <- sweep(train_original, 1, center, "+")
  
  test_original <- sweep(test_z, 1, scale, "*")
  test_original <- sweep(test_original, 1, center, "+")
  
  all_original <- cbind(train_original, test_original)
  
  wide <- as.data.table(t(all_original), keep.rownames = "id")
  
  old_names <- setdiff(names(wide), "id")
  
  if (block_name == "tx") {
    
    annot <- data.table(
      EntrezGeneID = old_names
    )
    
    annot <- merge(
      annot,
      tx_meta_small[, .(EntrezGeneID, GeneSymbol)],
      by = "EntrezGeneID",
      all.x = TRUE,
      sort = FALSE
    )
    
    annot[
      is.na(GeneSymbol) | GeneSymbol == "",
      GeneSymbol := EntrezGeneID
    ]
    
    new_names <- clean_names(annot$GeneSymbol)
    
  } else {
    
    new_names <- clean_names(old_names)
  }
  
  new_names <- make.unique(new_names, sep = "_")
  
  setnames(wide, old_names, new_names)
  
  wide[]
}

# ------------------------------------------------------------
# 4. Reconstruir bloques en WIDE sin prefijos
# ------------------------------------------------------------

tx_wide <- unscale_block_to_wide(
  block = model$tx,
  block_name = "tx",
  tx_meta_small = tx_meta_small
)

pr_wide <- unscale_block_to_wide(
  block = model$pr,
  block_name = "pr"
)

me_wide <- unscale_block_to_wide(
  block = model$me,
  block_name = "me"
)

cl_wide <- unscale_block_to_wide(
  block = model$cl,
  block_name = "cl"
)

# ------------------------------------------------------------
# 5. Revisar duplicados entre bloques
# ------------------------------------------------------------

all_var_names <- c(
  setdiff(names(tx_wide), "id"),
  setdiff(names(pr_wide), "id"),
  setdiff(names(me_wide), "id"),
  setdiff(names(cl_wide), "id")
)

duplicados_globales <- all_var_names[duplicated(all_var_names)]

if (length(duplicados_globales) > 0) {
  warning(
    "Hay nombres repetidos entre bloques. Se harán únicos al unir: ",
    paste(unique(duplicados_globales), collapse = ", ")
  )
}

# ------------------------------------------------------------
# 6. Unir todos los bloques por muestra
# ------------------------------------------------------------

vars_sin_integrar_wide <- Reduce(
  function(x, y) merge(x, y, by = "id", all = TRUE),
  list(tx_wide, pr_wide, me_wide, cl_wide)
)

cols <- names(vars_sin_integrar_wide)
cols[-1] <- make.unique(cols[-1], sep = "_")
setnames(vars_sin_integrar_wide, names(vars_sin_integrar_wide), cols)

# ------------------------------------------------------------
# 7. Añadir variable dependiente
# ------------------------------------------------------------

grupos <- rbindlist(
  list(
    as.data.table(model$grupo_train),
    as.data.table(model$grupo_test)
  ),
  fill = TRUE
)

grupos[, grupo := factor(grupo, levels = c("NP", "FD", "SO"))]

vars_sin_integrar_wide <- merge(
  grupos,
  vars_sin_integrar_wide,
  by = "id",
  all.x = TRUE
)

setcolorder(vars_sin_integrar_wide, c("id", "grupo"))
setorder(vars_sin_integrar_wide, grupo, id)

# ------------------------------------------------------------
# 8. Tabla de equivalencia para transcriptómica
# ------------------------------------------------------------

tx_annotation_used <- data.table(
  original_variable = model$tx$vars
)

tx_annotation_used <- merge(
  tx_annotation_used,
  tx_meta_small,
  by.x = "original_variable",
  by.y = "EntrezGeneID",
  all.x = TRUE,
  sort = FALSE
)

tx_annotation_used[
  is.na(GeneSymbol) | GeneSymbol == "",
  GeneSymbol := original_variable
]

tx_annotation_used[, final_column := clean_names(GeneSymbol)]
tx_annotation_used[, final_column := make.unique(final_column, sep = "_")]

# ------------------------------------------------------------
# 9. Revisión rápida
# ------------------------------------------------------------

dim(vars_sin_integrar_wide)
table(vars_sin_integrar_wide$grupo)
names(vars_sin_integrar_wide)


# ============================================================
# PARTE II. PANEL FINAL Y Z-TRANSFORMATION
# ============================================================

dt <- as.data.table(vars_sin_integrar_wide)
dt[, grupo := factor(grupo, levels = c("NP", "FD", "SO"))]

vars_panel <- c(
  "IL_2",
  "IL_8",
  "VEGF_A",
  "EIF1AD",
  "FAM98A",
  "neg_271",
  "neg_38",
  "Osteonectin",
  "pos_220",
  "ZNF239",
  "Dias",
  "kg_Ginoide",
  "Masa_grasa_0M",
  "neg_754",
  "Peso",
  "pos_935"
)

vars_presentes <- intersect(vars_panel, names(dt))
vars_faltantes <- setdiff(vars_panel, names(dt))

if (length(vars_faltantes) > 0) {
  warning(
    "Variables no encontradas en vars_sin_integrar_wide: ",
    paste(vars_faltantes, collapse = ", ")
  )
}

dt_panel <- dt[, c("id", "grupo", vars_presentes), with = FALSE]

vars_numericas <- vars_presentes[
  vapply(dt_panel[, ..vars_presentes], is.numeric, logical(1))
]

vars_no_numericas <- setdiff(vars_presentes, vars_numericas)

if (length(vars_no_numericas) > 0) {
  warning(
    "Variables no numéricas eliminadas: ",
    paste(vars_no_numericas, collapse = ", ")
  )
}

vars_constantes <- vars_numericas[
  vapply(dt_panel[, ..vars_numericas], function(x) {
    length(unique(na.omit(x))) <= 1
  }, logical(1))
]

if (length(vars_constantes) > 0) {
  warning(
    "Variables constantes eliminadas: ",
    paste(vars_constantes, collapse = ", ")
  )
}

vars_usar <- setdiff(vars_numericas, vars_constantes)

if (length(vars_usar) < 2) {
  stop("Hay menos de 2 variables válidas para el análisis.")
}

dt_z_final <- copy(dt_panel[, c("id", "grupo", vars_usar), with = FALSE])

for (v in vars_usar) {
  dt_z_final[, (v) := as.numeric(scale(get(v)))]
}

z_check <- dt_z_final[
  ,
  lapply(.SD, function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  }),
  .SDcols = vars_usar
]

z_check


# ============================================================
# PARTE III. FUNCIONES AUXILIARES
# ============================================================

calc_auc <- function(y, p) {
  
  y <- as.numeric(y)
  p <- as.numeric(p)
  
  ok <- complete.cases(y, p)
  y <- y[ok]
  p <- p[ok]
  
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  
  r <- rank(p)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  
  auc
}

safe_logloss <- function(y, p) {
  
  y <- as.numeric(y)
  p <- as.numeric(p)
  
  ok <- complete.cases(y, p)
  y <- y[ok]
  p <- p[ok]
  
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

summarise_beta_draws <- function(draws, delta = 0.20) {
  
  draws <- as.numeric(draws)
  
  beta_median <- median(draws, na.rm = TRUE)
  abs_beta_median <- median(abs(draws), na.rm = TRUE)
  
  prob_pos <- mean(draws > 0, na.rm = TRUE)
  prob_neg <- mean(draws < 0, na.rm = TRUE)
  
  prob_abs_gt_delta <- mean(abs(draws) > delta, na.rm = TRUE)
  
  score <- abs_beta_median * prob_abs_gt_delta
  
  data.table(
    beta_median = beta_median,
    abs_beta_median = abs_beta_median,
    prob_pos = prob_pos,
    prob_neg = prob_neg,
    prob_abs_gt_delta = prob_abs_gt_delta,
    score = score
  )
}

summarise_glm_beta <- function(beta, se, delta = 0.20) {
  
  if (is.na(beta) || is.na(se) || se <= 0) {
    return(data.table(
      beta_median = NA_real_,
      abs_beta_median = NA_real_,
      prob_pos = NA_real_,
      prob_neg = NA_real_,
      prob_abs_gt_delta = NA_real_,
      score = NA_real_,
      p_value = NA_real_
    ))
  }
  
  z_value <- beta / se
  p_value <- 2 * pnorm(-abs(z_value))
  
  # Aproximación normal al posterior/frecuentista:
  # beta_hat ~ Normal(beta, se)
  prob_pos <- pnorm(beta / se)
  prob_neg <- 1 - prob_pos
  
  prob_abs_gt_delta <- pnorm((-delta - beta) / se) +
    (1 - pnorm((delta - beta) / se))
  
  score <- abs(beta) * prob_abs_gt_delta
  
  data.table(
    beta_median = beta,
    abs_beta_median = abs(beta),
    prob_pos = prob_pos,
    prob_neg = prob_neg,
    prob_abs_gt_delta = prob_abs_gt_delta,
    score = score,
    p_value = p_value
  )
}


# ============================================================
# PARTE IV. MODELO MULTIVARIADO BAYESIANO POR CONDICIÓN
# ============================================================

run_multivariate_condition_bayes <- function(data,
                                             condition,
                                             vars,
                                             iter = 1500,
                                             warmup = 750,
                                             chains = 2,
                                             cores = 2,
                                             seed = 123,
                                             use_cmdstanr = FALSE) {
  
  d <- copy(data)
  d[, y := ifelse(grupo == condition, 1, 0)]
  
  d_model <- d[, c("y", vars), with = FALSE]
  d_model <- na.omit(d_model)
  
  if (length(unique(d_model$y)) < 2) {
    stop("La condición ", condition, " no tiene ambas clases.")
  }
  
  formula_txt <- paste("y ~", paste(vars, collapse = " + "))
  f <- as.formula(formula_txt)
  
  backend_arg <- if (use_cmdstanr) "cmdstanr" else NULL
  
  fit <- brm(
    formula = f,
    data = d_model,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(0, 1), class = "Intercept"),
      prior(normal(0, 0.30), class = "b")
    ),
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    backend = backend_arg,
    refresh = 0,
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 11
    )
  )
  
  draws <- as_draws_df(fit)
  
  coef_results <- rbindlist(lapply(vars, function(v) {
    
    par_name <- paste0("b_", v)
    
    if (!par_name %in% names(draws)) {
      return(data.table(
        condicion = condition,
        variable = v,
        beta_median_multi = NA_real_,
        abs_beta_median_multi = NA_real_,
        prob_pos_multi = NA_real_,
        prob_neg_multi = NA_real_,
        prob_abs_gt_delta_multi = NA_real_,
        score_multi = NA_real_
      ))
    }
    
    s <- summarise_beta_draws(draws[[par_name]], delta = 0.20)
    
    data.table(
      condicion = condition,
      variable = v,
      beta_median_multi = s$beta_median,
      abs_beta_median_multi = s$abs_beta_median,
      prob_pos_multi = s$prob_pos,
      prob_neg_multi = s$prob_neg,
      prob_abs_gt_delta_multi = s$prob_abs_gt_delta,
      score_multi = s$score
    )
  }))
  
  epred <- posterior_epred(fit, newdata = d_model)
  pred_prob <- colMeans(epred)
  
  auc_aparente <- calc_auc(d_model$y, pred_prob)
  
  pred_class <- ifelse(pred_prob >= 0.5, 1, 0)
  
  sensibilidad <- ifelse(
    sum(d_model$y == 1) > 0,
    sum(pred_class == 1 & d_model$y == 1) / sum(d_model$y == 1),
    NA_real_
  )
  
  especificidad <- ifelse(
    sum(d_model$y == 0) > 0,
    sum(pred_class == 0 & d_model$y == 0) / sum(d_model$y == 0),
    NA_real_
  )
  
  balanced_accuracy <- mean(c(sensibilidad, especificidad), na.rm = TRUE)
  logloss <- safe_logloss(d_model$y, pred_prob)
  
  metrics <- data.table(
    condicion = condition,
    n = nrow(d_model),
    n_caso = sum(d_model$y == 1),
    n_resto = sum(d_model$y == 0),
    auc_aparente = auc_aparente,
    sensibilidad_aparente = sensibilidad,
    especificidad_aparente = especificidad,
    balanced_accuracy_aparente = balanced_accuracy,
    logloss_aparente = logloss
  )
  
  list(
    fit = fit,
    coef_results = coef_results,
    metrics = metrics
  )
}


# ============================================================
# PARTE V. CONTRAPARTE UNIVARIANTE RÁPIDA CON GLM
# ============================================================

run_univariate_condition_glm <- function(data,
                                         condition,
                                         variable,
                                         delta = 0.20) {
  
  d <- copy(data)
  d[, y := ifelse(grupo == condition, 1, 0)]
  
  d_model <- d[, .(
    y = y,
    x = get(variable)
  )]
  
  d_model <- na.omit(d_model)
  
  if (length(unique(d_model$y)) < 2) {
    return(data.table(
      condicion = condition,
      variable = variable,
      n = nrow(d_model),
      n_caso = sum(d_model$y == 1),
      n_resto = sum(d_model$y == 0),
      beta_median_uni = NA_real_,
      abs_beta_median_uni = NA_real_,
      prob_pos_uni = NA_real_,
      prob_neg_uni = NA_real_,
      prob_abs_gt_delta_uni = NA_real_,
      score_uni = NA_real_,
      auc_uni = NA_real_,
      logloss_uni = NA_real_,
      p_value_uni = NA_real_,
      comentario = "una sola clase"
    ))
  }
  
  fit <- tryCatch(
    glm(
      y ~ x,
      data = d_model,
      family = binomial()
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.table(
      condicion = condition,
      variable = variable,
      n = nrow(d_model),
      n_caso = sum(d_model$y == 1),
      n_resto = sum(d_model$y == 0),
      beta_median_uni = NA_real_,
      abs_beta_median_uni = NA_real_,
      prob_pos_uni = NA_real_,
      prob_neg_uni = NA_real_,
      prob_abs_gt_delta_uni = NA_real_,
      score_uni = NA_real_,
      auc_uni = NA_real_,
      logloss_uni = NA_real_,
      p_value_uni = NA_real_,
      comentario = "glm falló"
    ))
  }
  
  sm <- summary(fit)$coefficients
  
  if (!"x" %in% rownames(sm)) {
    return(data.table(
      condicion = condition,
      variable = variable,
      n = nrow(d_model),
      n_caso = sum(d_model$y == 1),
      n_resto = sum(d_model$y == 0),
      beta_median_uni = NA_real_,
      abs_beta_median_uni = NA_real_,
      prob_pos_uni = NA_real_,
      prob_neg_uni = NA_real_,
      prob_abs_gt_delta_uni = NA_real_,
      score_uni = NA_real_,
      auc_uni = NA_real_,
      logloss_uni = NA_real_,
      p_value_uni = NA_real_,
      comentario = "coeficiente no disponible"
    ))
  }
  
  beta <- sm["x", "Estimate"]
  se <- sm["x", "Std. Error"]
  
  s <- summarise_glm_beta(beta, se, delta = delta)
  
  pred_prob <- as.numeric(predict(fit, type = "response"))
  
  auc_uni <- calc_auc(d_model$y, pred_prob)
  logloss_uni <- safe_logloss(d_model$y, pred_prob)
  
  data.table(
    condicion = condition,
    variable = variable,
    n = nrow(d_model),
    n_caso = sum(d_model$y == 1),
    n_resto = sum(d_model$y == 0),
    beta_median_uni = s$beta_median,
    abs_beta_median_uni = s$abs_beta_median,
    prob_pos_uni = s$prob_pos,
    prob_neg_uni = s$prob_neg,
    prob_abs_gt_delta_uni = s$prob_abs_gt_delta,
    score_uni = s$score,
    auc_uni = auc_uni,
    logloss_uni = logloss_uni,
    p_value_uni = s$p_value,
    comentario = "OK"
  )
}


# ============================================================
# PARTE VI. EJECUTAR MODELOS
# ============================================================

condiciones <- c("NP", "FD", "SO")

fits_multi <- list()
resultados_multi <- list()
metricas_panel <- list()

# Cambia a TRUE solo si ya tienes cmdstanr configurado.
use_cmdstanr <- FALSE

for (cc in condiciones) {
  
  message("============================================")
  message("Modelo bayesiano multivariado para condición: ", cc)
  message("============================================")
  
  res <- tryCatch(
    run_multivariate_condition_bayes(
      data = dt_z_final,
      condition = cc,
      vars = vars_usar,
      iter = 1500,
      warmup = 750,
      chains = 2,
      cores = 2,
      seed = 123,
      use_cmdstanr = use_cmdstanr
    ),
    error = function(e) {
      message("ERROR multivariado en ", cc, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(res)) {
    fits_multi[[cc]] <- res$fit
    resultados_multi[[cc]] <- res$coef_results
    metricas_panel[[cc]] <- res$metrics
  }
}

resultados_multi_dt <- rbindlist(resultados_multi, fill = TRUE)
metricas_panel_dt <- rbindlist(metricas_panel, fill = TRUE)

# ------------------------------------------------------------
# Univariantes rápidos
# ------------------------------------------------------------

resultados_uni_dt <- rbindlist(
  lapply(condiciones, function(cc) {
    
    rbindlist(lapply(vars_usar, function(v) {
      
      message("GLM univariante: ", cc, " ~ ", v)
      
      run_univariate_condition_glm(
        data = dt_z_final,
        condition = cc,
        variable = v,
        delta = 0.20
      )
      
    }), fill = TRUE)
    
  }),
  fill = TRUE
)


# ============================================================
# PARTE VII. INTEGRAR EVIDENCIA MULTIVARIADA + UNIVARIANTE
# ============================================================

validacion_biomarcadores <- merge(
  resultados_uni_dt,
  resultados_multi_dt,
  by = c("condicion", "variable"),
  all = TRUE
)

validacion_biomarcadores[, signo_consistente := fifelse(
  is.na(beta_median_uni) | is.na(beta_median_multi),
  NA_real_,
  fifelse(
    sign(beta_median_uni) == sign(beta_median_multi),
    1,
    0.5
  )
)]

validacion_biomarcadores[, score_consistencia := sqrt(score_uni * score_multi) * signo_consistente]

validacion_biomarcadores[
  !is.na(score_consistencia),
  score_consistencia_norm := score_consistencia / max(score_consistencia, na.rm = TRUE),
  by = condicion
]

setorder(validacion_biomarcadores, condicion, -score_consistencia_norm)


# ============================================================
# PARTE VIII. RANKING GLOBAL
# ============================================================

ranking_global_biomarcadores <- validacion_biomarcadores[
  ,
  .(
    score_max = suppressWarnings(max(score_consistencia_norm, na.rm = TRUE)),
    score_mean = suppressWarnings(mean(score_consistencia_norm, na.rm = TRUE)),
    n_condiciones_score_alto = sum(score_consistencia_norm >= 0.70, na.rm = TRUE),
    condicion_principal = condicion[which.max(score_consistencia_norm)],
    auc_uni_max = suppressWarnings(max(auc_uni, na.rm = TRUE)),
    abs_beta_uni_max = suppressWarnings(max(abs_beta_median_uni, na.rm = TRUE)),
    abs_beta_multi_max = suppressWarnings(max(abs_beta_median_multi, na.rm = TRUE)),
    p_value_uni_min = suppressWarnings(min(p_value_uni, na.rm = TRUE))
  ),
  by = variable
]

ranking_global_biomarcadores[
  is.infinite(score_max), score_max := NA_real_
]

ranking_global_biomarcadores[
  is.infinite(score_mean), score_mean := NA_real_
]

ranking_global_biomarcadores[
  is.infinite(auc_uni_max), auc_uni_max := NA_real_
]

ranking_global_biomarcadores[
  is.infinite(abs_beta_uni_max), abs_beta_uni_max := NA_real_
]

ranking_global_biomarcadores[
  is.infinite(abs_beta_multi_max), abs_beta_multi_max := NA_real_
]

ranking_global_biomarcadores[
  is.infinite(p_value_uni_min), p_value_uni_min := NA_real_
]

setorder(ranking_global_biomarcadores, -score_max, -score_mean)


# ============================================================
# PARTE IX. TOP POR CONDICIÓN
# ============================================================

top_por_condicion <- validacion_biomarcadores[
  order(condicion, -score_consistencia_norm)
][
  ,
  head(.SD, 10),
  by = condicion
]


# ============================================================
# PARTE X. RESULTADOS FINALES
# ============================================================

cat("\n==============================\n")
cat("Variables usadas en el panel:\n")
cat("==============================\n")
print(vars_usar)

cat("\n==============================\n")
cat("Distribución de grupos:\n")
cat("==============================\n")
print(table(dt_z_final$grupo))

cat("\n==============================\n")
cat("Métricas aparentes del panel bayesiano multivariado:\n")
cat("==============================\n")
print(metricas_panel_dt)

cat("\n==============================\n")
cat("Top por condición:\n")
cat("==============================\n")
print(top_por_condicion)

cat("\n==============================\n")
cat("Ranking global de biomarcadores:\n")
cat("==============================\n")
print(ranking_global_biomarcadores)

# Objetos finales disponibles:
# vars_sin_integrar_wide
# tx_annotation_used
# dt_z_final
# metricas_panel_dt
# resultados_multi_dt
# resultados_uni_dt
# validacion_biomarcadores
# ranking_global_biomarcadores
# top_por_condicion