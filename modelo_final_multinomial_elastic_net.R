# ============================================================
# MODELO FINAL: MULTINOMIAL CLÁSICO
# Panel concordante
# Mantiene partición original train/test de ready_for_modeling.rds
#
# No guarda archivos
# ============================================================

library(data.table)
library(nnet)

# ============================================================
# 0. Parámetros
# ============================================================

set.seed(123)

model_path <- "./ready_for_modeling.rds"

# Usar núcleo concordante:
FEATURE_THRESHOLD <- 0.70

# Si quieres usar todas las variables con score > 0:
# FEATURE_THRESHOLD <- 0

# ============================================================
# 1. Cargar modelo base
# ============================================================

model <- readRDS(model_path)

# ============================================================
# 2. Función para limpiar nombres
# ============================================================

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

# ============================================================
# 3. Metadata transcriptómica
# ============================================================

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

# ============================================================
# 4. Reconstruir valores originales SPLIT train/test
# ============================================================

unscale_block_to_wide_split <- function(block, block_name, tx_meta_small = NULL) {
  
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
  
  train_wide <- as.data.table(t(train_original), keep.rownames = "id")
  test_wide  <- as.data.table(t(test_original), keep.rownames = "id")
  
  old_names <- setdiff(names(train_wide), "id")
  
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
  
  setnames(train_wide, old_names, new_names)
  setnames(test_wide, old_names, new_names)
  
  list(
    train = train_wide[],
    test = test_wide[]
  )
}

tx_split <- unscale_block_to_wide_split(model$tx, "tx", tx_meta_small)
pr_split <- unscale_block_to_wide_split(model$pr, "pr")
me_split <- unscale_block_to_wide_split(model$me, "me")
cl_split <- unscale_block_to_wide_split(model$cl, "cl")

# ============================================================
# 5. Resolver duplicados globales igual en train/test
# ============================================================

all_block_names <- list(
  tx = setdiff(names(tx_split$train), "id"),
  pr = setdiff(names(pr_split$train), "id"),
  me = setdiff(names(me_split$train), "id"),
  cl = setdiff(names(cl_split$train), "id")
)

all_names <- unlist(all_block_names, use.names = FALSE)
all_names_unique <- make.unique(all_names, sep = "_")

rename_block_global <- function(split_obj, block_names, all_names_unique, idx_start) {
  
  n <- length(block_names)
  new_names <- all_names_unique[idx_start:(idx_start + n - 1)]
  
  setnames(split_obj$train, block_names, new_names)
  setnames(split_obj$test, block_names, new_names)
  
  list(
    split = split_obj,
    next_idx = idx_start + n
  )
}

idx <- 1

tmp <- rename_block_global(tx_split, all_block_names$tx, all_names_unique, idx)
tx_split <- tmp$split
idx <- tmp$next_idx

tmp <- rename_block_global(pr_split, all_block_names$pr, all_names_unique, idx)
pr_split <- tmp$split
idx <- tmp$next_idx

tmp <- rename_block_global(me_split, all_block_names$me, all_names_unique, idx)
me_split <- tmp$split
idx <- tmp$next_idx

tmp <- rename_block_global(cl_split, all_block_names$cl, all_names_unique, idx)
cl_split <- tmp$split
idx <- tmp$next_idx

# ============================================================
# 6. Unir bloques por separado en train/test
# ============================================================

vars_train_wide <- Reduce(
  function(x, y) merge(x, y, by = "id", all = TRUE),
  list(
    tx_split$train,
    pr_split$train,
    me_split$train,
    cl_split$train
  )
)

vars_test_wide <- Reduce(
  function(x, y) merge(x, y, by = "id", all = TRUE),
  list(
    tx_split$test,
    pr_split$test,
    me_split$test,
    cl_split$test
  )
)

# ============================================================
# 7. Añadir grupo manteniendo partición original
# ============================================================

grupo_train <- as.data.table(model$grupo_train)
grupo_test  <- as.data.table(model$grupo_test)

grupo_train[, grupo := factor(grupo, levels = c("NP", "FD", "SO"))]
grupo_test[, grupo := factor(grupo, levels = c("NP", "FD", "SO"))]

train_dt <- merge(
  grupo_train,
  vars_train_wide,
  by = "id",
  all.x = TRUE
)

test_dt <- merge(
  grupo_test,
  vars_test_wide,
  by = "id",
  all.x = TRUE
)

setcolorder(train_dt, c("id", "grupo"))
setcolorder(test_dt, c("id", "grupo"))

setorder(train_dt, grupo, id)
setorder(test_dt, grupo, id)

cat("\n==============================\n")
cat("Partición original\n")
cat("==============================\n")
print(dim(train_dt))
print(dim(test_dt))
print(table(train_dt$grupo))
print(table(test_dt$grupo))

# ============================================================
# 8. Seleccionar biomarcadores concordantes
# ============================================================

if (!exists("ranking_global_biomarcadores")) {
  stop("No existe ranking_global_biomarcadores en memoria.")
}

ranking_global_biomarcadores <- as.data.table(ranking_global_biomarcadores)

features_concordantes <- ranking_global_biomarcadores[
  !is.na(score_max) & score_max > FEATURE_THRESHOLD,
  variable
]

features_concordantes <- unique(features_concordantes)

features_presentes <- intersect(features_concordantes, names(train_dt))
features_faltantes <- setdiff(features_concordantes, names(train_dt))

if (length(features_faltantes) > 0) {
  warning(
    "Estas features no están en train/test: ",
    paste(features_faltantes, collapse = ", ")
  )
}

if (length(features_presentes) < 2) {
  stop("Hay menos de 2 features concordantes presentes.")
}

cat("\n==============================\n")
cat("Features usadas por el multinomial clásico\n")
cat("==============================\n")
print(features_presentes)

# ============================================================
# 9. Preparar matrices con z-score usando SOLO train
# ============================================================

X_train_raw <- as.matrix(train_dt[, ..features_presentes])
X_test_raw  <- as.matrix(test_dt[, ..features_presentes])

mode(X_train_raw) <- "numeric"
mode(X_test_raw)  <- "numeric"

y_train <- factor(train_dt$grupo, levels = c("NP", "FD", "SO"))
y_test  <- factor(test_dt$grupo, levels = c("NP", "FD", "SO"))

# ------------------------------------------------------------
# Imputación con mediana del train
# ------------------------------------------------------------

train_medians <- apply(
  X_train_raw,
  2,
  function(x) median(x, na.rm = TRUE)
)

for (j in seq_along(train_medians)) {
  
  if (is.na(train_medians[j])) {
    train_medians[j] <- 0
  }
  
  X_train_raw[is.na(X_train_raw[, j]), j] <- train_medians[j]
  X_test_raw[is.na(X_test_raw[, j]), j] <- train_medians[j]
}

# ------------------------------------------------------------
# Escalado train-only
# ------------------------------------------------------------

train_means <- colMeans(X_train_raw, na.rm = TRUE)
train_sds <- apply(X_train_raw, 2, sd, na.rm = TRUE)

valid_sd <- !is.na(train_sds) & train_sds > 0

X_train_raw <- X_train_raw[, valid_sd, drop = FALSE]
X_test_raw  <- X_test_raw[, valid_sd, drop = FALSE]

train_means <- train_means[valid_sd]
train_sds <- train_sds[valid_sd]

features_finales_modelo <- colnames(X_train_raw)

X_train <- scale(
  X_train_raw,
  center = train_means,
  scale = train_sds
)

X_test <- scale(
  X_test_raw,
  center = train_means,
  scale = train_sds
)

X_train <- as.matrix(X_train)
X_test <- as.matrix(X_test)

cat("\n==============================\n")
cat("Matriz final de modelado\n")
cat("==============================\n")
print(dim(X_train))
print(dim(X_test))
print(features_finales_modelo)

# ============================================================
# 10. Funciones métricas
# ============================================================

calc_auc_binary <- function(y01, score) {
  
  ok <- complete.cases(y01, score)
  y01 <- as.numeric(y01[ok])
  score <- as.numeric(score[ok])
  
  if (length(unique(y01)) < 2) {
    return(NA_real_)
  }
  
  r <- rank(score)
  n1 <- sum(y01 == 1)
  n0 <- sum(y01 == 0)
  
  auc <- (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  
  auc
}

safe_logloss_multiclass <- function(y_true, prob_mat, levels_y) {
  
  y_true <- factor(y_true, levels = levels_y)
  
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, levels_y, drop = FALSE]
  prob_mat <- pmin(pmax(prob_mat, 1e-8), 1 - 1e-8)
  
  idx <- cbind(seq_along(y_true), as.integer(y_true))
  
  -mean(log(prob_mat[idx]))
}

calc_metrics_multiclass <- function(y_true, pred_class, prob_mat, model_name) {
  
  levels_y <- levels(y_true)
  
  y_true <- factor(y_true, levels = levels_y)
  pred_class <- factor(pred_class, levels = levels_y)
  
  cm <- table(
    Truth = y_true,
    Pred = pred_class
  )
  
  acc <- mean(y_true == pred_class)
  
  recall_by_class <- diag(cm) / rowSums(cm)
  bal_acc <- mean(recall_by_class, na.rm = TRUE)
  ber <- 1 - bal_acc
  
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, levels_y, drop = FALSE]
  
  auc_by_class <- rbindlist(
    lapply(levels_y, function(cc) {
      
      y01 <- as.integer(y_true == cc)
      auc <- calc_auc_binary(y01, prob_mat[, cc])
      
      data.table(
        Model = model_name,
        Class = cc,
        AUC = auc
      )
    })
  )
  
  macro_auc <- mean(auc_by_class$AUC, na.rm = TRUE)
  logloss <- safe_logloss_multiclass(y_true, prob_mat, levels_y)
  
  metrics <- data.table(
    Model = model_name,
    Accuracy = acc,
    BalAccuracy = bal_acc,
    BER = ber,
    MacroAUC = macro_auc,
    LogLoss = logloss,
    n_train = length(y_train),
    n_test = length(y_true),
    n_features = ncol(X_train)
  )
  
  list(
    metrics = metrics,
    auc_by_class = auc_by_class,
    confusion = as.data.table(cm)
  )
}

# ============================================================
# 11. Ajustar multinomial clásico
# ============================================================

df_train_multinom <- data.frame(
  grupo = y_train,
  X_train,
  check.names = FALSE
)

df_test_multinom <- data.frame(
  X_test,
  check.names = FALSE
)

fit_multinom_final <- nnet::multinom(
  grupo ~ .,
  data = df_train_multinom,
  trace = FALSE,
  maxit = 1000
)

# ============================================================
# 12. Predicción en test
# ============================================================

prob_multinom <- predict(
  fit_multinom_final,
  newdata = df_test_multinom,
  type = "probs"
)

prob_multinom <- as.matrix(prob_multinom)
prob_multinom <- prob_multinom[, levels(y_train), drop = FALSE]

pred_multinom <- colnames(prob_multinom)[max.col(prob_multinom)]
pred_multinom <- factor(pred_multinom, levels = levels(y_train))

res_multinom_final <- calc_metrics_multiclass(
  y_true = y_test,
  pred_class = pred_multinom,
  prob_mat = prob_multinom,
  model_name = "Classic_multinomial"
)

metrics_multinom_final <- res_multinom_final$metrics
auc_multinom_final <- res_multinom_final$auc_by_class
confusion_multinom_final <- res_multinom_final$confusion

confusion_multinom_final[, Model := "Classic_multinomial"]

predictions_multinom_final <- data.table(
  id = test_dt$id,
  Truth = y_test,
  Pred = pred_multinom,
  prob_multinom
)

# ============================================================
# 13. Coeficientes del multinomial clásico
# ============================================================

coef_mat <- coef(fit_multinom_final)

coef_dt <- as.data.table(
  coef_mat,
  keep.rownames = "Class"
)

coef_long <- melt(
  coef_dt,
  id.vars = "Class",
  variable.name = "Feature",
  value.name = "Coef"
)

coef_long[
  ,
  Term := fifelse(Feature == "(Intercept)", "Intercept", "Biomarker")
]

coef_long[
  ,
  AbsCoef := abs(Coef)
]

coef_long[
  Term == "Biomarker",
  Direction := fifelse(Coef > 0, "Positive", "Negative")
]

setorder(coef_long, Class, -AbsCoef)

coef_nonzero <- coef_long[
  Term == "Biomarker" & Coef != 0
]

coef_wide <- dcast(
  coef_long[
    Term == "Biomarker",
    .(Feature, Class, Coef)
  ],
  Feature ~ Class,
  value.var = "Coef"
)

coef_wide[
  ,
  MaxAbsCoef := apply(
    abs(as.matrix(.SD)),
    1,
    max,
    na.rm = TRUE
  ),
  .SDcols = setdiff(names(coef_wide), "Feature")
]

setorder(coef_wide, -MaxAbsCoef)

coef_importance <- coef_long[
  Term == "Biomarker",
  .(
    MaxAbsCoef = max(AbsCoef, na.rm = TRUE),
    MeanAbsCoef = mean(AbsCoef, na.rm = TRUE),
    SumAbsCoef = sum(AbsCoef, na.rm = TRUE),
    DirectionSummary = paste0(
      Class,
      ":",
      ifelse(Coef > 0, "+", "-"),
      collapse = "; "
    )
  ),
  by = Feature
]

setorder(coef_importance, -MaxAbsCoef)

# ============================================================
# 14. Anotar coeficientes con coherencia interna
# ============================================================

coef_importance_annot <- merge(
  coef_importance,
  ranking_global_biomarcadores[
    ,
    .(
      Feature = variable,
      score_max,
      score_mean,
      condicion_principal,
      auc_uni_max,
      p_value_uni_min
    )
  ],
  by = "Feature",
  all.x = TRUE
)

setorder(coef_importance_annot, -MaxAbsCoef)

# ============================================================
# 15. Matriz de confusión y métricas por clase
# ============================================================

conf_mat <- table(
  Truth = predictions_multinom_final$Truth,
  Predicted = predictions_multinom_final$Pred
)

conf_mat_dt <- as.data.table(conf_mat)
setnames(conf_mat_dt, c("Truth", "Predicted", "N"))

conf_mat_pct <- prop.table(conf_mat, margin = 1) * 100
conf_mat_pct_dt <- as.data.table(conf_mat_pct)
setnames(conf_mat_pct_dt, c("Truth", "Predicted", "Percent"))

conf_mat_plot_dt <- merge(
  conf_mat_dt,
  conf_mat_pct_dt,
  by = c("Truth", "Predicted"),
  all = TRUE
)

conf_mat_plot_dt[
  ,
  label := paste0(N, "\n", round(Percent, 1), "%")
]

classes <- levels(y_test)

class_metrics <- rbindlist(lapply(classes, function(cl) {
  
  TP <- conf_mat[cl, cl]
  FN <- sum(conf_mat[cl, ]) - TP
  FP <- sum(conf_mat[, cl]) - TP
  TN <- sum(conf_mat) - TP - FN - FP
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  f1 <- ifelse(
    !is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0,
    2 * precision * sensitivity / (precision + sensitivity),
    NA_real_
  )
  
  data.table(
    Class = cl,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1
  )
}))

# ============================================================
# 16. Plot matriz de confusión
# ============================================================

library(ggplot2)

p_confusion_multinom <- ggplot(
  conf_mat_plot_dt,
  aes(x = Predicted, y = Truth, fill = N)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = label),
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "white",
    high = "gray40"
  ) +
  labs(
    title = "Confusion matrix - Classic multinomial model",
    subtitle = "Counts and row percentages",
    x = "Predicted class",
    y = "True class",
    fill = "N"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

# ============================================================
# 17. Imprimir resultados finales
# ============================================================

cat("\n============================================================\n")
cat("MODELO FINAL: MULTINOMIAL CLÁSICO\n")
cat("============================================================\n")

cat("\nFeatures finales:\n")
print(features_finales_modelo)

cat("\nDistribución train:\n")
print(table(y_train))

cat("\nDistribución test:\n")
print(table(y_test))

cat("\nMétricas finales en test:\n")
print(metrics_multinom_final)

cat("\nAUC one-vs-rest por clase:\n")
print(auc_multinom_final)

cat("\nMatriz de confusión:\n")
print(confusion_multinom_final)

cat("\nMétricas por clase desde matriz de confusión:\n")
print(class_metrics)

cat("\nPredicciones en test:\n")
print(predictions_multinom_final)

cat("\nCoeficientes no cero por clase:\n")
print(coef_nonzero)

cat("\nCoeficientes formato ancho:\n")
print(coef_wide)

cat("\nImportancia global por coeficientes anotada:\n")
print(coef_importance_annot)

print(p_confusion_multinom)

# ============================================================
# 18. Objetos finales disponibles
# ============================================================

# train_dt
# test_dt
# X_train
# X_test
# y_train
# y_test
# features_finales_modelo
# fit_multinom_final
# metrics_multinom_final
# auc_multinom_final
# confusion_multinom_final
# predictions_multinom_final
# coef_long
# coef_nonzero
# coef_wide
# coef_importance
# coef_importance_annot
# class_metrics
# conf_mat_plot_dt
# p_confusion_multinom