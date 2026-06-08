# ============================================================
# TABLA DE COHERENCIA DINÁMICA DEL PANEL FINAL
# Cruza:
# - selección final del modelo: sin_factores_65$features
# - validación interna: ranking_global_biomarcadores
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. Cargar / preparar objetos
# ------------------------------------------------------------
sin_factores_65 <- readRDS("./biomarcadores_65.rds")[[2]]
ranking_global_biomarcadores <- readRDS("./ranking_global_biomarcadores.rds")
# Si ya tienes sin_factores_65 cargado como lista:
features_final <- as.data.table(sin_factores_65)

# Si en cambio tienes solo el CSV, usa esto:
# features_final <- fread("features_65.csv")

ranking_internal <- as.data.table(ranking_global_biomarcadores)

# ------------------------------------------------------------
# 2. Variables finales a resumir
# ------------------------------------------------------------

vars_finales <- ranking_internal$variable

# ------------------------------------------------------------
# 3. Función para extraer clases separadas por ";"
# ------------------------------------------------------------

split_classes <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  
  if (length(x) == 0) return(character(0))
  
  out <- unlist(strsplit(x, ";", fixed = TRUE))
  out <- trimws(out)
  out <- out[!is.na(out) & out != "" & out != "NA"]
  out <- out[out %in% c("NP", "FD", "SO")]
  
  out
}

# ------------------------------------------------------------
# 4. Resumir clase esperada según selección final
# ------------------------------------------------------------

get_expected_class <- function(row) {
  
  clases <- character(0)
  
  # sPLS-DA
  if ("crit_splsda" %in% names(row) && !is.na(row$crit_splsda) && row$crit_splsda == 1) {
    clases <- c(clases, split_classes(row$Class_sPLSDA))
  }
  
  # ElasticNet / glmnet
  if ("crit_elastic" %in% names(row) && !is.na(row$crit_elastic) && row$crit_elastic == 1) {
    clases <- c(clases, split_classes(row$Class_elastic))
  }
  
  # PLS-DA clásico
  if ("crit_plsda" %in% names(row) && !is.na(row$crit_plsda) && row$crit_plsda == 1) {
    clases <- c(clases, split_classes(row$PLSDA_Class))
  }
  
  # Multinomial
  if ("crit_multinom" %in% names(row) && !is.na(row$crit_multinom) && row$crit_multinom == 1) {
    clases <- c(clases, split_classes(row$Multinom_Class))
  }
  
  # Fallback: si no hay criterios, usar cualquier clase disponible
  if (length(clases) == 0) {
    clases <- c(
      split_classes(row$Class_sPLSDA),
      split_classes(row$Class_elastic),
      split_classes(row$PLSDA_Class),
      split_classes(row$Multinom_Class)
    )
  }
  
  clases <- clases[clases %in% c("NP", "FD", "SO")]
  
  if (length(clases) == 0) {
    return(NA_character_)
  }
  
  tab <- sort(table(clases), decreasing = TRUE)
  top_n <- max(tab)
  top_classes <- names(tab)[tab == top_n]
  
  paste(top_classes, collapse = " / ")
}

# ------------------------------------------------------------
# 5. Aplicar a features_final
# ------------------------------------------------------------

features_sub <- features_final[
  Feature %in% vars_finales
]

features_sub[
  ,
  clase_esperada_seleccion := vapply(
    seq_len(.N),
    function(i) get_expected_class(.SD[i]),
    character(1)
  )
]

# ------------------------------------------------------------
# 6. Unir con ranking de validación interna
# ------------------------------------------------------------

tabla_coherencia <- merge(
  ranking_internal,
  features_sub[, .(
    Feature,
    Score,
    Candidate,
    Freq_sPLSDA,
    Class_sPLSDA,
    Class_elastic,
    PLSDA_Class,
    Multinom_Class,
    crit_splsda,
    crit_elastic,
    crit_plsda,
    crit_multinom,
    clase_esperada_seleccion
  )],
  by.x = "variable",
  by.y = "Feature",
  all.x = TRUE
)

# ------------------------------------------------------------
# 7. Clasificación de coherencia dinámica
# ------------------------------------------------------------

tabla_coherencia[
  ,
  tipo_coherencia := fifelse(
    is.na(clase_esperada_seleccion) | is.na(condicion_principal),
    "No evaluable",
    fifelse(
      clase_esperada_seleccion == condicion_principal,
      "Concordancia fuerte",
      fifelse(
        grepl(condicion_principal, clase_esperada_seleccion, fixed = TRUE),
        "Concordancia parcial / puente",
        fifelse(
          score_max >= 0.70,
          "Biomarcador puente",
          fifelse(
            score_max >= 0.50,
            "Discordancia dinámica moderada",
            "Señal débil o inestable"
          )
        )
      )
    )
  )
]

# ------------------------------------------------------------
# 8. Lectura dinámica automática
# ------------------------------------------------------------

make_lectura <- function(biomarcador,
                         clase_esperada,
                         condicion_dominante,
                         score_max,
                         score_mean,
                         tipo) {
  
  if (is.na(clase_esperada) || is.na(condicion_dominante)) {
    return("No evaluable: falta clase esperada o condición dominante.")
  }
  
  score_txt <- sprintf("%.3f", score_max)
  
  if (tipo == "Concordancia fuerte") {
    
    if (condicion_dominante == "NP") {
      return(paste0(
        "Concordancia fuerte. Biomarcador consistente para NP; mantiene coherencia entre la selección final y la validación interna. Score máximo = ",
        score_txt, "."
      ))
    }
    
    if (condicion_dominante == "FD") {
      return(paste0(
        "Concordancia fuerte. Biomarcador consistente para FD; mantiene señal coherente entre selección final y validación interna. Score máximo = ",
        score_txt, "."
      ))
    }
    
    if (condicion_dominante == "SO") {
      return(paste0(
        "Concordancia fuerte. Biomarcador consistente para SO; señal estable en la validación interna. Score máximo = ",
        score_txt, "."
      ))
    }
  }
  
  if (tipo == "Concordancia parcial / puente") {
    return(paste0(
      "Concordancia parcial. En la selección final aparece asociado a ",
      clase_esperada,
      " y en la validación interna domina ",
      condicion_dominante,
      ". Puede funcionar como biomarcador puente o de contraste entre condiciones. Score máximo = ",
      score_txt, "."
    ))
  }
  
  if (tipo == "Biomarcador puente") {
    return(paste0(
      "Biomarcador puente. La clase esperada fue ",
      clase_esperada,
      ", pero la coherencia interna se desplaza hacia ",
      condicion_dominante,
      ". No necesariamente debe descartarse; puede capturar contraste o transición entre grupos. Score máximo = ",
      score_txt, "."
    ))
  }
  
  if (tipo == "Discordancia dinámica moderada") {
    return(paste0(
      "Discordancia dinámica moderada. La selección final lo asocia a ",
      clase_esperada,
      ", pero la validación interna favorece ",
      condicion_dominante,
      ". Su señal puede ser contextual, redundante o dependiente del panel. Score máximo = ",
      score_txt, "."
    ))
  }
  
  if (tipo == "Señal débil o inestable") {
    return(paste0(
      "Señal débil o inestable. Aunque fue seleccionado como ",
      clase_esperada,
      ", la coherencia interna es baja y domina ",
      condicion_dominante,
      ". Interpretar como biomarcador secundario. Score máximo = ",
      score_txt, "."
    ))
  }
  
  "Lectura no clasificada."
}

tabla_coherencia[
  ,
  lectura_dinamica := mapply(
    make_lectura,
    biomarcador = variable,
    clase_esperada = clase_esperada_seleccion,
    condicion_dominante = condicion_principal,
    score_max = score_max,
    score_mean = score_mean,
    tipo = tipo_coherencia,
    SIMPLIFY = TRUE
  )
]

# ------------------------------------------------------------
# 9. Tabla final interpretable
# ------------------------------------------------------------

tabla_resumen_dinamica <- tabla_coherencia[
  ,
  .(
    Biomarcador = variable,
    Clase_esperada_segun_seleccion_final = clase_esperada_seleccion,
    Condicion_dominante_en_coherencia_interna = condicion_principal,
    Score_coherencia_max = round(score_max, 3),
    Score_coherencia_mean = round(score_mean, 3),
    Condiciones_con_score_alto = n_condiciones_score_alto,
    AUC_univariante_max = round(auc_uni_max, 3),
    Beta_univariante_abs_max = round(abs_beta_uni_max, 3),
    Beta_multivariado_abs_max = round(abs_beta_multi_max, 3),
    p_value_univariante_min = signif(p_value_uni_min, 3),
    Score_seleccion_modelo = Score,
    Candidate = Candidate,
    Tipo_coherencia = tipo_coherencia,
    Lectura_dinamica = lectura_dinamica
  )
]

setorder(tabla_resumen_dinamica, -Score_coherencia_max, -Score_coherencia_mean)

tabla_resumen_dinamica[]