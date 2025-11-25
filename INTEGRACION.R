library(MOFA2)

## ------------------------------------------------------------
## 0) Helpers
## ------------------------------------------------------------

.get_train_ids <- function(rfm) {
  for (nm in c("cl","tx","pr","me")) {
    if (!is.null(rfm[[nm]]) && !is.null(rfm[[nm]]$train)) {
      ids <- colnames(rfm[[nm]]$train)
      if (is.null(ids)) stop(paste0("Faltan colnames en ", nm, "$train"))
      return(ids)
    }
  }
  stop("No hay vistas con matriz $train y colnames definidos.")
}

.stack_view_train <- function(view) {
  Xtr <- view$train
  if (is.null(colnames(Xtr))) stop("La matriz train no tiene colnames (sujetos).")
  return(Xtr)
}

# Construcción de metadata usando grupo_df$train
.build_metadata <- function(rfm, extra_cols = NULL) {
  sample_ids <- .get_train_ids(rfm)
  
  if (!is.null(rfm$grupo_df) && !is.null(rfm$grupo_df$train)) {
    gdf <- as.data.frame(rfm$grupo_df$train, stringsAsFactors = FALSE)
    stopifnot(all(c("id","grupo") %in% colnames(gdf)))
    # asegurar alineación
    rownames(gdf) <- gdf$id
    gdf <- gdf[sample_ids, , drop=FALSE]
    fac <- factor(gdf$grupo)
    names(fac) <- gdf$id
  } else {
    # fallback: prefijo
    fac <- factor(sub("^([^_]+)_.*$", "\\1", sample_ids))
    names(fac) <- sample_ids
  }
  
  meta <- data.frame(
    sample = sample_ids,
    grupo  = fac,
    stringsAsFactors = FALSE,
    row.names = sample_ids
  )
  
  if (!is.null(extra_cols)) {
    for (nm in names(extra_cols)) {
      v <- extra_cols[[nm]]
      if (!is.null(names(v))) v <- v[sample_ids]
      else if (length(v) == length(sample_ids)) names(v) <- sample_ids
      else stop(paste0("extra_metadata '", nm, "' no alinea con muestras."))
      meta[[nm]] <- v
    }
  }
  
  colnames(meta) <- gsub("[^[:alnum:]]", "_", colnames(meta))
  colnames(meta) <- gsub("_+", "_", colnames(meta))
  colnames(meta) <- sub("_$", "", colnames(meta))
  return(meta)
}

## ------------------------------------------------------------
## 1) Función central intacta
## ------------------------------------------------------------
mofa_componentes_original <- function(ncomp,
                                      semilla,
                                      directorio_modelo,
                                      mofa.obj,
                                      sample_metadata_) {
  samples_metadata(mofa.obj) <- sample_metadata_
  
  if (!dir.exists(directorio_modelo)) {
    dir.create(directorio_modelo, recursive = TRUE)
  }
  
  data_opts <- get_default_data_options(mofa.obj)
  data_opts$scale_views <- FALSE
  
  model_opts <- get_default_model_options(mofa.obj)
  model_opts$num_factors <- ncomp
  
  train_opts <- get_default_training_options(mofa.obj)
  train_opts$seed <- semilla
  train_opts$convergence_mode <- "slow"
  
  MOFAobject <- prepare_mofa(
    object = mofa.obj,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile <- file.path(directorio_modelo, paste0("modelo_", ncomp, ".hdf5"))
  if (file.exists(outfile)) {
    message("Eliminando archivo previo: ", outfile)
    file.remove(outfile)
  }
  gc(); reticulate::py_run_string("import gc; gc.collect()")
  
  run_mofa(MOFAobject, outfile,
           use_basilisk = TRUE,
           save_data = FALSE)
}

## ------------------------------------------------------------
## 2) Construcción de vistas SOLO con train
## ------------------------------------------------------------
build_mofa_from_rfm <- function(read_for_modeling) {
  vistas_raw <- list(
    transcriptomica = read_for_modeling$tx,
    proteomica      = read_for_modeling$pr,
    metabolomica    = read_for_modeling$me,
    clinical        = read_for_modeling$cl
  )
  vistas_raw <- vistas_raw[!vapply(vistas_raw, is.null, logical(1))]
  
  omicas <- lapply(vistas_raw, .stack_view_train)
  
  omicas_MOFA <- create_mofa_from_matrix(omicas)
  list(omicas = omicas, mofa = omicas_MOFA)
}

## ------------------------------------------------------------
## 3) Pipeline
## ------------------------------------------------------------
run_mofa_grid_from_rfm <- function(read_for_modeling,
                                   factores = 2:6,
                                   semilla = 23355,
                                   directorio_modelo = "./modelos",
                                   extra_metadata = NULL,
                                   seleccionar = TRUE,
                                   plot = TRUE,
                                   save_dir = "../../scripts_Data") {
  
  built <- build_mofa_from_rfm(read_for_modeling)
  omicas_MOFA <- built$mofa
  
  aux <- .build_metadata(read_for_modeling, extra_cols = extra_metadata)
  
  # check
  common_ids <- colnames(built$omicas[[1]])
  if (!identical(rownames(aux), common_ids))
    stop("rownames(metadata) deben coincidir con colnames de las ómicas (train).")
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  saveRDS(built$omicas, file.path(save_dir, "omicas_train.rds"))
  
  lista_modelos <- lapply(factores, function(x) {
    mofa_componentes_original(
      ncomp = x,
      semilla = semilla,
      directorio_modelo = directorio_modelo,
      mofa.obj = omicas_MOFA,
      sample_metadata_ = aux
    )
  })
  
  MOFA2::compare_elbo(lista_modelos)
  modelo <- if (seleccionar) select_model(lista_modelos, plot = plot) else lista_modelos[[length(lista_modelos)]]
  
  if (plot) {
    plot_factor_cor(modelo)
    plot_variance_explained(modelo)
    plot_factors(modelo, color_by = "grupo")
  }
  
  saveRDS(modelo, file.path(save_dir, "modelo.rds"))
  
  invisible(list(
    modelos = lista_modelos,
    modelo_seleccionado = modelo,
    metadata = samples_metadata(modelo)
  ))
}


# ------------------------------------------------------------
# 4) Uso
# ------------------------------------------------------------
rfm <- readRDS("./ready_for_modeling.rds")
out <- run_mofa_grid_from_rfm(
  rfm,
  factores = 2:4,
  semilla = 23355,
  directorio_modelo = "./modelos",
  extra_metadata = NULL,
  seleccionar = TRUE,
  plot = TRUE,
  save_dir = "modelos"
)


modelos <- list.files("./modelos", full.names = TRUE,pattern = "hdf5")
modelo  <- lapply(modelos, load_model)
plot_models <- MOFA2::compare_elbo(modelo)  # gráfico ELBO

modelo_final <- select_model(modelo)        # asume función disponible (MOFA2 >= 1.4)


saveRDS(modelo_final,"./modelo.rds")

