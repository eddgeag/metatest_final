# Librerias======
library(reticulate)
use_condaenv("r-tf", required = TRUE)
library(tensorflow)
# tf$config$experimental$set_memory_growth(tf$config$list_physical_devices('GPU')[[1]], TRUE)  # Si usas GPU
# options(tensorflow.memory.allocator = "gpu_allocator")  # Opcional para GPU
library(keras3)
library(caret)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(glmnet)
library(dplyr)
library(purrr)
library(caret)
library(ggplot2)
library(MOFA2)
library(ParamHelpers)
library(mlrMBO)
library(DiceKriging)
library(ParamHelpers)





options(tensorflow.extract.disallow_out_of_bounds = FALSE)
tf$debugging$set_log_device_placement(TRUE)






modelo <- readRDS("./modelos/modelo.rds")

# Normalización
minmax_scale <- function(X) {
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  X_scaled <- sweep(X, 2, mins, FUN = "-")
  X_scaled <- sweep(X_scaled, 2, maxs - mins, FUN = "/")
  attr(X_scaled, "min") <- mins
  attr(X_scaled, "max") <- maxs
  return(X_scaled)
}
minmax_descale <- function(X_scaled, mins, maxs) {
  X_unscaled <- sweep(X_scaled, 2, maxs - mins, FUN = "*")
  X_unscaled <- sweep(X_unscaled, 2, mins, FUN = "+")
  return(X_unscaled)
}
clasifica_lista_omicas <- function(vec) {
  lista <- list(
    transcriptoma = vec[grepl("^X[0-9]+$", vec)],
    metaboloma    = vec[grepl("^(pos_|neg_)", vec)],
    proteoma      = vec[!(grepl("^X[0-9]+$", vec) |
                            grepl("^(pos_|neg_)", vec))]
    
  )
  return(lista)
}

# 1. SamplingLayer - KL en add_loss
SamplingLayer <- new_layer_class(
  classname = "SamplingLayer",
  initialize = function(self) {
    super()$`__init__`()
  },
  call = function(self, inputs) {
    z_mean <- inputs[[1]]
    z_log_var <- inputs[[2]]
    epsilon <- tf$random$normal(shape = tf$shape(z_mean))
    z <- z_mean + keras$ops$exp(0.5 * z_log_var) * epsilon
    z  # Solo retorna z, sin añadir kl_loss
  }
)
CustomLossLayer <- keras::new_layer_class(
  classname = "CustomLossLayer",
  
  initialize = function(self, params, ...) {
    tf <- tensorflow::tf
    k  <- tf$keras$backend
    super$initialize(...)
    
    self$params     <- params
    self$latent_dim <- as.integer(params$latent_dim)
    self$proj_dim   <- as.integer(params$proj_dim %||% 2)  # Default a 6, puedes cambiarlo al crear la capa
    
    # Matriz R en el espacio proyectado
    self$R <- tf$Variable(tf$eye(self$proj_dim, dtype = "float32"),
                          trainable = FALSE,
                          name = "R_var")
    # Pesos de pérdida
    self$params$mse_weight        <- params$mse_weight        %||% 1.0
    self$kl_weight_var            <- k$variable(
      as.numeric(params$kl_weight %||% params$beta %||% 0.001),
      dtype = "float32",
      name  = "kl_weight_var"
    )
    self$params$procrustes_weight <- params$procrustes_weight %||% 0.5
    self$params$ortho_weight      <- params$ortho_weight      %||% 0.01
    # Métricas
    self$mse_metric        <- keras::metric_mean(name = "mse_loss")
    self$kl_metric         <- keras::metric_mean(name = "kl_loss")
    self$procrustes_metric <- keras::metric_mean(name = "procrustes_loss")
    self$ortho_metric      <- keras::metric_mean(name = "ortho_loss")
    self$current_epoch <- k$variable(0L, dtype = "int32", name = "current_epoch")
    self$projection <- layer_dense(
      units = self$proj_dim,
      # AHORA parámetro general
      activation = NULL,
      use_bias = FALSE,
      kernel_initializer = "orthogonal"
    )
  },
  
  call = function(self, inputs) {
    tf <- tensorflow::tf
    k  <- tf$keras$backend
    
    # Inputs
    z_mean       <- inputs[[1]]             # (batch, latent_dim)
    z_log_var    <- inputs[[2]]
    y_true       <- inputs[[3]]
    y_pred       <- inputs[[4]]
    z_orig_batch <- inputs[[5]]             # (batch, proj_dim)
    
    projected_z <- self$projection(z_mean)  # (batch, proj_dim)
    
    # 1) MSE + KL
    mse_loss <- k$mean(k$square(y_true - y_pred))
    kl_loss  <- -0.5 * k$mean(1 + z_log_var - k$square(z_mean) - k$exp(z_log_var))
    
    # 2) Procrustes en espacio de proyección
    mc <- projected_z - tf$reduce_mean(projected_z, axis = 0L)
    z_o_centered <- z_orig_batch - tf$reduce_mean(z_orig_batch, axis = 0L)
    
    cov_matrix <- tf$linalg$matmul(tf$transpose(mc), z_o_centered)  # (proj_dim, proj_dim)
    svd_res    <- tf$linalg$svd(cov_matrix,
                                compute_uv = TRUE,
                                full_matrices = FALSE)
    u          <- svd_res[[2]]
    v          <- svd_res[[3]]
    newR       <- tf$linalg$matmul(u, tf$transpose(v))
    self$R$assign(newR)
    
    aligned         <- tf$linalg$matmul(mc, self$R)
    procrustes_loss <- k$mean(k$square(aligned - z_o_centered))
    
    # 3) Ortogonalidad en espacio proyectado
    ortho_loss <- self$params$ortho_weight * k$mean(k$square(
      k$dot(self$R, k$transpose(self$R)) - tf$eye(self$proj_dim, dtype = "float32")
    ))
    
    # Actualizar métricas
    self$mse_metric$update_state(mse_loss)
    self$kl_metric$update_state(kl_loss)
    self$procrustes_metric$update_state(procrustes_loss)
    self$ortho_metric$update_state(ortho_loss)
    
    # Pérdida total
    total_loss <- self$params$mse_weight        * mse_loss +
      self$kl_weight_var            * kl_loss +
      self$params$procrustes_weight * procrustes_loss +
      ortho_loss
    self$add_loss(total_loss)
    
    # Devolver la reconstrucción
    k$set_value(self$current_epoch, k$get_value(self$current_epoch) + 1L)
    y_pred
  },
  
  compute_output_shape = function(self, input_shape) {
    input_shape[[4]]
  }
)

# 1) Capa de rotación


# 2) Callback para detener si aparecen NaNs
RotationLayer <- new_layer_class(
  classname = "RotationLayer",
  initialize = function(self, rotation_matrix, ...) {
    super$initialize(...)
    self$rotation_matrix <- tf$constant(rotation_matrix, dtype = "float32")
  },
  call = function(self, inputs) {
    tf$linalg$matmul(inputs, self$rotation_matrix)
  }
)

# 2) Callback para detener si aparecen NaNs
callback_terminate_on_naan <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (any(is.na(unlist(logs)))) {
      message("NaN detected - Stopping training")
      keras3:::k_stop_training()
    }
  }
)
def_cvae_conditional <- function(params, input_dim, n_classes, z_original_) {
  # Extract scalar parameters
  latent_dim       <- as.integer(params$latent_dim)
  intermediate_dim <- as.integer(params$intermediate_dim)
  n_layers         <- as.integer(params$n_layers)
  l2_reg           <- params$l2_reg
  learning_rate    <- params$learning_rate
  epochs           <- as.integer(params$epochs)
  init_mse_weight        <- params$init_mse_weight
  init_procrustes_weight <- params$init_procrustes_weight
  final_procrustes_weight <- params$final_procrustes_weight
  phase_switch_epoch     <- as.integer(params$phase_switch_epoch)
  beta                   <- params$beta
  
  # Trim reference latent for Procrustes
  z_original_trim <- as.matrix(z_original_)
  
  # === Encoder inputs ===
  inp_x      <- layer_input(shape = input_dim, name = "features")
  inp_cls    <- layer_input(shape = n_classes, name = "class")
  inp_z_orig <- layer_input(shape = 2, name = "z_orig_batch")
  
  # Encode features + class
  h <- layer_concatenate(list(inp_x, inp_cls))
  for (i in seq_len(n_layers)) {
    h <- h %>%
      layer_dense(
        intermediate_dim,
        activation = "relu",
        kernel_regularizer = regularizer_l2(l2_reg),
        kernel_initializer = initializer_he_normal()
      ) %>%
      layer_batch_normalization()
  }
  z_mean    <- h %>% layer_dense(latent_dim,
                                 name = "z_mean",
                                 kernel_initializer = initializer_he_normal())
  z_log_var <- h %>% layer_dense(latent_dim,
                                 name = "z_log_var",
                                 kernel_initializer = initializer_he_normal())
  z         <- SamplingLayer()(list(z_mean, z_log_var))
  
  # === Decoder definition ===
  dec_z   <- layer_input(shape = latent_dim, name = "decoder_latent")
  dec_cls <- layer_input(shape = n_classes, name = "decoder_class")
  d <- layer_concatenate(list(dec_z, dec_cls))
  for (i in seq_len(n_layers)) {
    d <- d %>%
      layer_dense(
        intermediate_dim,
        activation = "relu",
        kernel_regularizer = regularizer_l2(l2_reg),
        kernel_initializer = initializer_he_normal()
      ) %>%
      layer_batch_normalization()
  }
  dec_out <- d %>% layer_dense(input_dim,
                               activation = "linear",
                               kernel_initializer = initializer_he_normal())
  decoder <- keras_model(inputs = list(dec_z, dec_cls), outputs = dec_out)
  
  # Reconstruct within CVAE graph
  cvae_recon <- decoder(list(z, inp_cls))
  
  # === Custom loss layer for Procrustes + VAE ===
  outputs <- CustomLossLayer(
    params = list(
      latent_dim              = latent_dim,
      mse_weight              = init_mse_weight,
      kl_weight               = beta,
      procrustes_weight       = init_procrustes_weight,
      final_procrustes_weight = final_procrustes_weight,
      phase_switch_epoch      = phase_switch_epoch,
      epochs                  = epochs,
      proj_dim = 2 # aquí defines la proyección que quieras
      
    ),
    name = "procrustes_loss_layer"
  )(list(z_mean, z_log_var, inp_x, cvae_recon, inp_z_orig))
  
  # === Full CVAE model ===
  cvae <- keras_model(inputs = list(inp_x, inp_cls, inp_z_orig),
                      outputs = outputs) %>% compile(
                        optimizer = optimizer_adam(learning_rate = learning_rate, clipnorm = 1.0),
                        loss      = NULL
                      )
  
  # === Encoder for downstream use ===
  encoder <- keras_model(
    inputs  = list(inp_x, inp_cls, inp_z_orig),
    outputs = list(z_mean, z_log_var, z)
  )
  
  list(model   = cvae,
       encoder = encoder,
       decoder = decoder)
}


build_and_evaluate_cvae_conditional <- function(params,
                                                Z_data,
                                                grupo,
                                                clase_levels,
                                                n_folds = 2) {
  # Asegurar reproducibilidad en TF
  options(tensorflow.extract.disallow_out_of_bounds = FALSE)
  
  
  
  
  # Helper para extraer valor escalar de params
  extract_param <- function(p) {
    if (is.list(p) && !is.null(p$values))
      p$values[[1]]
    else
      p
  }
  
  # Extraer parámetros como escalares
  latent_dim_raw   <- as.integer(extract_param(params$latent_dim))
  intermediate_dim <- as.integer(extract_param(params$intermediate_dim))
  n_layers         <- as.integer(extract_param(params$n_layers))
  learning_rate    <- extract_param(params$learning_rate)
  beta             <- extract_param(params$beta)
  l2_reg           <- extract_param(params$l2_reg)
  batch_size       <- as.integer(extract_param(params$batch_size))
  epochs           <- as.integer(extract_param(params$epochs))
  
  # Nuevos parámetros para el entrenamiento por fases (optimizados por MBO)
  init_mse_weight         <- extract_param(params$init_mse_weight)
  init_procrustes_weight  <- extract_param(params$init_procrustes_weight)
  final_procrustes_weight <- extract_param(params$final_procrustes_weight)
  phase_switch_epoch      <- as.integer(extract_param(params$phase_switch_epoch))
  
  # Preparar Z_data
  if (!is.matrix(Z_data))
    Z_data <- as.matrix(Z_data)
  storage.mode(Z_data) <- "numeric"
  
  # Validar latent_dim
  # if (latent_dim_raw > ncol(Z_data)) {
  #   stop(sprintf(
  #     "latent_dim (%d) excede columnas de Z_data (%d)",
  #     latent_dim_raw,
  #     ncol(Z_data)
  #   ))
  # }
  latent_dim <- latent_dim_raw
  
  # Recortar Z_data para Procrustes
  Z_original_subset <- Z_data
  
  # Wrapper para reconstrucción segura usando encoder+decoder
  # NUEVO
  safe_reconstruct <- function(encoder,
                               decoder,
                               x,
                               labels,
                               z_orig_batch,
                               target_shape) {
    tryCatch({
      # 1) pasa x, labels y z_orig_batch al encoder
      enc_out <- predict(encoder, list(x, labels, z_orig_batch))
      z_mean  <- enc_out[[1]]
      # 2) decodifica a partir de z_mean y labels
      recon   <- predict(decoder, list(z_mean, labels))
      if (!all(dim(recon) == target_shape)) {
        stop(sprintf(
          "Shape mismatch: got %s, expected %s",
          paste(dim(recon), collapse = "x"),
          paste(target_shape, collapse = "x")
        ))
      }
      recon
    }, error = function(e) {
      warning("Reconstrucción falló: ", e$message)
      matrix(0, nrow = target_shape[1], ncol = target_shape[2])
    })
  }
  
  
  # Configurar validación cruzada
  safe_n_folds <- function(grupo) {
    min_class_size <- min(table(grupo))
    if (min_class_size < 4) return(2)   # forzar como mínimo 2 folds
    max_folds <- floor(min_class_size / 2)
    max(2, min(n_folds, max_folds))
  }
  n_folds_actual <- safe_n_folds(grupo)
  
  
  # 3. Crear folds estratificados
  library(splitTools)
  train_folds <- create_folds(grupo, k = n_folds_actual, type = "stratified")
  val_folds <- lapply(train_folds, function(idx)
    setdiff(seq_along(grupo), idx))
  
  metrics <- list()
  input_dim <- ncol(Z_data)
  n_classes <- length(clase_levels)
  R_list <- list()
  
  for (i in seq_along(train_folds)) {
    val_idx   <- val_folds[[i]]
    train_idx <- train_folds[[i]]
    
    x_train <- Z_data[train_idx, , drop = FALSE]
    x_val   <- Z_data[val_idx, , drop = FALSE]
    f_train <- factor(grupo[train_idx], levels = clase_levels)
    f_val   <- factor(grupo[val_idx], levels = clase_levels)
    
    # Saltar pliegues sin todas las clases
    if (any(table(f_train) == 0) || any(table(f_val) == 0))
      next
    
    labels_train <- to_categorical(as.integer(f_train) - 1, num_classes = n_classes)
    labels_val   <- to_categorical(as.integer(f_val) - 1, num_classes = n_classes)
    
    # Chequeos de dimensión
    stopifnot(
      ncol(x_train) == input_dim,
      ncol(x_val) == input_dim,
      nrow(labels_train) == nrow(x_train),
      ncol(labels_train) == n_classes
    )
    
    # Dentro del bucle, justo después de calcular train_idx / val_idx…
    Z_original_train <- Z_original_subset[train_idx, , drop = FALSE]
    Z_original_val   <- Z_original_subset[val_idx, , drop = FALSE]
    procrustes_dim <- 2 # Cambia este valor si necesitas otro
    
    # Dentro del loop por fold:
    Z_original_train_proj <- Z_original_train[, 1:procrustes_dim, drop =
                                                FALSE]
    Z_original_val_proj   <- Z_original_val[, 1:procrustes_dim, drop = FALSE]
    
    train_inputs <- list(
      features     = array(x_train, c(nrow(x_train), input_dim)),
      class        = array(labels_train, c(nrow(labels_train), n_classes)),
      z_orig_batch = array(Z_original_train_proj, c(nrow(x_train), procrustes_dim))
    )
    
    val_inputs <- list(
      features     = array(x_val, c(nrow(x_val), input_dim)),
      class        = array(labels_val, c(nrow(labels_val), n_classes)),
      z_orig_batch = array(Z_original_val_proj, c(nrow(x_val), procrustes_dim))
    )
    
    batch_size <- min(batch_size, nrow(x_train))
    params      <- list(
      latent_dim              = latent_dim,
      intermediate_dim        = intermediate_dim,
      n_layers                = n_layers,
      learning_rate           = learning_rate,
      beta                    = beta,
      l2_reg                  = l2_reg,
      batch_size              = batch_size,
      epochs                  = epochs,
      init_mse_weight         = init_mse_weight,
      init_procrustes_weight  = init_procrustes_weight,
      final_procrustes_weight = final_procrustes_weight,
      phase_switch_epoch      = phase_switch_epoch
    )
    cvae_list <- def_cvae_conditional(
      params
      ,
      input_dim    = input_dim,
      n_classes    = n_classes,
      z_original_  = Z_original_train
    )
    cvae    <- cvae_list$model
    encoder <- cvae_list$encoder
    decoder <- cvae_list$decoder
    
    # Callbacks… (igual que antes)
    # dentro de build_and_evaluate_cvae_conditional(...)
    callbacks <- list(
      # Llama al constructor directamente
      BetaAnnealingCallback(
        epochs_anneal = min(40, epochs),
        beta_final    = beta
      ),
      # Igual para el ResetEpochCallback si lo definiste con new_callback_class()
      ResetEpochCallback(),
      # Callbacks estándar de Keras ya vienen como funciones
      callback_early_stopping(
        monitor               = "val_loss",
        patience              = 10,
        restore_best_weights  = TRUE
      ),
      callback_reduce_lr_on_plateau(
        monitor = "val_loss",
        factor  = 0.5,
        patience = 3,
        min_lr  = 1e-6
      )
    )
    
    
    # Ejemplo completo dentro de build_and_evaluate:
    history <- cvae %>% fit(
      x = train_inputs,
      y = x_train,
      validation_data = list(val_inputs, x_val),
      epochs    = epochs,
      batch_size = batch_size,
      callbacks = callbacks,
      verbose   = 1
    )
    custom_layer <- cvae$get_layer("procrustes_loss_layer")
    R_list[[i]] <- tensorflow::tf$keras$backend$get_value(cvae$get_layer("procrustes_loss_layer")$R)
    
    R_fold <- R_list[[i]]
    # 2. Define projection_layer y buildéalo
    
    # Reconstrucciones seguras
    preds_train <- safe_reconstruct(
      encoder,
      decoder,
      x_train,
      labels_train,
      Z_original_train,
      # <--- tu matriz de referencia recortada
      dim(x_train)
    )
    preds_val   <- safe_reconstruct(
      encoder,
      decoder,
      x_val,
      labels_val,
      Z_original_val,
      # <--- idem aquí
      dim(x_val)
    )
    
    # Extraemos el diccionario de Python como lista R
    py_hist   <- history$history
    loss_vec  <- unlist(py_hist[["loss"]])
    
    # 2) si 'val_loss' no existe, lo evaluamos directamente
    vloss_vec <- if ("val_loss" %in% names(py_hist)) {
      unlist(py_hist[["val_loss"]])
    } else {
      as.numeric(cvae$evaluate(val_inputs, x_val, verbose = 0)[[1]])
    }
    
    if (length(loss_vec) > 0 && length(vloss_vec) > 0) {
      metrics[[i]] <- list(
        train_loss = min(loss_vec),
        val_loss   = min(vloss_vec),
        train_mse  = mean((x_train - preds_train)^2),
        val_mse    = mean((x_val   - preds_val)^2),
        train_cor  = cor(as.vector(x_train), as.vector(preds_train)),
        val_cor    = cor(as.vector(x_val), as.vector(preds_val)),
        epochs     = length(loss_vec),
      )
    } else {
      warning(sprintf("Fold %d omitido: sin historia de pérdidas", i))
    }
    
  }
  
  if (length(metrics) == 0)
    stop("No valid folds (all partitions had empty classes).")
  
  avg_metrics <- list(
    train_loss = mean(map_dbl(metrics, "train_loss")),
    val_loss   = mean(map_dbl(metrics, "val_loss")),
    train_mse  = mean(map_dbl(metrics, "train_mse")),
    val_mse    = mean(map_dbl(metrics, "val_mse")),
    train_cor  = mean(map_dbl(metrics, "train_cor")),
    val_cor    = mean(map_dbl(metrics, "val_cor")),
    epochs     = mean(map_dbl(metrics, "epochs"))
  )
  
  return(list(
    params    = params,
    metrics   = avg_metrics,
    histories = map(metrics, "history"),
    R_matrices = R_list  # Nuevo: matrices de rotación por fold
  ))
}


reconstruct_omics <- function(Z_matrix, W_list) {
  synthetic_omics <- list()
  for (view_name in names(W_list)) {
    W <- W_list[[view_name]]
    if (ncol(Z_matrix) != ncol(W)) {
      stop(paste(
        "Dimensiones incompatibles para",
        view_name,
        "Z:",
        ncol(Z_matrix),
        "W:",
        ncol(W)
      ))
    }
    synthetic_omics[[view_name]] <- Z_matrix %*% t(W)
  }
  synthetic_omics
}
filtra_por_omicas <- function(W, omicas) {
  out <- list()
  for (nombre in names(W)) {
    out[[nombre]] <- W[[nombre]][intersect(rownames(W[[nombre]]), omicas[[nombre]]), drop = FALSE, ]
  }
  return(out)
}


extrae_por_omica <- function(datos, lista_clasificada) {
  lapply(lista_clasificada, function(biomarcadores) {
    # Para cada tipo de ómica, filtra columnas de cada dataframe en la lista original
    lapply(datos, function(df) {
      df[, intersect(biomarcadores, colnames(df)), drop = FALSE]
    })
  })
}



##====code=====

model <- readRDS("../../scripts_Data/modelo.rds")

##SO14 ## outlier del proteoma.

model@samples_metadata$sample

plot_factors(model, color_by = "grupo", factors = c(1, 2))
plot_factor(model, factors = c(1:2), color_by = "grupo")
plot_variance_explained(model)
# Obtener factores y pesos
factors <- get_expectations(model, "Z")
weights <- get_expectations(model, "W")
Z <- factors$group1
W_list <- weights
detach(package:MOFA2)



grupo <- model@samples_metadata$grupo

Z_scaled <- minmax_scale(Z)
# One-hot encoding
clase_levels <- levels(as.factor(model@samples_metadata$grupo))
n_classes <- length(clase_levels)
grupo_onehot <- to_categorical(as.integer(as.factor(grupo)) - 1, num_classes = n_classes)


##======CVAE=====


options(tensorflow.extract.disallow_out_of_bounds = FALSE)

# Annealing de beta global
current_beta <- 0
BetaAnnealingCallback <- keras3::new_callback_class(
  classname = "BetaAnnealingCallback",
  
  initialize = function(self, epochs_anneal, beta_final) {
    self$epochs_anneal <- epochs_anneal
    self$beta_final   <- beta_final
  },
  
  on_epoch_begin = function(self, epoch, logs = NULL) {
    # 1) calculamos el nuevo beta
    frac     <- min(1, epoch / self$epochs_anneal)
    new_beta <- frac * self$beta_final
    
    # 2) obtenemos la capa de pérdida
    layer <- self$model$get_layer("procrustes_loss_layer")
    
    # 3) recuperamos el backend de TF *dentro* del callback
    k <- tensorflow::tf$keras$backend
    
    # 4) actualizamos la variable de kl_weight
    k$set_value(layer$kl_weight_var, as.numeric(new_beta))
    
    cat(sprintf("\n-> Annealing beta: %0.5f (epoch %d)\n", new_beta, epoch))
  }
)



ResetEpochCallback <- keras::new_callback_class(
  classname = "ResetEpochCallback",
  on_epoch_begin = function(self, epoch, logs = NULL) {
    layer <- self$model$get_layer("procrustes_loss_layer")
    keras::k_set_value(layer$current_epoch, as.integer(0))
    cat(sprintf("\nReset Procrustes epoch → 0 at epoch %d\n", epoch))
  }
)


# Prueba con parámetros simples
## Parámetros de prueba (test_params) con nombres consistentes

tf$keras$backend$clear_session()
gc()
test_result <- def_cvae_conditional(
  params = list(
    latent_dim = 6,
    intermediate_dim = 64,
    n_layers = 2,
    learning_rate = 0.001,
    beta = 0.1,
    l2_reg = 1e-4,
    batch_size = 32,
    epochs = 10,
    init_mse_weight = 0.5,
    init_procrustes_weight = 0.5,
    final_procrustes_weight = 0.2,
    phase_switch_epoch = 5
  ),
  input_dim = 2,
  n_classes = 3,
  z_original_ = matrix(rnorm(60), ncol = 6)  # Ejemplo con 10 muestras
)

Z_original_matrix <- Z_scaled
# 1. Reducir la dimensionalidad del espacio de parámetros
# max_latent_dim <- min(6, ncol(Z_original_matrix)) # No mayor que el número de características

# ===== ParamSet para la optimización =====
param_space <- makeParamSet(
  makeIntegerParam("latent_dim", lower = 6, upper = 10),
  # Explora >6 dimensiones
  makeIntegerParam("intermediate_dim", lower = 64, upper = 256),
  makeIntegerParam("n_layers", lower = 3, upper = 5),
  makeNumericParam("learning_rate", lower = 1e-4, upper = 1e-2),
  # Rango más estrecho
  makeNumericParam("beta", lower = 0.01, upper = 0.3),
  makeNumericParam("l2_reg", lower = 0, upper = 1e-4),
  makeIntegerParam("batch_size", lower = 64, upper = 256),
  makeIntegerParam("epochs", lower = 150, upper = 500),
  makeNumericParam("init_mse_weight", lower = 0.5, upper = 1.5),
  makeNumericParam("init_procrustes_weight", lower = 0.1, upper = 1.0),
  makeNumericParam(
    "final_procrustes_weight",
    lower = 0.1,
    upper = 0.5
  ),
  # Reducido
  makeIntegerParam("phase_switch_epoch", lower = 30, upper = 200)
)
# ===== Función objetivo robusta =====
objective_function <- makeSingleObjectiveFunction(
  name    = "CVAE_Optimization_Procrustes",
  fn      = function(x) {
    # Convertir el vector x en lista de parámetros
    prm <- as.list(x)
    cat(">>> Evaluando: ",
        paste(
          names(prm),
          round(unlist(prm), 4),
          sep = "=",
          collapse = "; "
        ),
        "\n")
    
    # Ejecutar build_and_evaluate dentro de un try para capturar errores
    res <- tryCatch({
      # Validaciones básicas
      if (any(sapply(prm[c("latent_dim",
                           "init_mse_weight",
                           "init_procrustes_weight")], is.null)))
        stop("Param NULL detectado")
      if ((prm$init_mse_weight + (1 - prm$init_mse_weight)) > 1.1)
        stop("Pesos inválidos")
      # if (prm$latent_dim > ncol(Z_original_matrix))                   stop("latent_dim demasiado grande")
      
      # Invocar la función de evaluación de CVAE
      build_res <- build_and_evaluate_cvae_conditional(
        params       = list(
          latent_dim              = prm$latent_dim,
          intermediate_dim        = prm$intermediate_dim,
          n_layers                = prm$n_layers,
          learning_rate           = prm$learning_rate,
          beta                    = prm$beta,
          l2_reg                  = prm$l2_reg,
          batch_size              = prm$batch_size,
          epochs                  = prm$epochs,
          init_mse_weight         = prm$init_mse_weight,
          init_procrustes_weight  = 1 - prm$init_mse_weight,
          final_procrustes_weight = prm$final_procrustes_weight,
          phase_switch_epoch      = prm$phase_switch_epoch
        ),
        Z_data       = Z_original_matrix,
        grupo        = grupo,
        clase_levels = clase_levels,
        n_folds      = 3
      )
      
      val <- build_res$metrics$val_loss
      if (is.null(val) ||
          !is.finite(val))
        stop("val_loss inválido o NULL")
      cat("    -> val_loss real =", round(val, 4), "\n")
      val
    }, error = function(e) {
      # Imprimir detalles del error para debugging
      cat("    !!! ERROR (closure-debug):", e$message, "\n")
      # Mostrar llamada que falló (ayuda a localizar el closure)
      
      # Retornar alto valor penalizado
      1e6
    })
    return(res)
  },
  par.set  = param_space,
  minimize = TRUE,
  noisy    = TRUE
)


surrogate_model <- makeLearner(
  "regr.randomForest",
  predict.type = "se"   # si quieres poder usar EI
)


# ===== Configuración de control para MBO =====
ctrl <- makeMBOControl(propose.points = 1, y.name = "val_loss")
ctrl <- setMBOControlTermination(ctrl, iters = 15, time.budget = 1800)
ctrl <- setMBOControlInfill(
  ctrl,
  crit                    = makeMBOInfillCritEI(se.threshold = 0.5),
  opt                     = "focussearch",
  opt.focussearch.maxit   = 10,
  opt.focussearch.points  = 50
)
# ===== Diseño inicial de MBO =====
design <- generateDesign(n       = 12,
                         par.set = param_space,
                         fun     = lhs::maximinLHS)

# ===== Ejecución de la optimización =====
set.seed(123)
result <- mbo(
  fun     = objective_function,
  design  = design,
  learner = surrogate_model,
  control = ctrl,
  show.info = TRUE
)

# Guardar resultado intermedio
saveRDS(result, "mbo_intermediate_result.rds")
# Asegurarse de que el entrenami
# Asegurarse de que el entrenamiento final se ejecute correctamente
callback_terminate_on_naan <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (any(is.na(unlist(logs)))) {
      message("NaN detected in epoch ", epoch)
    }
  }
)

safe_reconstruct <- function(encoder,
                             decoder,
                             x,
                             labels,
                             z_orig_batch,
                             target_shape) {
  tryCatch({
    # 1) pasa x, labels y z_orig_batch al encoder
    enc_out <- predict(encoder, list(x, labels, z_orig_batch))
    z_mean  <- enc_out[[1]]
    # 2) decodifica a partir de z_mean y labels
    recon   <- predict(decoder, list(z_mean, labels))
    if (!all(dim(recon) == target_shape)) {
      stop(sprintf(
        "Shape mismatch: got %s, expected %s",
        paste(dim(recon), collapse = "x"),
        paste(target_shape, collapse = "x")
      ))
    }
    recon
  }, error = function(e) {
    warning("Reconstrucción falló: ", e$message)
    matrix(0, nrow = target_shape[1], ncol = target_shape[2])
  })
}

##==== FINAL TRAINING ====
if (!is.null(result)) {
  best_params <- result$x
  latent_dim  <- best_params$latent_dim
  
  # Hold-out 80/20
  set.seed(123)
  n         <- nrow(Z_scaled)
  train_idx <- sample(n, floor(0.8 * n))
  val_idx   <- setdiff(seq_len(n), train_idx)
  
  x_train    <- Z_scaled[train_idx, , drop = FALSE]
  x_val      <- Z_scaled[val_idx, , drop = FALSE]
  labels_train <- to_categorical(as.integer(factor(grupo[train_idx], levels =
                                                     clase_levels)) - 1, num_classes = length(clase_levels))
  labels_val <- to_categorical(as.integer(factor(grupo[val_idx], levels =
                                                   clase_levels)) - 1, num_classes = length(clase_levels))
  
  # Z_orig para Procrustes (primeras 6 dim latentes)
  z_orig_train <- Z_scaled[train_idx, 1:2, drop = FALSE]
  z_orig_val   <- Z_scaled[val_idx, 1:2, drop = FALSE]
  
  # Reconstruye CVAE con best_params
  final_cvae_list <- def_cvae_conditional(
    params      = best_params,
    input_dim   = ncol(Z_scaled),
    n_classes   = length(clase_levels),
    z_original_ = Z_scaled
  )
  final_cvae <- final_cvae_list$model
  encoder    <- final_cvae_list$encoder
  decoder    <- final_cvae_list$decoder
  
  # Prepara inputs
  train_inputs <- list(features     = x_train,
                       class        = labels_train,
                       z_orig_batch = z_orig_train)
  val_inputs <- list(features     = x_val,
                     class        = labels_val,
                     z_orig_batch = z_orig_val)
  
  # Entrena el CVAE final
  history <- final_cvae %>% fit(
    x               = train_inputs,
    y               = x_train,
    validation_data = list(val_inputs, x_val),
    epochs          = best_params$epochs,
    batch_size      = best_params$batch_size,
    callbacks       = list(
      callback_early_stopping(
        "val_loss",
        patience = 15,
        restore_best_weights = TRUE
      ),
      callback_reduce_lr_on_plateau(
        "val_loss",
        factor = 0.2,
        patience = 5,
        min_lr = 1e-6
      ),
      BetaAnnealingCallback(min(50, best_params$epochs), best_params$beta),
      callback_terminate_on_naan
    ),
    verbose = 2
  )
  
  # Reconstrucción “pura” para métricas
  preds_train <- safe_reconstruct(encoder,
                                  decoder,
                                  x_train,
                                  labels_train,
                                  z_orig_train,
                                  dim(x_train))
  preds_val   <- safe_reconstruct(encoder, decoder, x_val, labels_val, z_orig_val, dim(x_val))
  
  # Limpieza
  tf$keras$backend$clear_session()
  gc()
}  # —— Ahora, UNA SOLA VEZ: creamos el aligned_encoder ——
message("Creando aligned_encoder a partir del modelo final entrenado…")


##==== CREA ALIGNED_ENCODER (UNA VEZ) ====
message("Creando aligned_encoder…")
# Extrae R_final
custom_layer <- final_cvae$get_layer("procrustes_loss_layer")
R_final      <- tf$keras$backend$get_value(custom_layer$R)

# Capa de proyección (latent_dim → 6), inicial identidad
projection_layer <- layer_dense(
  units = 2,
  use_bias = FALSE,
  input_shape = c(latent_dim)
)
projection_layer(matrix(0, 1, latent_dim))
set_weights(projection_layer, list(diag(latent_dim)[, 1:2]))

# Función de rotación
rotation_layer <- function(x, name) {
  layer_lambda(
    f = function(z)
      tf$matmul(z, R_final),
    name = name
  )(x)
}

aligned_encoder <- keras_model(
  inputs = final_cvae_list$encoder$inputs,
  outputs = list(
    projection_layer(final_cvae_list$encoder$outputs[[1]]) %>% rotation_layer("aligned_z_mean"),
    final_cvae_list$encoder$outputs[[2]],
    # z_log_var
    projection_layer(final_cvae_list$encoder$outputs[[3]]) %>% rotation_layer("aligned_z_sample")
  )
)
message("✅ aligned_encoder listo.")
# —— Verificación y creación del encoder alineado ——

# proj_dim usado en tu modelo (coincide con inp_z_orig y CustomLossLayer$proj_dim)
proj_dim <- 2

generate_synthetic_samples_conditional_simple <- function(
    n_samples,
    clase_name,
    base_encoder,      # <--- usar el encoder original (NO el aligned)
    decoder,
    Z_scaled,
    grupo,
    clase_levels,
    n_classes,
    latent_dim,        # = best_params$latent_dim
    proj_dim = 2,      # <--- NUEVO
    noise_level = 0.05,
    adaptive_noise = TRUE
) {
  if (!clase_name %in% clase_levels) stop("Clase no encontrada: ", clase_name)
  idx <- which(grupo == clase_name)
  if (!length(idx)) stop("No hay muestras para la clase: ", clase_name)
  
  Z_class <- Z_scaled[idx, , drop = FALSE]
  
  # one-hot
  class_mat <- matrix(0, nrow = nrow(Z_class), ncol = n_classes)
  class_mat[, which(clase_levels == clase_name)] <- 1
  
  # z_orig_dummy con proj_dim columnas (NO latent_dim)
  z_orig_dummy <- matrix(0, nrow = nrow(Z_class), ncol = proj_dim)
  
  # Obtener z_mean en dimensiones latentes completas usando el encoder base
  enc_out   <- predict(base_encoder, list(Z_class, class_mat, z_orig_dummy))
  z_mean_rl <- enc_out[[1]]                       # shape: (n, latent_dim)
  
  # Estadística de mu y Sigma en latent_dim
  mu <- colMeans(z_mean_rl)
  Sigma <- if (nrow(z_mean_rl) < 3) diag(pmax(apply(z_mean_rl, 2, var), 1e-8)) else {
    s <- cov(z_mean_rl); diag(s)[diag(s) < 1e-8] <- 1e-8; as.matrix(Matrix::nearPD(s)$mat)
  }
  
  z_samples <- tryCatch(MASS::mvrnorm(n = n_samples, mu = mu, Sigma = Sigma),
                        error = function(e) matrix(mu, nrow = n_samples, byrow = TRUE) +
                          matrix(rnorm(n_samples * length(mu)), nrow = n_samples) * 0.1)
  
  # Decodificar: decoder espera (z_samples con latent_dim, class_gen)
  class_gen <- matrix(0, nrow = n_samples, ncol = n_classes)
  class_gen[, which(clase_levels == clase_name)] <- 1
  synth_scaled <- predict(decoder, list(z_samples, class_gen))
  
  # Ruido controlado
  noise_levels <- if (adaptive_noise) noise_level * sqrt(pmax(apply(Z_class, 2, var), 1e-12))
  else                noise_level * apply(Z_scaled, 2, sd)
  synth_scaled <- synth_scaled + matrix(rnorm(n_samples * ncol(synth_scaled)), nrow = n_samples) %*% diag(noise_levels)
  synth_scaled[!is.finite(synth_scaled)] <- 0
  synth_scaled
}


generate_balanced_synthetic_data <- function(
    n_per_class,
    base_encoder,     # <--- encoder original
    decoder,
    Z_scaled,
    grupo,
    clase_levels,
    latent_dim,
    proj_dim = 2,     # <--- NUEVO
    noise_level = 0.08,
    adaptive_noise = FALSE,
    align_with_real = TRUE,
    return_noise_info = FALSE
) {
  n_classes <- length(clase_levels)
  Z_synth_list <- lapply(clase_levels, function(cl) {
    tryCatch({
      generate_synthetic_samples_conditional_simple(
        n_samples       = n_per_class,
        clase_name      = cl,
        base_encoder    = base_encoder,      # <---
        decoder         = decoder,
        Z_scaled        = Z_scaled,
        grupo           = grupo,
        clase_levels    = clase_levels,
        n_classes       = n_classes,
        latent_dim      = latent_dim,
        proj_dim        = proj_dim,          # <---
        noise_level     = noise_level,
        adaptive_noise  = adaptive_noise
      )
    }, error = function(e) { warning("Error generando ", cl, ": ", e$message); NULL })
  })
  valid_idx <- !sapply(Z_synth_list, is.null)
  if (!any(valid_idx)) stop("No se pudo generar ninguna muestra sintética válida")
  Z_synth_final <- do.call(rbind, Z_synth_list[valid_idx])
  y_synth_final <- factor(rep(clase_levels[valid_idx], each = n_per_class)[seq_len(nrow(Z_synth_final))],
                          levels = clase_levels)
  
  # Escalado consistente
  real_mins   <- apply(Z_scaled, 2, min)
  real_ranges <- pmax(apply(Z_scaled, 2, function(x) diff(range(x))), 1)
  Z_scaled_final <- sweep(Z_synth_final, 2, real_mins, "-")
  Z_scaled_final <- sweep(Z_scaled_final, 2, real_ranges, "/")
  
  if (align_with_real) {
    centroids_real <- t(sapply(clase_levels, function(cl) if (sum(grupo == cl)) colMeans(Z_scaled[grupo == cl, , drop = FALSE]) else rep(NA, ncol(Z_scaled))))
    centroids_real <- centroids_real[complete.cases(centroids_real), , drop = FALSE]
    centroids_synth <- t(sapply(levels(y_synth_final), function(cl) colMeans(Z_scaled_final[y_synth_final == cl, , drop = FALSE])))
    if (nrow(centroids_real) >= 2 && nrow(centroids_synth) >= 2) {
      proc <- vegan::procrustes(centroids_real, centroids_synth, scale = TRUE)
      Z_scaled_final <- scale(Z_scaled_final, center = TRUE, scale = FALSE) %*% proc$rotation * proc$scale
      Z_scaled_final <- sweep(Z_scaled_final, 2, proc$translation, "+")
    } else warning("No se pudo aplicar Procrustes")
  }
  
  out <- list(features = Z_scaled_final, labels = y_synth_final,
              scaling_params = list(min = real_mins, range = real_ranges))
  if (return_noise_info) out$noise_info <- list(noise_level = noise_level, adaptive_noise = adaptive_noise)
  out
}



### Generación de muestras sintéticas usando el encoder alineado

n_bal <- 30
classes <- clase_levels  # usa los niveles definidos previamente

synthetic_data <- generate_balanced_synthetic_data(
  n_per_class   = 30,
  base_encoder  = final_cvae_list$encoder,   # <---
  decoder       = decoder,
  Z_scaled      = Z_scaled,
  grupo         = grupo,
  clase_levels  = clase_levels,
  latent_dim    = best_params$latent_dim,    # p.ej. 10
  proj_dim      = 2,                         # coincide con inp_z_orig
  noise_level   = 0.08,
  adaptive_noise = FALSE
)
# a tienes Z_scaled2 listo para downstream…
# Asignar nombres a las filas de synthetic_data$features
labels_ <- synthetic_data$labels

# Contador por clase
name_vec <- unlist(lapply(levels(labels_), function(cl) {
  idx <- which(labels_ == cl)
  n <- length(idx)
  paste0("synth_", cl, seq_len(n))
}))

# Reordena para que el largo coincida
rownames(synthetic_data$features) <- name_vec

# Ahora tus filas tendrán nombres como synth_FD1, synth_FD2, synth_NP1, etc.
head(synthetic_data$features)
colnames(synthetic_data$features) <- c("Factor1",
                                       "Factor2")



# ——— 1) Invierte el scaling de Z_scaled para los sintéticos ———

# Recupera mins y maxs guardados en el atributo de Z_scaled
mins   <- attr(Z_scaled, "min")
maxs   <- attr(Z_scaled, "max")
ranges <- maxs - mins

# synthetic_data$features está en [0,1] por minmax sobre Z_scaled;
# ahora lo llevamos a la escala original de Z:
Z_synth_orig <- sweep(synthetic_data$features, 2, ranges, FUN = "*")
Z_synth_orig <- sweep(Z_synth_orig, 2, mins,   FUN = "+")

library(expm)

# 1. Alineamiento GLOBAL de rotación y escala
mu_real_global <- colMeans(Z)
mu_synth_global <- colMeans(Z_synth_orig)
C_real_global <- cov(Z)
C_synth_global <- cov(Z_synth_orig)

# Transformación global
A_global <- sqrtm(C_real_global) %*% solve(sqrtm(C_synth_global))
Z_synth_global <- sweep(Z_synth_orig, 2, mu_synth_global, "-") %*% A_global

# 2. Ajuste por grupo (solo traslación)
Z_synth_hybrid <- matrix(NA, nrow = nrow(Z_synth_global), ncol = ncol(Z_synth_global))

for(g in levels(labels_)) {
  idx_real <- which(grupo == g)
  idx_synth <- which(labels_ == g)
  
  if (length(idx_real) > 0 && length(idx_synth) > 0) {
    mu_real_g <- colMeans(Z[idx_real, , drop = FALSE])
    mu_synth_g <- colMeans(Z_synth_global[idx_synth, , drop = FALSE])
    
    # Solo corrección de media (sin tocar covarianza)
    Z_synth_hybrid[idx_synth, ] <- sweep(
      Z_synth_global[idx_synth, , drop = FALSE],
      2, 
      mu_synth_g - mu_real_g, 
      "-"
    )
  }
}

# 3. Preparar datos para gráfico
# Preparar el data.frame con nombres correctos
# Preparar el data.frame con nombres correctos

colnames(Z_synth_hybrid) <- colnames(Z)
df_all <- rbind(
  data.frame(Z, Tipo = "Real", Grupo = grupo) %>% 
    `rownames<-`(rownames(Z)),
  data.frame(Z_synth_hybrid, Tipo = "Sintético", Grupo = labels_) %>% 
    `rownames<-`(rownames(Z_synth_orig))
)
colnames(df_all) <- gsub("X", "Factor", colnames(df_all))

# Asegurar que los nombres de los factores sean consistentes
df_all$Grupo <- factor(df_all$Grupo)

library(ggplot2)
library(scales)  # para scale_color_hue()

ggplot(df_all, aes(
  x      = Factor1,
  y      = Factor2,
  color  = Grupo,
  shape  = Tipo,
  weight = 1 / as.numeric(table(Grupo)[Grupo])
)) +
  geom_point(alpha = 0.7, size = 3, stroke = 1) +
  stat_ellipse(
    aes(linetype = Tipo),
    level       = 0.8,
    linewidth   = 0.8,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = c(16, 4)) +     # Círculo para Real, X para Sintético
  scale_color_hue() +                          # Deja que ggplot asigne automáticamente
  scale_linetype_manual(values = c("solid", "longdash")) +
  labs(
    title    = "Comparación de Factores: Datos Reales vs. Sintéticos",
    subtitle = "Elipses al 80% de confianza",
    x        = "Factor 1",
    y        = "Factor 2",
    color    = "Grupo",
    shape    = "Tipo de dato",
    linetype = "Tipo de dato"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.spacing   = unit(0.5, "cm"),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(nrow = 1, override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  )


###====diagnostico-=====
set.seed(123)
train_idx <- createDataPartition(df_all$Tipo, p = 0.8, list = FALSE)
train_df <- df_all[train_idx, ]
test_df <- df_all[-train_idx, ]

# Pesos inversamente proporcionales al desbalance
class_weights <- ifelse(train_df$Tipo == "Real",
                        nrow(train_df)/sum(train_df$Tipo == "Real"),
                        nrow(train_df)/sum(train_df$Tipo == "Sintético"))

# Entrenamiento con pesos
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

rf_model <- caret::train(
  Tipo ~ .,
  data = train_df,
  method = "rf",
  weights = class_weights,  # Aplicación de pesos
  metric = "ROC",
  trControl = ctrl
)

# Evaluación
test_pred <- predict(rf_model, test_df, type = "prob")
roc_obj <- pROC::roc(test_df$Tipo, test_pred$Sintético)

# Plot ROC
plot(roc_obj, print.auc = TRUE, main = "ROC con Pesos de Clase")




library(umap)

# Preparación de datos
set.seed(123)
umap_data <- umap(df_all[, grep("Factor", names(df_all))])

# Dataframe para ggplot
plot_data <- data.frame(
  UMAP1 = umap_data$layout[,1],
  UMAP2 = umap_data$layout[,2],
  Type = df_all$Tipo,
  grupo = df_all$Grupo
)

# Visualización
ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color =grupo, shape = Type)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.8) +
  labs(title = "Separación Real vs. Sintético en UMAP",
       subtitle = "Si los grupos se solapan, los sintéticos son realistas") +
  theme_minimal()


# Test de Kolmogorov-Smirnov
ks_test <- ks.test(
  df_all$Factor1[df_all$Tipo == "Real"],
  df_all$Factor1[df_all$Tipo == "Sintético"]
)

# Gráfico de densidad
ggplot(df_all, aes(x = Factor1, fill = Tipo)) +
  geom_density(alpha = 0.5) +
  labs(title = paste("Factor1: D =", round(ks_test$statistic, 3),
                     "p-value =", signif(ks_test$p.value, 3))) +
  theme_minimal()



write.csv(df_real,"./df_real.csv")



salvar <- list(pesos = weights,
               factores = df_all)

saveRDS(salvar,"./datos_para_modelar.rds")







