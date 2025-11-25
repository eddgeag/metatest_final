##====Carga librerias=====
library(readxl)
library(dplyr)
library(ggplot2)
library(zoo)
library(ggrepel)
# library(mixOmics)
# library(tidyr)
##=== Funciones====
calcular_mahalanobis_general <- function(datos,
                                         variables = NULL,
                                         umbral_chi2 = 0.975) {
  #' @param datos: Dataframe con las variables numéricas.
  #' @param variables: Vector con nombres de columnas a usar (si NULL, usa todas las numéricas).
  #' @param umbral_chi2: Percentil para el umbral de outlier (ej. 0.975 para 97.5%).
  
  # 1. Seleccionar variables
  if (is.null(variables)) {
    datos <- datos %>% select_if(is.numeric)  # Usar solo columnas numéricas
  } else {
    datos <- datos %>% select(all_of(variables))  # Usar variables especificadas
  }
  
  # 2. Verificar que hay datos válidos
  if (ncol(datos) == 0)
    stop("No hay variables numéricas para calcular Mahalanobis.")
  if (nrow(datos) < 2)
    stop("Se requieren al menos 2 observaciones.")
  
  # 3. Calcular media y matriz de covarianza (con manejo de NA/singularidades)
  media <- tryCatch(
    colMeans(datos, na.rm = TRUE),
    error = function(e)
      rep(NA, ncol(datos))
  )
  
  cov_matrix <- tryCatch(
    cov(datos, use = "complete.obs"),
    error = function(e) {
      warning("Matriz de covarianza singular. Se usará una matriz diagonal.")
      diag(ncol(datos)) * 1e6  # Matriz diagonal grande como fallback
    }
  )
  
  # 4. Calcular distancias
  dist_mahal <- tryCatch(
    mahalanobis(datos, center = media, cov = cov_matrix),
    error = function(e) {
      warning("Error al calcular distancias. Verifique los datos.")
      rep(NA, nrow(datos))
    }
  )
  
  # 5. Identificar outliers
  umbral <- qchisq(umbral_chi2, df = ncol(datos))
  outliers <- dist_mahal > umbral
  
  # 6. Resultado
  resultados <- data.frame(
    indice = rownames(datos),
    dist_mahal = dist_mahal,
    outlier = outliers,
    umbral_chi2 = umbral
  )
  
  return(resultados)
}
split_matrix <- function(mat, block = 6) {
  # índices de inicio de cada bloque
  starts <- seq(1, ncol(mat), by = block)
  # para cada inicio, extrae columnas i:(i+block-1)
  lapply(starts, function(i) {
    mat[, i:min(i + block - 1, ncol(mat)), drop = FALSE]
  })
}

fun_process_transcriptomics <- function(path) {
  raw <- suppressMessages(read_excel(path, sheet = 1, col_names = FALSE))
  
  # 1) Metadata de features
  feat_meta_names <- raw %>% slice(7) %>% dplyr::select(1:18) %>% unlist() %>% as.character()
  feat_meta <- raw %>%
    slice(8:n()) %>%
    dplyr::select(1:18) %>%
    as.data.frame(stringsAsFactors = FALSE)
  colnames(feat_meta) <- feat_meta_names
  
  # 2) Metadata de muestras
  sample_meta_names <- raw %>% slice(1:6) %>% pull(19) %>% as.character()
  sample_meta_vals <- raw %>%
    slice(1:6) %>%
    dplyr::select(20:60) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(sample_meta_vals) <- sample_meta_names
  sample_meta_vals <- t(sample_meta_vals)
  sample_meta_vals <- as.data.frame(sample_meta_vals, stringsAsFactors = FALSE)
  
  code_row <- which(tolower(sample_meta_names) %in% c("codigo final", "código final"))
  if (length(code_row) == 0) stop("No se encontró columna 'CÓDIGO FINAL' en metadata de muestras.")
  rownames(sample_meta_vals) <- sample_meta_vals[, code_row[1]]
  col_code <- colnames(sample_meta_vals)[code_row[1]]
  
  # 3) Matriz de expresión por feature
  data_mat <- raw %>% slice(8:n()) %>% select(20:60) %>% as.data.frame()
  colnames(data_mat) <- sample_meta_vals[[col_code]]
  rownames(data_mat) <- rownames(feat_meta)
  
  # Depuración según tu flujo original
  data_mat <- bind_cols(GeneSymbol = feat_meta$GeneSymbol, data_mat)
  auxlog <- complete.cases(data_mat)
  data_mat <- data_mat[auxlog, ]
  feat_meta <- feat_meta[auxlog, ]
  data_mat$GeneSymbol <- feat_meta$PrimaryAccession
  rownames(data_mat) <- feat_meta$FeatureNum
  expr_feat <- data_mat[, -1, drop = FALSE]
  
  grupos <- as.factor(gsub("[-_]\\d+$", "", sample_meta_vals[[col_code]]))
  sample_meta_vals$grupo <- grupos
  
  # 4) Anotación y best_id
  nz <- function(x) { x <- trimws(as.character(x)); x[x == ""] <- NA; x }
  anno <- feat_meta
  anno$EnsemblID        <- nz(anno$EnsemblID)
  anno$EntrezGeneID     <- nz(anno$EntrezGeneID)
  anno$RefSeqAccession  <- nz(anno$RefSeqAccession)
  anno$GeneSymbol       <- nz(anno$GeneSymbol)
  
  anno$best_id <- dplyr::coalesce(anno$EnsemblID,
                                  anno$EntrezGeneID,
                                  anno$RefSeqAccession,
                                  anno$GeneSymbol)
  anno$id_type <- dplyr::case_when(
    !is.na(anno$EnsemblID)       ~ "EnsemblID",
    !is.na(anno$EntrezGeneID)    ~ "EntrezGeneID",
    !is.na(anno$RefSeqAccession) ~ "RefSeqAccession",
    !is.na(anno$GeneSymbol)      ~ "GeneSymbol",
    TRUE                         ~ "NA"
  )
  anno <- anno[match(rownames(expr_feat), anno$FeatureNum), ]
  
  # Excluir controles si existieran
  if ("ControlType" %in% colnames(anno)) {
    keep <- is.na(anno$ControlType) | anno$ControlType == 0
    expr_feat <- expr_feat[keep, , drop = FALSE]
    anno <- anno[keep, , drop = FALSE]
  }
  
  coverage <- list(
    by_id_type = table(anno$id_type, useNA = "ifany"),
    prop_annotated = mean(!is.na(anno$best_id))
  )
  
  # 5) Colapso por gen/transcrito usando mayor IQR
  valid <- !is.na(anno$best_id)
  expr_valid <- expr_feat[valid, , drop = FALSE]
  anno_valid <- anno[valid, , drop = FALSE]
  
  if (nrow(expr_valid) > 0) {
    iqr_feat <- apply(expr_valid, 1, IQR, na.rm = TRUE)
    idx_list <- split(seq_len(nrow(expr_valid)), anno_valid$best_id)
    pick <- vapply(idx_list, function(ix) ix[which.max(iqr_feat[ix])], integer(1))
    
    expr_gene <- expr_valid[pick, , drop = FALSE]
    rn_best <- anno_valid$best_id[pick]
    rownames(expr_gene) <- rn_best
    
    # columnas legibles
    expr_gene <- cbind(
      BestID     = rn_best,
      GeneSymbol = anno_valid$GeneSymbol[pick],
      expr_gene
    )
    
    anno_gene <- anno_valid[pick, c("best_id","id_type","GeneSymbol","GeneName",
                                    "EnsemblID","EntrezGeneID","RefSeqAccession",
                                    "GenbankAccession","Cytoband"), drop = FALSE]
  } else {
    expr_gene <- expr_valid
    anno_gene <- anno_valid[, c("best_id","id_type","GeneSymbol","GeneName",
                                "EnsemblID","EntrezGeneID","RefSeqAccession",
                                "GenbankAccession","Cytoband"), drop = FALSE]
  }
  entrez_vec <- anno_gene$EntrezGeneID
  okE <- !is.na(entrez_vec) & entrez_vec != ""
  
  dg <- expr_gene[okE, -(1:2), drop = FALSE]  # expresión (sin BestID/GeneSymbol)
  ag <- anno_gene[okE, , drop = FALSE]        # anotación alineada
  
  iqr_rows <- apply(dg, 1, IQR, na.rm = TRUE)
  idx_list <- split(seq_len(nrow(dg)), entrez_vec[okE])
  pickE <- vapply(idx_list, function(ix) ix[which.max(iqr_rows[ix])], integer(1))
  
  data_by_entrez <- dg[pickE, , drop = FALSE]
  anno_by_entrez <- ag[pickE, , drop = FALSE]
  rownames(data_by_entrez) <- anno_by_entrez$EntrezGeneID
  rownames(anno_by_entrez) <- anno_by_entrez$EntrezGeneID
  # 6) Salida
  list(
    feature_metadata = feat_meta,
    samples_metadata = sample_meta_vals,
    data_by_feature  = expr_feat,
    anno_by_feature  = anno,
    data_by_gene     = expr_gene[,-c(1:2)],
    anno_by_gene     = anno_gene,
    # >>> NUEVO <<<
    data_by_entrez   = data_by_entrez,   # filas = EntrezID únicos
    anno_by_entrez   = anno_by_entrez,   # filas = EntrezID únicos
    coverage         = coverage,
    id_priority      = c("EnsemblID","EntrezGeneID","RefSeqAccession","GeneSymbol"),
    code_field       = col_code
  )
  
}


prepare_transcriptomics_for_model <- function(
    trans_obj, metadata, y,
    fixed_keep_ids = NULL, fixed_train_ids = NULL, fixed_test_ids = NULL,
    p_train = 0.7, seed = 123, min_test_per_class = 1,
    validate_with_bitr = FALSE, orgdb = "org.Hs.eg.db"   # <- opcional
){
  # --- 0) Construir matriz por Entrez (una fila por Entrez con mayor IQR) ---
  # Si el objeto ya trae data_by_entrez, úsalo; si no, constrúyelo desde data_by_gene/anno_by_gene.
  if (!is.null(trans_obj$data_by_entrez) && !is.null(trans_obj$anno_by_entrez)) {
    M_entrez <- as.matrix(trans_obj$data_by_entrez)  # filas=Entrez, cols=muestras
    anno_entrez <- trans_obj$anno_by_entrez
  } else {
    M0   <- as.matrix(trans_obj$data_by_gene)        # filas = best_id, cols = muestras
    anno <- trans_obj$anno_by_gene
    entrez_vec <- as.character(anno$EntrezGeneID)[match(rownames(M0), anno$best_id)]
    keep <- !is.na(entrez_vec) & nzchar(entrez_vec)
    M0   <- M0[keep, , drop = FALSE]
    entrez_vec <- entrez_vec[keep]
    iqr_rows  <- apply(M0, 1, IQR, na.rm = TRUE)
    idx_list  <- split(seq_len(nrow(M0)), entrez_vec)
    pick_rows <- vapply(idx_list, function(ix) ix[which.max(iqr_rows[ix])], integer(1))
    M_entrez  <- M0[pick_rows, , drop = FALSE]
    rownames(M_entrez) <- names(idx_list)            # rownames = Entrez únicos
    anno_entrez <- anno[match(rownames(M_entrez), anno$EntrezGeneID),
                        c("best_id","id_type","GeneSymbol","GeneName",
                          "EnsemblID","EntrezGeneID","RefSeqAccession",
                          "GenbankAccession","Cytoband"), drop = FALSE]
    rownames(anno_entrez) <- anno_entrez$EntrezGeneID
  }
  
  # --- 0b) Mapa de etiquetas Entrez -> GeneSymbol (rellena con Entrez si falta) ---
  id_to_label <- {
    gs <- anno_entrez$GeneSymbol
    names(gs) <- rownames(anno_entrez)              # nombres = Entrez
    gs[is.na(gs) | gs == ""] <- names(gs)[is.na(gs) | gs == ""]
    gs
  }
  
  # --- 0c) Validación/relleno con clusterProfiler::bitr (opcional) ---
  if (isTRUE(validate_with_bitr) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    # Resolver OrgDb: aceptar objeto o nombre de paquete
    OrgDb_obj <- NULL
    if (inherits(orgdb, "OrgDb")) {
      OrgDb_obj <- orgdb
    } else if (is.character(orgdb) && length(orgdb) == 1L) {
      # cargar el paquete si está disponible
      if (!requireNamespace(orgdb, quietly = TRUE)) {
        warning(sprintf("No se pudo cargar '%s'; se omite bitr.", orgdb))
      } else {
        # obtener el objeto OrgDb exportado con el mismo nombre
        OrgDb_obj <- try(getExportedValue(orgdb, orgdb), silent = TRUE)
        if (inherits(OrgDb_obj, "try-error") || !inherits(OrgDb_obj, "OrgDb")) {
          # fallback: buscar cualquier objeto OrgDb exportado
          nms <- getNamespaceExports(orgdb)
          cand <- Filter(function(x) inherits(getExportedValue(orgdb, x), "OrgDb"), nms)
          if (length(cand)) {
            OrgDb_obj <- getExportedValue(orgdb, cand[[1]])
          } else {
            warning(sprintf("No se encontró un objeto OrgDb válido en '%s'; se omite bitr.", orgdb))
            OrgDb_obj <- NULL
          }
        }
      }
    }
    if (!is.null(OrgDb_obj) && inherits(OrgDb_obj, "OrgDb")) {
      # ids a validar (ENTREZ únicamente)
      entrez_ids <- names(id_to_label)
      entrez_ids <- unique(entrez_ids[nzchar(entrez_ids)])
      if (length(entrez_ids)) {
        bit <- try(
          suppressMessages(
            clusterProfiler::bitr(entrez_ids,
                                  fromType = "ENTREZID",
                                  toType   = "SYMBOL",
                                  OrgDb    = OrgDb_obj)
          ),
          silent = TRUE
        )
        if (!inherits(bit, "try-error") && !is.null(bit) && nrow(bit) > 0) {
          map_sym <- setNames(as.character(bit$SYMBOL), as.character(bit$ENTREZID))
          # solo rellenar donde aún el "símbolo" es el propio entrez
          missing <- names(id_to_label)[ id_to_label == names(id_to_label) ]
          repl <- map_sym[missing]
          repl[is.na(repl) | !nzchar(repl)] <- missing
          id_to_label[missing] <- repl
        }
      }
    }
  }
  # --- 1) Pasar a muestras x genes (columnas = Entrez) ---
  Xg <- t(M_entrez)
  mode(Xg) <- "numeric"
  
  # --- 2) Alinear a metadata / y ---
  ids0 <- intersect(rownames(metadata), rownames(Xg))
  if (!length(ids0)) stop("No hay IDs comunes en transcriptómica.")
  metadata <- metadata[ids0, , drop = FALSE]
  Xg <- Xg[ids0, , drop = FALSE]
  if (is.null(names(y))) names(y) <- rownames(metadata)
  y <- droplevels(y[rownames(metadata)])
  if (anyNA(y)) stop("y contiene NA tras alinear.")
  
  # --- 3) Filtrar por IDs fijos (si llegan) ---
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(Xg))
    Xg <- Xg[ids_use, , drop = FALSE]
    metadata <- metadata[ids_use, , drop = FALSE]
    y <- droplevels(y[ids_use])
  }
  
  # --- 4) Split + escalado coherente ---
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- intersect(fixed_train_ids, rownames(Xg))
    te <- intersect(fixed_test_ids,  rownames(Xg))
    X_tr <- Xg[tr, , drop = FALSE]; X_te <- Xg[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]);       y_te <- droplevels(y[te])
    X_tr_sc <- scale(X_tr)
    ctr <- attr(X_tr_sc, "scaled:center"); scl <- attr(X_tr_sc, "scaled:scale")
    X_te_sc <- scale(X_te, center = ctr, scale = scl)
    ss <- list(mode="split_scaled_fixed",
               X_train=X_tr, X_test=X_te,
               X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
               tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
               y_train=y_tr, y_test=y_te, center=ctr, scale=scl)
  } else {
    ss <- split_scale(X = Xg, y = y, p = p_train, seed = seed,
                      remove_outliers = FALSE, do_split = TRUE,
                      min_test_per_class = min_test_per_class)
  }
  
  # --- 5) Helper TEST: exige mismas columnas (Entrez) en el mismo orden ---
  transform_test <- function(X_new){
    req_cols <- colnames(ss$X_train_scaled)
    if (!all(req_cols %in% colnames(X_new)))
      stop("X_new no contiene todos los Entrez del train.")
    X_new <- X_new[, req_cols, drop = FALSE]
    X_new <- as.matrix(X_new); mode(X_new) <- "numeric"
    X_new_scaled <- scale(X_new, center = ss$center, scale = ss$scale)
    list(X_new_scaled = X_new_scaled,
         tX_new_scaled = t(X_new_scaled))
  }
  
  list(
    X_all     = Xg,
    split     = ss,
    center    = ss$center,
    scale     = ss$scale,
    genes     = colnames(Xg),     # EntrezGeneID
    metadata  = metadata,
    id_to_label = id_to_label,    # <- mapa Entrez -> símbolo (validado si se pidió)
    transform_test = transform_test
  )
}

process_proteoma <- function(
    path_citoquina,
    path_miokina,
    sep_cito = "",
    sep_mio = "",
    metadata, y,
    fixed_keep_ids = NULL,
    fixed_train_ids = NULL,
    fixed_test_ids  = NULL,
    use_feats_for_outliers = NULL,
    remove_outliers = TRUE,
    p_train = 0.7,
    seed = 123,
    min_test_per_class = 1,
    epsilon = 0
){
  # 1) Carga
  citoquina <- read.delim(path_citoquina, sep = sep_cito, row.names = 1, check.names = FALSE)
  miokina   <- read.delim(path_miokina,   sep = sep_mio,  row.names = 1, check.names = FALSE)
  proteoma  <- cbind(citoquina, miokina)
  
  # 2) Alinear a metadata
  ids0 <- intersect(rownames(metadata), rownames(proteoma))
  if (!length(ids0)) stop("No hay IDs comunes.")
  metadata <- metadata[ids0, , drop = FALSE]
  proteoma <- proteoma[ids0, , drop = FALSE]
  if (is.null(names(y))) names(y) <- rownames(metadata)
  y <- droplevels(y[rownames(metadata)])
  
  # 3) Filtro por IDs fijos o outliers
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(proteoma))
    proteoma <- proteoma[ids_use, , drop = FALSE]
    metadata <- metadata[ids_use, , drop = FALSE]
    y        <- droplevels(y[ids_use])
  } else if (isTRUE(remove_outliers)) {
    proteoma <- fun_remove_outliers(proteoma, use_feats_for_outliers)
    ids1 <- intersect(rownames(metadata), rownames(proteoma))
    proteoma <- proteoma[ids1, , drop = FALSE]
    metadata <- metadata[ids1, , drop = FALSE]
    y        <- droplevels(y[ids1])
  }
  if (anyNA(y)) stop("y contiene NA tras filtrar/alinear.")
  
  # 4) Normalización composicional + log2
  rs <- rowSums(proteoma, na.rm = TRUE)
  if (any(rs == 0)) { if (epsilon <= 0) stop("Filas con suma 0."); rs <- rs + epsilon }
  X_log2 <- log2(sweep(proteoma, 1, rs, "/"))
  
  # 5) Split + escalado
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- intersect(fixed_train_ids, rownames(X_log2))
    te <- intersect(fixed_test_ids,  rownames(X_log2))
    X_tr <- X_log2[tr, , drop = FALSE]; X_te <- X_log2[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]);           y_te <- droplevels(y[te])
    X_tr_sc <- scale(X_tr)
    ctr <- attr(X_tr_sc, "scaled:center"); scl <- attr(X_tr_sc, "scaled:scale")
    X_te_sc <- scale(X_te, center = ctr, scale = scl)
    ss <- list(mode="split_scaled_fixed",
               X_train=X_tr, X_test=X_te,
               X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
               tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
               y_train=y_tr, y_test=y_te, center=ctr, scale=scl)
  } else {
    ss <- split_scale(X = X_log2, y = y, p = p_train, seed = seed,
                      remove_outliers = FALSE, do_split = TRUE,
                      min_test_per_class = min_test_per_class)
  }
  
  # 6) Helper para TEST futuro
  transform_test <- function(X_new){
    # mismas columnas y orden que el train
    req_cols <- colnames(ss$X_train_scaled)
    if (!all(req_cols %in% colnames(X_new)))
      stop("X_new no contiene todas las columnas del train.")
    X_new <- X_new[, req_cols, drop = FALSE]
    # misma normalización
    rsn <- rowSums(X_new, na.rm = TRUE)
    if (any(rsn == 0)) { if (epsilon <= 0) stop("X_new filas con suma 0."); rsn <- rsn + epsilon }
    X_new_log2 <- log2(sweep(X_new, 1, rsn, "/"))
    # mismo escalado
    X_new_scaled <- scale(X_new_log2, center = ss$center, scale = ss$scale)
    list(X_new_log2 = X_new_log2,
         X_new_scaled = X_new_scaled,
         tX_new_scaled = t(X_new_scaled))
  }
  
  list(X_log2_all = X_log2,
       split = ss,
       center = ss$center, scale = ss$scale,
       metadata = metadata,
       transform_test = transform_test)
}


split_scale <- function(
    X,                      # matriz/data.frame: n x p (obs x features)
    y = NULL,               # factor/clase para split estratificado (requerido si do_split=TRUE)
    p = 0.7,                # proporción train
    seed = 123,
    remove_outliers = FALSE,
    outlier_patterns = NULL,# vector de patrones regex a buscar en rownames(X)
    outlier_ids = NULL,     # vector exacto de IDs (rownames) a remover
    do_split = NULL,        # si NULL: se infiere (TRUE si remove_outliers, FALSE si no)
    min_test_per_class = 1
) {
  # Coerciones
  rnames <- rownames(X)
  clnames <- colnames(X)
  X <- as.matrix(as.data.frame(lapply(as.data.frame(X),FUN = as.numeric)))
  colnames(X) <- clnames
  rownames(X) <- rnames
  rn <- rownames(X)
  if (is.null(do_split)) do_split <- isTRUE(remove_outliers)
  
  # --- 1) Remoción de outliers (opcional) ---
  to_drop <- logical(nrow(X))
  if (isTRUE(remove_outliers)) {
    if (!is.null(outlier_patterns) && length(outlier_patterns) > 0) {
      if (is.null(rn)) stop("Se requieren rownames(X) para usar outlier_patterns.")
      rx <- paste(outlier_patterns, collapse = "|")
      to_drop <- to_drop | grepl(rx, rn)
    }
    if (!is.null(outlier_ids) && length(outlier_ids) > 0) {
      if (is.null(rn)) stop("Se requieren rownames(X) para usar outlier_ids.")
      to_drop <- to_drop | rn %in% outlier_ids
    }
  }
  removed_ids <- if (!is.null(rn)) rn[to_drop] else as.character(which(to_drop))
  Xf <- X[!to_drop, , drop = FALSE]
  yf <- if (!is.null(y)) y[!to_drop] else NULL
  
  # --- 2) Sin split: solo escalar todo con sus propios parámetros ---
  if (!isTRUE(do_split)) {
    X_scaled <- scale(Xf)
    ctr <- attr(X_scaled, "scaled:center")
    scl <- attr(X_scaled, "scaled:scale")
    tX_scaled <- t(X_scaled) # "traspuesta de X" ya escalada por columnas originales
    
    return(list(
      mode = "nosplit_scaled",
      removed_ids = removed_ids,
      X_scaled = X_scaled,
      tX_scaled = tX_scaled,
      center = ctr,
      scale = scl,
      n = nrow(Xf),
      p = ncol(Xf)
    ))
  }
  
  # --- 3) Con split estratificado ---
  if (is.null(yf)) stop("y es requerido para split estratificado.")
  if (anyNA(yf)) stop("y contiene NA.")
  yf <- as.factor(as.character(yf))

  set.seed(seed)
  idx_train <- integer(0)
  for (lev in levels(yf)) {
    idx <- which(yf == lev)
    n <- length(idx)
    if (n <= min_test_per_class) {
      stop(sprintf("Clase '%s' tiene %d obs (< min_test_per_class)", lev, n))
    }
    n_train <- max(1, min(n - min_test_per_class, round(p * n)))
    idx_train <- c(idx_train, sample(idx, n_train, replace = FALSE))
  }
  idx_train <- sort(idx_train)
  idx_test  <- setdiff(seq_len(nrow(Xf)), idx_train)
  
  X_train <- Xf[idx_train, , drop = FALSE]
  X_test  <- Xf[idx_test,  , drop = FALSE]
  y_train <- yf[idx_train]
  y_test  <- yf[idx_test]
  
  # --- 4) Escalado: train con sus propios parámetros, test con atributos de train ---
  X_train_sc <- scale(X_train)
  ctr <- attr(X_train_sc, "scaled:center")
  scl <- attr(X_train_sc, "scaled:scale")
  
  X_test_sc <- scale(X_test, center = ctr, scale = scl)
  
  # Traspuestas escaladas coherentes: usar las versiones ya escaladas y luego transponer
  tX_train_sc <- t(X_train_sc)
  tX_test_sc  <- t(X_test_sc)
  
  # --- 5) Salida ---
  list(
    mode = "split_scaled",
    removed_ids = removed_ids,
    idx_train = idx_train,
    idx_test = idx_test,
    X_train = X_train,           # sin escalar
    X_test  = X_test,            # sin escalar
    X_train_scaled = X_train_sc, # escalado con center/scale de train
    X_test_scaled  = X_test_sc,  # escalado con center/scale de train
    tX_train_scaled = tX_train_sc,
    tX_test_scaled  = tX_test_sc,
    y_train = y_train,
    y_test  = y_test,
    counts = list(
      total = table(yf),
      train = table(y_train),
      test  = table(y_test)
    ),
    center = ctr,
    scale  = scl,
    n_train = nrow(X_train),
    n_test  = nrow(X_test),
    p = ncol(Xf)
  )
}




# Congela IDs y split global para todos los bloques
define_split_ids <- function(
    samples_metadata,            # data.frame con rownames = IDs y columna de grupo
    grupo_col = "grupo",
    outlier_patterns = c("NP_02","NP_03","SO_14","6M"),
    outlier_ids = NULL,
    p = 0.6,
    seed = 42,
    min_test_per_class = 1
){
  if (is.null(rownames(samples_metadata)))
    stop("samples_metadata necesita rownames=IDs de muestra.")
  if (!grupo_col %in% colnames(samples_metadata))
    stop(sprintf("No existe la columna '%s' en samples_metadata.", grupo_col))
  
  ids_master <- rownames(samples_metadata)
  
  # y nombrado
  y_master <- setNames(as.factor(as.character(samples_metadata[[grupo_col]])), ids_master)
  if (anyNA(y_master)) stop("y_master contiene NA.")
  
  # outliers por patrón + lista explícita
  rm_pat <- if (length(outlier_patterns) > 0)
    ids_master[grepl(paste(outlier_patterns, collapse="|"), ids_master)] else character(0)
  rm_ids <- unique(c(rm_pat, intersect(ids_master, outlier_ids %||% character(0))))
  
  keep_ids <- setdiff(ids_master, rm_ids)
  if (length(keep_ids) == 0) stop("No quedan IDs tras remover outliers.")
  
  # split estratificado
  set.seed(seed)
  y_keep <- droplevels(y_master[keep_ids])
  
  idx_tr <- integer(0)
  for (lev in levels(y_keep)) {
    idx <- which(y_keep == lev)
    n <- length(idx)
    if (n <= min_test_per_class)
      stop(sprintf("Clase '%s' tiene %d obs (< min_test_per_class).", lev, n))
    n_tr <- max(1, min(n - min_test_per_class, round(p * n)))
    idx_tr <- c(idx_tr, sample(idx, n_tr, replace = FALSE))
  }
  idx_tr <- sort(idx_tr)
  idx_te <- setdiff(seq_along(keep_ids), idx_tr)
  
  train_ids <- keep_ids[idx_tr]
  test_ids  <- keep_ids[idx_te]
  
  list(
    ids_master = ids_master,
    y_master   = y_master,
    removed_ids = rm_ids,
    keep_ids   = keep_ids,
    train_ids  = train_ids,
    test_ids   = test_ids,
    p = p, seed = seed
  )
}

# helper para operador coalescencia
`%||%` <- function(a,b) if (!is.null(a)) a else b

process_metaboloma <- function(
    path_pos, path_neg,
    metadata, y,
    sep_pos = "", sep_neg = "",
    prefix_pos = "pos_", prefix_neg = "neg_",
    fixed_keep_ids = NULL,
    fixed_train_ids = NULL,
    fixed_test_ids  = NULL,
    align_by_suffix = FALSE,                 # TRUE si necesitas casar por sufijo numérico
    suffix_pattern  = ".*_(\\d+)",           # como en tu script
    p_train = 0.7, seed = 123, min_test_per_class = 1
){
  # --- leer POS ---
  df_raw <- readxl::read_excel(path_pos, col_names = FALSE)
  sample_names <- as.character(df_raw[1, -1])
  data_matrix <- as.data.frame(df_raw[-c(1, 2), ])
  rownames(data_matrix) <- data_matrix[[1]]
  data_matrix <- data_matrix[, -1, drop = FALSE]
  colnames(data_matrix) <- sample_names
  data_matrix[] <- lapply(data_matrix, as.numeric)
  df_final <- as.data.frame(t(data_matrix))
  colnames(df_final) <- paste0(prefix_pos, colnames(df_final))
  
  # --- leer NEG ---
  df_raw2 <- readxl::read_excel(path_neg, col_names = FALSE)
  sample_names2 <- as.character(df_raw2[1, -1])
  data_matrix2 <- as.data.frame(df_raw2[-c(1, 2), ])
  rownames(data_matrix2) <- data_matrix2[[1]]
  data_matrix2 <- data_matrix2[, -1, drop = FALSE]
  colnames(data_matrix2) <- sample_names2
  data_matrix2[] <- lapply(data_matrix2, as.numeric)
  df_final2 <- as.data.frame(t(data_matrix2))
  colnames(df_final2) <- paste0(prefix_neg, colnames(df_final2))
  
  # --- combinar ---
  metaboloma <- cbind(df_final, df_final2)
  if (anyNA(metaboloma)) {
    keep_cols <- which(colSums(is.na(metaboloma)) < nrow(metaboloma))
    metaboloma <- metaboloma[, keep_cols, drop = FALSE]
  }
  rownames(metaboloma) <- gsub("-", "_", rownames(df_final))
  
  # --- alinear a metadata / IDs fijos ---
  if (is.null(names(y))) names(y) <- rownames(metadata)
  if (!align_by_suffix) {
    ids0 <- intersect(rownames(metadata), rownames(metaboloma))
    if (!length(ids0)) stop("No hay IDs comunes entre metadata y metaboloma.")
    metadata <- metadata[ids0, , drop = FALSE]
    metaboloma <- metaboloma[ids0, , drop = FALSE]
    y <- droplevels(y[ids0])
  } else {
    ref_ids <- if (!is.null(fixed_keep_ids)) fixed_keep_ids else rownames(metadata)
    suf_ref <- as.numeric(gsub(suffix_pattern, "\\1", ref_ids))
    suf_met <- as.numeric(gsub(suffix_pattern, "\\1", rownames(metaboloma)))
    idx <- match(suf_ref, suf_met)
    if (anyNA(idx)) stop("No se pudieron casar sufijos entre metaboloma y referencia.")
    metaboloma <- metaboloma[idx, , drop = FALSE]
    rownames(metaboloma) <- ref_ids
    metadata <- metadata[ref_ids, , drop = FALSE]
    y <- droplevels(y[ref_ids])
  }
  if (anyNA(y)) stop("y contiene NA tras la alineación.")
  
  # --- filtrar por IDs fijos si se dan ---
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(metaboloma))
    metaboloma <- metaboloma[ids_use, , drop = FALSE]
    metadata   <- metadata[ids_use, , drop = FALSE]
    y          <- droplevels(y[ids_use])
  }
  
  # --- split + escalado coherente con otros bloques ---
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- intersect(fixed_train_ids, rownames(metaboloma))
    te <- intersect(fixed_test_ids,  rownames(metaboloma))
    X_tr <- metaboloma[tr, , drop = FALSE]; X_te <- metaboloma[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]);              y_te <- droplevels(y[te])
    X_tr_sc <- scale(X_tr)
    ctr <- attr(X_tr_sc, "scaled:center"); scl <- attr(X_tr_sc, "scaled:scale")
    X_te_sc <- scale(X_te, center = ctr, scale = scl)
    ss <- list(
      mode="split_scaled_fixed",
      X_train=X_tr, X_test=X_te,
      X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
      tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
      y_train=y_tr, y_test=y_te,
      center=ctr, scale=scl
    )
  } else {
    ss <- split_scale(
      X = metaboloma, y = y, p = p_train, seed = seed,
      remove_outliers = FALSE, do_split = TRUE,
      min_test_per_class = min_test_per_class
    )
  }
  
  # --- helper para TEST futuro ---
  transform_test <- function(X_new){
    if (!all(colnames(ss$X_train) %in% colnames(X_new)))
      stop("X_new no contiene todas las variables del train metabolómico.")
    X_new <- X_new[, colnames(ss$X_train), drop = FALSE]
    X_new <- as.matrix(X_new); mode(X_new) <- "numeric"
    X_new_scaled <- scale(X_new, center = ss$center, scale = ss$scale)
    list(
      X_new_scaled   = X_new_scaled,
      tX_new_scaled  = t(X_new_scaled)
    )
  }
  
  list(
    X_all = metaboloma,
    split = ss,
    center = ss$center,
    scale  = ss$scale,
    metadata = metadata,
    transform_test = transform_test
  )
}

process_clinicos <- function(
    path_excel,
    metadata, y,
    fixed_keep_ids = NULL,
    fixed_train_ids = NULL,
    fixed_test_ids  = NULL,
    p_train = 0.7, seed = 123, min_test_per_class = 1
){
  raw <- readxl::read_excel(path_excel, col_names = FALSE)
  
  # subset y nombres
  dat <- raw[4:32, 4:44]
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  nombres_sujetos <- as.character(raw[4:32, 3][[1]])
  col_names <- as.character(unlist(raw[3, 4:44]))
  rownames(dat) <- nombres_sujetos
  names(dat) <- col_names
  
  # tipos
  num_cols <- setdiff(names(dat), c("Grupo","Género","Edad"))
  dat[num_cols] <- lapply(dat[num_cols], as.numeric)
  dat <- dplyr::mutate(dat,
                       Grupo  = factor(Grupo),
                       Género = factor(Género),
                       Edad   = factor(Edad))
  
  # reemplazo filas 20:22
  idx_fix <- intersect(20:22, seq_len(nrow(dat)))
  if (length(idx_fix) > 0) {
    dat[idx_fix, c("Grupo","Género","Edad")] <-
      lapply(dat[idx_fix, c("Grupo","Género","Edad")], function(x)
        factor(replace(as.character(x), x == "sobrepeso", "obeso")))
  }
  
  # selección de variables
  vars <- c(
    "Grupo","Género","Edad","GGT","Glu","Col total","LDL","HDL","TG","PCR",
    "Insulina","HOMA-IR","QUICKI","TG-INDEX","FLI","TG/HDL","Peso","Altura","IMC",
    "Ratio cintura_0M","Ratio-cc","%BF-DEXA","CUN-BAE","VAT","VAT (%)","Masa_grasa_0M",
    "%VAT_grasa_Total_0M","kg Ginoide","% Ginoide","% Ginoide_grasa_Total",
    "kg Androide","% Androide","% Androide_grasa_Total","Magra (%)","Magra (kg)",
    "Ratio and/gin","FFM- 0M (masa magra total + masa ósea total)","Sis","Dias"
  )
  dat_final <- as.data.frame(dat[, vars, drop = FALSE])
  rownames(dat_final) <- nombres_sujetos
  
  # alineación con metadata/y
  ids0 <- intersect(rownames(metadata), rownames(dat_final))
  if (!length(ids0)) stop("No hay IDs comunes entre metadata y clínicos.")
  metadata <- metadata[ids0, , drop = FALSE]
  dat_final <- dat_final[ids0, , drop = FALSE]
  if (is.null(names(y))) names(y) <- rownames(metadata)
  y <- droplevels(y[rownames(metadata)])
  if (anyNA(y)) stop("y contiene NA tras alinear clínicos.")
  
  # numéricas
  num_cols2 <- setdiff(colnames(dat_final), c("Grupo","Género","Edad"))
  num_df <- as.data.frame(lapply(dat_final[, num_cols2, drop=FALSE], as.numeric))
  rownames(num_df) <- rownames(dat_final)
  
  # imputación de Altura si falta
  if ("Altura" %in% colnames(num_df)) {
    idx_bad <- !is.finite(num_df[["Altura"]]) | is.na(num_df[["Altura"]])
    if (any(idx_bad)) {
      med_alt <- stats::median(num_df[["Altura"]][is.finite(num_df[["Altura"]]) & num_df[["Altura"]]>0], na.rm = TRUE)
      if (!is.finite(med_alt) || med_alt <= 0) stop("No hay valores válidos para imputar Altura.")
      num_df[["Altura"]][idx_bad] <- med_alt
    }
  }
  
  # keep fijos
  if (!is.null(fixed_keep_ids)) {
    ids_use <- intersect(fixed_keep_ids, rownames(num_df))
    num_df   <- num_df[ids_use, , drop = FALSE]
    metadata <- metadata[ids_use, , drop = FALSE]
    y        <- droplevels(y[ids_use])
  }
  
  # log + escalado global
  if (any(num_df <= 0, na.rm = TRUE))
    stop("Hay valores <= 0 en clínicos numéricos; no se puede aplicar log() directo.")
  X_log <- log(num_df)
  X_scaled <- scale(X_log)
  ctr <- attr(X_scaled, "scaled:center"); scl <- attr(X_scaled, "scaled:scale")
  
  # split: respeta IDs fijos
  if (!is.null(fixed_train_ids) && !is.null(fixed_test_ids)) {
    tr <- fixed_train_ids[fixed_train_ids %in% rownames(X_scaled)]
    te <- fixed_test_ids [fixed_test_ids  %in% rownames(X_scaled)]
    if (!length(tr) || !length(te)) stop("IDs fijos no encontrados en clínicos.")
    
    X_tr_sc <- X_scaled[tr, , drop = FALSE]
    X_te_sc <- X_scaled[te, , drop = FALSE]
    y_tr <- droplevels(y[tr]); y_te <- droplevels(y[te])
    
    X_tr <- X_log[tr, , drop = FALSE]
    X_te <- X_log[te, , drop = FALSE]
    
    ss <- list(
      mode="split_scaled_fixed_global",
      X_train=X_tr, X_test=X_te,
      X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
      tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
      y_train=y_tr, y_test=y_te,
      center=ctr, scale=scl
    )
  } else {
    set.seed(seed)
    y_fac <- droplevels(y)
    idx_train <- integer(0)
    for (lev in levels(y_fac)) {
      idx <- which(y_fac == lev)
      n <- length(idx)
      if (n <= min_test_per_class)
        stop(sprintf("Clase '%s' tiene %d obs (< min_test_per_class)", lev, n))
      n_tr <- max(1, min(n - min_test_per_class, round(p_train * n)))
      idx_train <- c(idx_train, sample(idx, n_tr, replace = FALSE))
    }
    idx_train <- sort(idx_train)
    idx_test  <- setdiff(seq_len(nrow(X_scaled)), idx_train)
    
    X_tr_sc <- X_scaled[idx_train, , drop = FALSE]
    X_te_sc <- X_scaled[idx_test,  , drop = FALSE]
    y_tr <- y_fac[idx_train]; y_te <- y_fac[idx_test]
    X_tr <- X_log[idx_train, , drop = FALSE]
    X_te <- X_log[idx_test,  , drop = FALSE]
    
    ss <- list(
      mode="split_scaled_global",
      idx_train=idx_train, idx_test=idx_test,
      X_train=X_tr, X_test=X_te,
      X_train_scaled=X_tr_sc, X_test_scaled=X_te_sc,
      tX_train_scaled=t(X_tr_sc), tX_test_scaled=t(X_te_sc),
      y_train=y_tr, y_test=y_te,
      center=ctr, scale=scl
    )
  }
  
  transform_test <- function(X_new){
    X_new <- X_new[, colnames(X_scaled), drop = FALSE]
    if (any(X_new <= 0, na.rm=TRUE)) stop("X_new contiene <=0")
    X_new_log <- log(as.matrix(X_new))
    X_new_scaled <- scale(X_new_log, center=ctr, scale=scl)
    list(X_new_log=X_new_log,
         X_new_scaled=X_new_scaled,
         tX_new_scaled=t(X_new_scaled))
  }
  
  list(
    X_log_all = X_scaled,
    split     = ss,
    center    = ctr, scale = scl,
    metadata  = metadata,
    cols_numeric = colnames(X_scaled),
    transform_test = transform_test
  )
}

save_block_plots <- function(
    res, out_dir, name,
    group_col = "grupo",
    label_col = "CÓDIGO.FINAL",
    palette = c(FD="red", SO="blue", NP="black"),
    width = 7, height = 5, dpi = 300, quality = 95,
    print_plots = TRUE
){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if ("X_log2_all" %in% names(res)) X0 <- res$X_log2_all
  else if ("X_log_all" %in% names(res)) X0 <- res$X_log_all
  else if ("X_all" %in% names(res)) X0 <- res$X_all
  else stop("El objeto no contiene X_log2_all, X_log_all ni X_all.")
  
  X0 <- as.matrix(X0)
  meta <- res$metadata
  stopifnot(all(rownames(X0) %in% rownames(meta)))
  meta_all <- meta[rownames(X0), , drop = FALSE]
  if (!group_col %in% colnames(meta_all)) stop("group_col no existe en metadata.")
  
  ctr <- res$center; scl <- res$scale
  if (is.null(ctr) || is.null(scl)) stop("Faltan center/scale en res.")
  Xz_all <- scale(X0[, names(ctr), drop = FALSE], center = ctr, scale = scl)
  
  tr_ids <- rownames(res$split$X_train)
  X0_tr  <- X0[tr_ids, , drop = FALSE]
  Xz_tr  <- scale(X0_tr[, names(ctr), drop = FALSE], center = ctr, scale = scl)
  meta_tr <- meta[tr_ids, , drop = FALSE]
  
  # => filename primero, luego subset/proc
  .save <- function(p, fn, subset_type, proc_type){
    dir_target <- file.path(out_dir, name, subset_type, proc_type)
    if (!dir.exists(dir_target)) dir.create(dir_target, recursive = TRUE)
    ggsave(
      filename = fn, path = dir_target, plot = p,
      width = width, height = height, units = "in",
      device = "jpeg", dpi = dpi, quality = quality
    )
    if (isTRUE(print_plots)) print(p)
  }
  
  .dens <- function(X, meta_df, ttl){
    df <- cbind(meta_df[, group_col, drop = FALSE], as.data.frame(X))
    colnames(df)[1] <- "grupo"
    mlt <- reshape2::melt(df, id.vars = "grupo")
    ggplot(mlt, aes(value, color = grupo)) +
      geom_density(aes(y = after_stat(density)), alpha = 0.25, linewidth = 0.7) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
      labs(title = ttl, x = "Valor", y = "Densidad", color = "Grupo") +
      scale_color_manual(values = palette)
  }
  
  .pca <- function(X, meta_df, ttl){
    pc <- prcomp(X, center = FALSE, scale. = FALSE)
    var <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 2)
    lbl <- if (label_col %in% colnames(meta_df)) meta_df[[label_col]] else rownames(meta_df)
    plt <- cbind(as.data.frame(pc$x[, 1:2, drop=FALSE]),
                 grupo = as.factor(meta_df[[group_col]]),
                 label = lbl)
    ggplot(plt, aes(PC1, PC2, color = grupo, label = label)) +
      geom_point(size = 2, alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
      labs(title = ttl,
           x = paste0("PC1: ", var[1], "%"),
           y = paste0("PC2: ", var[2], "%"),
           color = "Grupo") +
      scale_color_manual(values = palette)
  }
  
  .plsda <- function(X, meta_df, ttl){
    if (!requireNamespace("mixOmics", quietly = TRUE)) stop("mixOmics no instalado.")
    y <- as.factor(meta_df[[group_col]])
    mdl <- mixOmics::plsda(X, y, ncomp = 2, scale = FALSE)
    var <- round(100 * mdl$prop_expl_var$X, 2)
    sc <- as.data.frame(mdl$variates$X)
    sc$grupo <- y
    sc$label <- if (label_col %in% colnames(meta_df)) meta_df[[label_col]] else rownames(meta_df)
    ggplot(sc, aes(comp1, comp2, color = grupo, label = label)) +
      geom_point(size = 2, alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
      labs(title = ttl,
           x = paste0("PC1: ", var[1], "%"),
           y = paste0("PC2: ", var[2], "%"),
           color = "Grupo") +
      scale_color_manual(values = palette)
  }
  
  # ALL
  .save(.dens(X0,     meta_all, paste0(name, " · Densidad · ALL · sin procesar")),  "density_all.jpg",  "ALL",  "raw")
  .save(.dens(Xz_all, meta_all, paste0(name, " · Densidad · ALL · procesado")),     "density_all.jpg",  "ALL",  "processed")
  .save(.pca(X0,      meta_all, paste0(name, " · PCA · ALL · sin procesar")),       "pca_all.jpg",      "ALL",  "raw")
  .save(.pca(Xz_all,  meta_all, paste0(name, " · PCA · ALL · procesado")),          "pca_all.jpg",      "ALL",  "processed")
  .save(.plsda(X0,    meta_all, paste0(name, " · PLS-DA · ALL · sin procesar")),    "plsda_all.jpg",    "ALL",  "raw")
  .save(.plsda(Xz_all,meta_all, paste0(name, " · PLS-DA · ALL · procesado")),       "plsda_all.jpg",    "ALL",  "processed")
  
  # TRAIN
  .save(.dens(X0_tr,   meta_tr, paste0(name, " · Densidad · TRAIN · sin procesar")), "density_train.jpg","TRAIN","raw")
  .save(.dens(Xz_tr,   meta_tr, paste0(name, " · Densidad · TRAIN · procesado")),    "density_train.jpg","TRAIN","processed")
  .save(.pca(X0_tr,    meta_tr, paste0(name, " · PCA · TRAIN · sin procesar")),      "pca_train.jpg",    "TRAIN","raw")
  .save(.pca(Xz_tr,    meta_tr, paste0(name, " · PCA · TRAIN · procesado")),         "pca_train.jpg",    "TRAIN","processed")
  .save(.plsda(X0_tr,  meta_tr, paste0(name, " · PLS-DA · TRAIN · sin procesar")),   "plsda_train.jpg",  "TRAIN","raw")
  .save(.plsda(Xz_tr,  meta_tr, paste0(name, " · PLS-DA · TRAIN · procesado")),      "plsda_train.jpg",  "TRAIN","processed")
}


# --- 1) PLS-DA SOLO CON TRAIN: loadings + histogramas (comp1 y comp2) ---
run_plsda_importance_train <- function(
    res, name, out_dir = "plots",
    id_to_label = NULL,
    palette = c(FD="red", SO="blue", NP="black")
){
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("mixOmics no instalado.")
  Xtr <- res$split$X_train_scaled
  ytr <- res$split$y_train
  if (is.null(Xtr) || is.null(ytr)) stop("Falta X_train_scaled o y_train.")
  
  mdl <- mixOmics::plsda(Xtr, ytr, ncomp = 2, scale = FALSE)
  
  L_raw <- data.frame(
    variable = rownames(mdl$loadings$X),
    comp1 = as.numeric(mdl$loadings$X[,1]),
    comp2 = as.numeric(mdl$loadings$X[,2]),
    stringsAsFactors = FALSE
  )
  
  get_group_from_plotLoadings <- function(m, comp) {
    aux <- mixOmics::plotLoadings(m, comp = comp, contrib = "max", plot = FALSE)
    mat <- NULL
    if (!is.null(aux$mat) && is.matrix(aux$mat)) mat <- aux$mat
    else if (!is.null(aux$X) && is.data.frame(aux$X)) {
      bool_cols <- vapply(aux$X, function(z) is.logical(z) || all(z %in% c(0,1,TRUE,FALSE), na.rm=TRUE), TRUE)
      if (any(bool_cols)) mat <- as.matrix(aux$X[, bool_cols, drop=FALSE])
    }
    if (is.null(mat)) stop("plotLoadings no devolvió matriz lógica de clases.")
    if (is.null(rownames(mat))) rownames(mat) <- rownames(m$loadings$X)
    grp <- colnames(mat)[max.col(mat, ties.method = "first")]
    data.frame(variable = rownames(mat), comp = comp, grupo_raw = grp, row.names = NULL)
  }
  G <- rbind(get_group_from_plotLoadings(mdl, 1), get_group_from_plotLoadings(mdl, 2))
  
  # limpiar nombres de grupo: "Contrib.FD" -> "FD"
  clean_group <- function(x){
    x <- as.character(x)
    x <- gsub("^(?i)(contrib\\.|class\\.|group\\.)", "", x, perl = TRUE)
    toupper(trimws(x))
  }
  G$grupo <- clean_group(G$grupo_raw)
  
  win <- L_raw |>
    dplyr::mutate(abs1 = abs(comp1), abs2 = abs(comp2)) |>
    dplyr::transmute(variable,
                     comp_win = ifelse(abs1 >= abs2, 1L, 2L),
                     loading_win = ifelse(abs1 >= abs2, comp1, comp2),
                     abs_loading_win = pmax(abs1, abs2)) |>
    dplyr::left_join(G[, c("variable","comp","grupo")], by = c("variable","comp_win" = "comp"))
  
  # Etiquetas legibles
  if (!is.null(id_to_label)) {
    lab <- id_to_label[win$variable]
    lab[is.na(lab) | lab == ""] <- win$variable
    win$label_var <- lab
  } else {
    win$label_var <- win$variable
  }
  
  # Importancia (para la siguiente fase)
  importance_train <- setNames(win$abs_loading_win, win$variable)
  
  # Colores por grupo (FD/SO/NP)
  grupos_presentes <- unique(win$grupo)
  pal_use <- palette[names(palette) %in% grupos_presentes]
  # si aparece algún grupo “raro”, color gris
  extra <- setdiff(grupos_presentes, names(palette))
  if (length(extra) > 0) pal_use <- c(pal_use, setNames(rep("grey70", length(extra)), extra))
  
  # Top-40
  top40 <- win |>
    dplyr::arrange(dplyr::desc(abs_loading_win)) |>
    dplyr::slice(1:40) |>
    dplyr::arrange(loading_win)
  
  p_top40 <- ggplot(top40,
                    aes(x = loading_win,
                        y = factor(label_var, levels = top40$label_var),
                        fill = factor(grupo, levels = names(pal_use)))) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = pal_use, breaks = names(pal_use), labels = names(pal_use), drop = FALSE) +
    theme_minimal(base_size = 11) +
    labs(title = paste0("Top-40 loadings (con signo) · TRAIN · ", name),
         x = "Loading", y = NULL, fill = "Grupo") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  
  dir_base <- file.path(out_dir, name, "feature_prefilter_train")
  if (!dir.exists(dir_base)) dir.create(dir_base, recursive = TRUE)
  f_top40 <- file.path(dir_base, "top40_importancia_piramidal_signo.png")
  ggsave(f_top40, p_top40, width = 8, height = 8, dpi = 300)
  print(p_top40)
  
  invisible(list(
    model_train = mdl,
    importance_train = importance_train,
    top40_table = top40[, c("variable","label_var","loading_win","abs_loading_win","grupo","comp_win")],
    paths = list(top40 = f_top40)
  ))
}



# --- helper para obtener la matriz BASE "normalizada sin escalar" por ómica ---
.get_base_matrix_unscaled <- function(res) {
  if ("X_log2_all" %in% names(res)) {
    return(as.matrix(res$X_log2_all))                 # proteoma (CLR+log2)
  } else if ("X_log_all" %in% names(res)) {
    # clínicos: X_log_all en res_cli es scale(log(datos)) -> desescalar a log crudo
    if (is.null(res$center) || is.null(res$scale))
      stop("Faltan center/scale para desescalar X_log_all.")
    X_scaled <- as.matrix(res$X_log_all)
    X_log <- sweep(sweep(X_scaled, 2, res$scale, `*`), 2, res$center, `+`)
    return(X_log)                                      # log(datos) sin escalar
  } else if ("X_all" %in% names(res)) {
    return(as.matrix(res$X_all))                       # metaboloma: intensidades ya limpias
  } else {
    stop("No encuentro X_log2_all / X_log_all / X_all en 'res'.")
  }
}

# --- 2) UMBRAL SOLO CON TRAIN, re-normaliza por ómica, aplica a TEST, refit y plots ---
# === Extensión: devolver también muestras×vars y un "paquete final" coherente ===
apply_importance_threshold_and_refit_from_train <- function(
    res, importance_train, q = 0.9, name, out_dir = "plots") {
  
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("mixOmics no instalado.")
  if (is.null(names(importance_train)))
    stop("importance_train debe estar nombrado por variable (nombres = variables).")
  
  # Umbral y selección SOLO usando TRAIN
  thr <- as.numeric(stats::quantile(importance_train, q))
  sel_vars <- names(importance_train)[importance_train >= thr]
  if (!length(sel_vars)) stop("Umbral demasiado alto: 0 variables seleccionadas.")
  
  # Matriz base sin escalar (según ómica)
  X_base <- .get_base_matrix_unscaled(res)            # obs x vars
  keep <- intersect(colnames(X_base), sel_vars)
  if (!length(keep)) stop("Las variables seleccionadas no están en la matriz base.")
  X_base_sel <- X_base[, keep, drop = FALSE]
  
  # IDs
  tr_ids <- rownames(res$split$X_train_scaled)
  te_ids <- rownames(res$split$X_test_scaled)
  if (is.null(tr_ids) || is.null(te_ids))
    stop("Faltan IDs de train/test en res.")
  
  # Re-escalado NUEVO con SOLO TRAIN
  X_tr_base <- X_base_sel[tr_ids, , drop = FALSE]
  X_tr_sc   <- scale(X_tr_base)
  ctr_new   <- attr(X_tr_sc, "scaled:center")
  scl_new   <- attr(X_tr_sc, "scaled:scale")
  X_all_sc  <- scale(X_base_sel, center = ctr_new, scale = scl_new)
  X_te_sc   <- X_all_sc[te_ids, , drop = FALSE]
  
  # Refit PLS-DA SOLO en TRAIN filtrado
  y_tr <- res$split$y_train
  mdl2 <- mixOmics::plsda(X_tr_sc, y_tr, ncomp = 2, scale = FALSE)
  var_pc <- round(100 * mdl2$prop_expl_var$X, 2)
  
  # ---- Salidas en ambas orientaciones ----
  # 1) vars×muestras (como ya usabas para gráficos/diagnóstico)
  X_vs_train <- t(X_tr_sc)
  X_vs_test  <- t(X_te_sc)
  
  # 2) muestras×vars (para detectar par y alimentar modelos)
  X_sv_train <- X_tr_sc
  X_sv_test  <- X_te_sc
  stopifnot(identical(colnames(X_sv_train), colnames(X_sv_test)))
  stopifnot(length(intersect(rownames(X_sv_train), rownames(X_sv_test))) == 0)
  
  # ---- Plots de umbral y scores (opcional, igual que antes) ----
  dir_base <- file.path(out_dir, name, "feature_filtering_from_train")
  if (!dir.exists(dir_base)) dir.create(dir_base, recursive = TRUE)
  
  df_imp <- data.frame(variable = names(importance_train),
                       abs_loading = as.numeric(importance_train))
  p_hist_thr <- ggplot(df_imp, aes(abs_loading)) +
    geom_histogram(bins = 30, color = "black") +
    geom_vline(xintercept = thr, linetype = 2) +
    theme_minimal(11) +
    labs(title = paste0("Umbral ", q, " (", name, ", TRAIN)"),
         x = "|loading| (max comp1–2 en TRAIN)", y = "Frecuencia")
  ggsave(file.path(dir_base, sprintf("hist_importancia_train_q%.2f.jpg", q)),
         p_hist_thr, width = 7, height = 5, dpi = 300)
  print(p_hist_thr)
  
  scores_tr <- as.data.frame(mdl2$variates$X)[, 1:2, drop = FALSE]
  colnames(scores_tr) <- c("comp1","comp2")
  df_tr <- cbind(scores_tr, grupo = as.factor(y_tr), id = rownames(scores_tr))
  p_scores_tr <- ggplot(df_tr, aes(comp1, comp2, color = grupo, label = id)) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    theme_minimal(11) +
    labs(title = paste0("PLS-DA final (", name, ") · TRAIN"),
         x = paste0("PC1: ", var_pc[1], "%"),
         y = paste0("PC2: ", var_pc[2], "%"),
         color = "Grupo")+ scale_color_manual(values = c(FD="red", SO="blue", NP="black"))

  ggsave(file.path(dir_base, "plsda_scores_train.jpg"), p_scores_tr, width = 7, height = 5, dpi = 300)
  print(p_scores_tr)
  
  invisible(list(
    selected_vars = keep,
    threshold = thr,
    center = ctr_new, scale = scl_new,
    # -> orientación vars×muestras
    X_vs_train = X_vs_train,
    X_vs_test  = X_vs_test,
    # -> orientación muestras×vars  **NUEVO**
    X_sv_train = X_sv_train,
    X_sv_test  = X_sv_test,
    model_train_final = mdl2,
    train_ids = tr_ids,
    test_ids  = te_ids,
    paths = list(
      hist_threshold_train = file.path(dir_base, sprintf("hist_importancia_train_q%.2f.jpg", q)),
      scores_train         = file.path(dir_base, "plsda_scores_train.jpg")
    )
  ))
}

# === Extractor estándar para modelado: devuelve siempre muestras×vars alineado ===
extract_for_model <- function(fx_obj) {
  if (!all(c("X_sv_train","X_sv_test") %in% names(fx_obj))) {
    stop("No se hallaron matrices finales en fx_obj. Falta X_sv_train/X_sv_test.")
  }
  Xtr <- fx_obj$X_sv_train
  Xte <- fx_obj$X_sv_test
  if (!identical(colnames(Xtr), colnames(Xte)))
    stop("Columnas de train/test no coinciden.")
  list(
    X_train = Xtr,
    X_test  = Xte,
    tX_train = t(Xtr),
    tX_test  = t(Xte),
    center = fx_obj$center,
    scale  = fx_obj$scale,
    vars   = colnames(Xtr),
    train_ids = rownames(Xtr),
    test_ids  = rownames(Xte)
  )
}



##==codigo====

# 1) Metadata base
transcriptomica <- fun_process_transcriptomics("../data/transcriptomica.xlsx")

# 2) Split global único
sp <- define_split_ids(
  samples_metadata = transcriptomica$samples_metadata,
  grupo_col = "grupo",
  outlier_patterns = c("NP_02","NP_03","SO_14","6M"),
  p = 0.6, seed = 42
)

# 3) Proteoma (NO volver a quitar outliers aquí)
res_prot <- process_proteoma(
  path_citoquina = "../data/citokina.csv",
  path_miokina   = "../data/miokina.csv",
  metadata = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  remove_outliers = FALSE
)

# 4) Transcriptómica
res_tx <- prepare_transcriptomics_for_model(
  trans_obj = transcriptomica,
  metadata  = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  validate_with_bitr = TRUE,
  orgdb = "org.Hs.eg.db"
)

# 5) Metaboloma
res_met <- process_metaboloma(
  path_pos = "../data/5_03. METAHEALTH_POS_3 grupos_filtered_normalized.xlsx",
  path_neg = "../data/5_04. METAHEALTH_NEG_3 grupos_filtered_normalized.xlsx",
  metadata = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids,
  align_by_suffix = TRUE,
  suffix_pattern  = ".*_(\\d+)"
)

# 6) Clínicos
res_cli <- process_clinicos(
  path_excel = "../data/clinicos.xlsx",
  metadata   = transcriptomica$samples_metadata,
  y = sp$y_master,
  fixed_keep_ids  = sp$keep_ids,
  fixed_train_ids = sp$train_ids,
  fixed_test_ids  = sp$test_ids
)

stopifnot(
  identical(rownames(res_prot$split$X_train_scaled), rownames(res_tx$split$X_train_scaled)),
  identical(rownames(res_prot$split$X_train_scaled), rownames(res_met$split$X_train_scaled)),
  identical(rownames(res_prot$split$X_train_scaled), rownames(res_cli$split$X_train_scaled)),
  identical(rownames(res_prot$split$X_test_scaled),  rownames(res_tx$split$X_test_scaled)),
  identical(rownames(res_prot$split$X_test_scaled),  rownames(res_met$split$X_test_scaled)),
  identical(rownames(res_prot$split$X_test_scaled),  rownames(res_cli$split$X_test_scaled))
)

# === Nuevos TEST (ejemplos) ===
# Proteoma:
# new_prot <- cbind(citoquina_new, miokina_new); rownames(new_prot) <- new_ids
# Xt_prot <- res_prot$transform_test(new_prot)$X_new_scaled
# Transcriptómica:
# new_tx <- t(new_data_by_gene_new)  # muestras x genes (Entrez como columnas)
# Xt_tx  <- res_tx$transform_test(new_tx)$X_new_scaled

# === Plots: usar shim para clínicos (desescalar X_log_all antes de pasarlo) ===
save_block_plots(res_tx,   out_dir = "plots/Transcriptoma", name = "Transcriptoma")
save_block_plots(res_prot, out_dir = "plots/Proteoma",      name = "Proteoma")
save_block_plots(res_met,  out_dir = "plots/Metaboloma",     name = "Metaboloma")

res_cli_plot <- res_cli
res_cli_plot$X_log_all <- sweep(sweep(res_cli$X_log_all, 2, res_cli$scale, `*`), 2, res_cli$center, `+`)
save_block_plots(res_cli_plot, out_dir = "plots/Clinicos",  name = "Clinicos")

sapply(res_cli$split[c("X_train","X_test","X_train_scaled","X_test_scaled")], function(m){
  M <- as.matrix(m)
  c(nNA = sum(is.na(M)), nInf = sum(is.infinite(M)))
})

# === Importancia SOLO con TRAIN por ómica ===
# Mapeo Entrez -> símbolo para transcriptómica
id_to_label_tx <- {
  a <- transcriptomica$anno_by_entrez
  setNames(ifelse(is.na(a$GeneSymbol) | a$GeneSymbol == "", rownames(a), a$GeneSymbol),
           rownames(a))
}

tx_imp <- run_plsda_importance_train(res_tx,  name = "Transcriptoma", id_to_label = id_to_label_tx)
pr_imp <- run_plsda_importance_train(res_prot, name = "Proteoma")
me_imp <- run_plsda_importance_train(res_met,  name = "Metaboloma")
cl_imp <- run_plsda_importance_train(res_cli,  name = "Clinicos")

imp_tx <- tx_imp$importance_train
imp_pr <- pr_imp$importance_train
imp_me <- me_imp$importance_train
imp_cl <- cl_imp$importance_train

# === Umbral por cuantil en TRAIN, refit y plots ===
fx_tx <- apply_importance_threshold_and_refit_from_train(res_tx,  imp_tx, q = 0.97, name = "Transcriptoma")
fx_pr <- apply_importance_threshold_and_refit_from_train(res_prot, imp_pr, q = 0.2, name = "Proteoma")
fx_me <- apply_importance_threshold_and_refit_from_train(res_met,  imp_me, q = 0.97, name = "Metaboloma")
fx_cl <- apply_importance_threshold_and_refit_from_train(res_cli,  imp_cl, q = 0.4, name = "Clinicos")




# ========= VERIFICACIÓN + EMPAQUETADO (sin re-procesar ni re-filtrar) =========

# --- Helpers robustos ---
is_dfm <- function(x) inherits(x, c("data.frame","matrix"))

# Localiza dentro de un objeto los DF/matrices de train/test sin asumir nombres exactos
.guess_pair <- function(obj){
  flat <- unlist(obj, recursive = FALSE, use.names = TRUE)
  # Mantener solo data.frame/matrix
  flat <- flat[vapply(flat, is_dfm, logical(1))]
  nms  <- names(flat)
  # Heurísticas de nombres comunes
  train_keys <- grep("(^|_)X?_*train($|_)|(^|\\.)train($|\\.)", nms, ignore.case = TRUE, value = TRUE)
  test_keys  <- grep("(^|_)X?_*test($|_)|(^|\\.)test($|\\.)",  nms, ignore.case = TRUE, value = TRUE)
  # Si hay varias, prioriza las escaladas/seleccionadas
  ord <- function(v){
    score <- 0L +
      grepl("scaled", v, ignore.case=TRUE) * 4L +
      grepl("sel|select", v, ignore.case=TRUE) * 8L +
      grepl("X_", v) * 1L
    -score
  }
  tk <- if (length(train_keys)) train_keys[order(ord(train_keys))][1] else NA_character_
  te <- if (length(test_keys))  test_keys[order(ord(test_keys))][1]   else NA_character_
  list(train_key=tk, test_key=te, train=if (!is.na(tk)) flat[[tk]] else NULL, test=if (!is.na(te)) flat[[te]] else NULL)
}

# Dado el objeto "fx_*" después de apply_importance_threshold..., recupera la versión final para modelado

.get_fx_pair <- function(fx){
  # 1) prioridad: nuevas salidas del apply_* que ya dejaste en muestras×vars
  if (is_dfm(fx$X_sv_train) && is_dfm(fx$X_sv_test)) {
    return(list(train = fx$X_sv_train, test = fx$X_sv_test))
  }
  # 2) fallback: si solo tienes vars×muestras, transpón
  if (is_dfm(fx$X_vs_train) && is_dfm(fx$X_vs_test)) {
    return(list(train = t(fx$X_vs_train), test = t(fx$X_vs_test)))
  }
  # 3) otros nombres comunes en split
  cand <- list(
    train = list(
      fx$split$X_train_selected_scaled, fx$split$X_train_sel_scaled,
      fx$split$X_train_scaled, fx$split$X_train_sel, fx$split$X_train
    ),
    test = list(
      fx$split$X_test_selected_scaled, fx$split$X_test_sel_scaled,
      fx$split$X_test_scaled, fx$split$X_test_sel, fx$split$X_test
    )
  )
  pick <- function(lst){
    for (m in lst) if (is_dfm(m) && !is.null(rownames(m)) && !is.null(colnames(m))) return(m)
    NULL
  }
  list(train = pick(cand$train), test = pick(cand$test))
}

# Construye resumen para una ómica
verify_block <- function(imp_obj, importance_train, q, fx_obj, label){
  # 1) localizar pares train/test usados en la etapa de importancia (sin re-procesar)
  g  <- .guess_pair(imp_obj)
  Xt <- g$train; Xv <- g$test
  
  if (is.null(Xt) || is.null(Xv))
    return(list(
      label=label, ok=FALSE, note="No se localizaron matrices _train/_test dentro del objeto de importancia.",
      details=NULL
    ))
  
  # 2) dimensiones
  dims <- list(
    train = c(nrow=nrow(Xt), ncol=ncol(Xt)),
    test  = c(nrow=nrow(Xv), ncol=ncol(Xv))
  )
  
  # 3) IDs de filas
  rn_tr <- rownames(Xt); rn_te <- rownames(Xv)
  inter_ids <- intersect(rn_tr, rn_te)
  test_subset_of_all <- all(rn_te %in% c(rn_tr, rn_te)) # trivialmente TRUE; nos interesa disyunción con train
  disjoint <- length(inter_ids) == 0
  
  # 4) Columnas iguales entre train/test
  same_cols <- identical(colnames(Xt), colnames(Xv))
  
  # 5) Correspondencia con filtro por importancia (sin re-filtrar datos, solo cálculo del conjunto esperado)
  thr <- as.numeric(stats::quantile(importance_train, q))
  vars_selected_by_q <- names(importance_train)[importance_train >= thr]
  # En la función real se usa intersect con la base disponible; verificamos igualdad con lo efectivamente presente
  cols_tr <- colnames(Xt)
  cols_te <- colnames(Xv)
  expected_keep <- vars_selected_by_q[vars_selected_by_q %in% union(cols_tr, cols_te)]
  match_filter_train <- setequal(cols_tr, expected_keep)
  match_filter_test  <- setequal(cols_te, expected_keep)
  
  # 6) Matrices finales post-apply_* para empaquetado
  fxp <- .get_fx_pair(fx_obj)
  
  list(
    label=label, ok=all(disjoint, same_cols, match_filter_train, match_filter_test, !is.null(fxp$train), !is.null(fxp$test)),
    note=if (!disjoint) sprintf("Advertencia: %d IDs en común entre train y test.", length(inter_ids)) else "OK",
    details=list(
      dims=dims,
      rows=list(
        n_train=dims$train["nrow"], n_test=dims$test["nrow"],
        intersect_ids=inter_ids, disjoint=disjoint
      ),
      cols=list(
        same_train_test=same_cols,
        n_vars_train=dims$train["ncol"], n_vars_test=dims$test["ncol"],
        match_filter_train=match_filter_train,
        match_filter_test=match_filter_test,
        n_selected_by_q=length(vars_selected_by_q)
      ),
      fx_available = list(train=!is.null(fxp$train), test=!is.null(fxp$test))
    )
  )
}

# --- Ejecutar verificación para las cuatro ómicas ---
summ_tx <- verify_block(tx_imp, imp_tx, q=0.97, fx_tx, label="Transcriptómica")
summ_pr <- verify_block(pr_imp, imp_pr, q=0.2, fx_pr, label="Proteoma")
summ_me <- verify_block(me_imp, imp_me, q=0.97, fx_me, label="Metaboloma")
summ_cl <- verify_block(cl_imp, imp_cl, q=0.40, fx_cl, label="Clínicos")

verification_summary <- list(tx=summ_tx, pr=summ_pr, me=summ_me, cl=summ_cl)

# --- Impresión compacta del summary ---
print_summary <- function(s){
  cat(sprintf("\n[%s]\n", s$label))
  if (!isTRUE(s$ok)) cat("Estado: CHECK con observaciones\n") else cat("Estado: OK\n")
  d <- s$details
  cat(sprintf("Dims: train=%dx%d | test=%dx%d\n",
              d$dims$train["nrow"], d$dims$train["ncol"],
              d$dims$test["nrow"],  d$dims$test["ncol"]))
  cat(sprintf("Filas: disjoint=%s | intersección=%d\n",
              d$rows$disjoint, length(d$rows$intersect_ids)))
  cat(sprintf("Columnas: iguales train/test=%s | match filtro (train)=%s | match filtro (test)=%s | n_selected_by_q=%d\n",
              d$cols$same_train_test, d$cols$match_filter_train, d$cols$match_filter_test,
              d$cols$n_selected_by_q))
  cat(sprintf("fx_obj disponibles: train=%s, test=%s\n",
              d$fx_available$train, d$fx_available$test))
  if (!isTRUE(s$ok)) cat(sprintf("Nota: %s\n", s$note))
}

cat("\n================= SUMMARY DE VERIFICACIÓN =================\n")
print_summary(summ_tx); print_summary(summ_pr); print_summary(summ_me); print_summary(summ_cl)
cat("\n===========================================================\n")

# --- Empaquetado final: ready_for_modeling ---
# Extrae matrices finales del objeto fx_* y las formatea como: variables en filas, sujetos en columnas
extract_for_model <- function(fx_obj){
  fxp <- .get_fx_pair(fx_obj)
  if (is.null(fxp$train) || is.null(fxp$test)) stop("No se hallaron matrices finales en fx_obj.")
  # Verifica columnas idénticas y filas disjuntas
  stopifnot(identical(colnames(fxp$train), colnames(fxp$test)))
  stopifnot(length(intersect(rownames(fxp$train), rownames(fxp$test))) == 0)
  # Formato final requerido: variables en filas, sujetos en columnas
  list(
    train = t(as.matrix(fxp$train)),
    test  = t(as.matrix(fxp$test)),
    vars  = colnames(fxp$train),
    ids_train = rownames(fxp$train),
    ids_test  = rownames(fxp$test)
  )
}

# === Extrae matrices finales (muestras×vars) de cada ómica ===
ex_tx <- extract_for_model(fx_tx)
ex_pr <- extract_for_model(fx_pr)
ex_me <- extract_for_model(fx_me)
ex_cl <- extract_for_model(fx_cl)

make_grupo_df <- function(obj) {
  extract_pref <- function(ids) sub("_.*", "", ids)
  list(
    train = data.frame(id = obj$ids_train,
                       grupo = factor(extract_pref(obj$ids_train),
                                      levels = c("NP","FD","SO"))),
    test  = data.frame(id = obj$ids_test,
                       grupo = factor(extract_pref(obj$ids_test),
                                      levels = c("NP","FD","SO")))
  )
}

grupo_df <- make_grupo_df(ex_cl)

# grupo_df$train
# grupo_df$test


lista <- list(tx = ex_tx,
              pr = ex_pr,
              me = ex_me,
              cl = ex_cl,
              grupo_train =grupo_df$train,
              grupo_test = grupo_df$test,
              features_metadata = transcriptomica$feature_metadata)

saveRDS(lista,"./ready_for_modeling.rds")


# Ajusta estos si quieres cambiar estética/etiquetas
group_col <- "grupo"
label_col <- "ID"  # no existe en meta_df; la función usará los rownames
palette   <- c(FD="red", SO="blue", NP="black")

# vector de grupos con nombres = ID de sujeto
grupo_vec <- setNames(
  as.factor(transcriptomica$samples_metadata$grupo),
  rownames(transcriptomica$samples_metadata)
)

plsda_plot <- function(X, meta_df, ttl,
                       group_col = "grupo",
                       label_col = NULL,
                       palette = c(FD="red", SO="blue", NP="black")) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) stop("mixOmics no instalado.")
  
  y <- as.factor(meta_df[[group_col]])
  mdl <- mixOmics::plsda(X, y, ncomp = 2, scale = FALSE)
  var <- round(100 * mdl$prop_expl_var$X, 2)
  
  sc <- as.data.frame(mdl$variates$X)
  sc$grupo <- y
  if (!is.null(label_col) && label_col %in% colnames(meta_df)) {
    sc$label <- meta_df[[label_col]]
  } else {
    sc$label <- rownames(meta_df)
  }
  
  ggplot(sc, aes(comp1, comp2, color = grupo, label = label)) +
    geom_point(size = 2, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")) +
    labs(title = ttl,
         x = paste0("PC1: ", var[1], "%"),
         y = paste0("PC2: ", var[2], "%"),
         color = "Grupo") +
    scale_color_manual(values = palette)
}

plot_plsda_split <- function(lst, split = c("train","test")) {
  split <- match.arg(split)
  for (omic in c("tx","pr","me","cl")) {
    M <- lst[[omic]][[split]]            # vars × sujetos
    ids <- colnames(M)                   # sujetos
    meta_df <- data.frame(grupo = grupo_vec[ids], row.names = ids)
    X <- t(M)                            # muestras × vars
    ttl <- sprintf("PLS-DA · %s · %s", toupper(omic), toupper(split))
    print(plsda_plot(X, meta_df, ttl))   # <-- usa la versión global
  }
}

plot_plsda_split(lista, "train")
plot_plsda_split(lista, "test")
