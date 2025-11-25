# =======================================================
# Paquetes
# =======================================================
suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(forcats)
  library(psych)        # corr.test
  library(corpcor)      # cor.shrink
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
})

# =======================================================
# Carga de modelos y selección
# =======================================================
modelos <- list.files("./modelos", full.names = TRUE)
modelo  <- lapply(modelos, load_model)
plot_models <- MOFA2::compare_elbo(modelo)  # gráfico ELBO
print(plot_models)
rfm <- readRDS("./ready_for_modeling.rds")

modelo_final <- select_model(modelo)        # asume función disponible (MOFA2 >= 1.4)

# Diagnósticos estándar del usuario
plot_factor_cor(modelo_final)
plot_factors(modelo_final, color_by = "grupo")
plot_factor(modelo_final, color_by = "grupo", factors = c(1,2))
plot_variance_explained(modelo_final)

# =======================================================
# Helpers genéricos
# =======================================================
get_Z <- function(m) {
  Z <- get_factors(m, factors = "all", as.data.frame = FALSE)
  if (is.list(Z)) Z <- Z[[1]]
  as.matrix(Z)
}

get_W_list <- function(m) {
  Ws <- get_weights(m, views = "all", factors = "all", as.data.frame = FALSE)
  lapply(Ws, function(W) as.matrix(W))  # features x factors
}

get_A_list <- function(m) {
  # Interceptos por vista (si existen)
  out <- try(get_intercepts(m, as.data.frame = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

get_metadata <- function(m) {
  md <- try(get_sample_metadata(m), silent = TRUE)
  if (inherits(md, "try-error")) return(NULL)
  md
}

safe_group_col <- function(md, col = "grupo") {
  if (is.null(md) || !(col %in% names(md))) return(factor(rep("NA", nrow(md))))
  factor(md[[col]])
}

# Reconstrucción por vista: X_hat = Z %*% t(W) + intercepto
reconstruct_views <- function(m) {
  Z  <- get_Z(m)                                # n x k
  Ws <- get_W_list(m)                           # lista de (p_v x k)
  As <- get_A_list(m)                           # lista de interceptos por vista (p_v)
  lapply(names(Ws), function(v) {
    Wv <- Ws[[v]]                               # p x k
    X  <- Z %*% t(Wv)                           # n x p
    if (!is.null(As) && v %in% names(As)) {
      a <- As[[v]]
      if (length(a) == ncol(X)) X <- sweep(X, 2, a, "+")
    }
    colnames(X) <- rownames(Wv)                 # features en columnas
    X
  }) |> setNames(names(Ws))
}

# PCA genérico
do_pca <- function(X, scale. = TRUE) {
  # quita columnas con var=0 o NA
  keep <- apply(X, 2, function(z) is.finite(sd(z, na.rm=TRUE)) && sd(z, na.rm=TRUE) > 0)
  X2 <- X[, keep, drop=FALSE]
  prcomp(X2, center = TRUE, scale. = scale.)
}

plot_scores <- function(pca, meta=NULL, color_col="grupo", title="Score plot") {
  sc <- as.data.frame(pca$x)[, 1:2, drop=FALSE]
  sc$sample <- rownames(sc)
  if (!is.null(meta)) {
    grp <- safe_group_col(meta, color_col)
    sc[[color_col]] <- grp
  } else {
    sc[[color_col]] <- factor("NA")
  }
  ggplot(sc, aes(PC1, PC2, color=.data[[color_col]])) +
    geom_point(size=2, alpha=0.9) +
    labs(title=title, x=paste0("PC1 (", percent(summary(pca)$importance[2,1]), ")"),
         y=paste0("PC2 (", percent(summary(pca)$importance[2,2]), ")"),
         color=color_col) +
    theme_minimal()
}

plot_loadings <- function(pca, n_top=20, title="Loading plot") {
  ld <- as.data.frame(pca$rotation)[, 1:2, drop=FALSE]
  ld$feature <- rownames(ld)
  ld$absPC1  <- abs(ld$PC1)
  top <- ld |> arrange(desc(absPC1)) |> slice_head(n=n_top)
  ggplot(top, aes(reorder(feature, absPC1), PC1)) +
    geom_col() +
    coord_flip() +
    labs(title=title, x=NULL, y="Loading PC1") +
    theme_minimal()
}

plot_scree <- function(pca, title="Scree plot") {
  var_exp <- summary(pca)$importance[2,]
  df <- data.frame(PC = seq_along(var_exp), VarExp = var_exp)
  ggplot(df, aes(PC, VarExp)) +
    geom_col() +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels=scales::percent) +
    labs(title=title, x="PC", y="% Var explicada") +
    theme_minimal()
}

# =======================================================
# Análisis de pesos y selección por umbral
# =======================================================
W_long <- get_weights(modelo_final, views="all", factors="all", as.data.frame=TRUE) %>%
  rename(weight=value) %>%
  mutate(abs_w = abs(weight))
W_long <- W_long %>% 
  mutate(feature = if_else(
    view == "transcriptomica",
    coalesce(
      rfm$features_metadata$GeneSymbol[match(feature, rfm$features_metadata$EntrezGeneID)],
      feature
    ),
    feature  # Valor para cuando view != "transcriptomica"
  ))
# Umbral por cuantil dentro de cada (vista,factor)
thr_q <- 0.97
thr_tbl <- W_long %>%
  group_by(view, factor) %>%
  summarise(thr = quantile(abs_w, thr_q, na.rm=TRUE), .groups="drop")



W_thr <- W_long %>%
  left_join(thr_tbl, by=c("view","factor")) %>%
  mutate(selected = abs_w >= thr)

# Gráfico de |peso| con umbral por factor
plot_abs_weights_with_thr <- function(W_thr_df, view_sel, factor_sel, n_max=2000) {
  df <- W_thr_df %>% filter(view==view_sel, factor==factor_sel) %>%
    arrange(desc(abs_w)) %>%
    slice_head(n=n_max) %>%
    mutate(feature = fct_reorder(feature, abs_w))
   p <- ggplot(df, aes(feature, abs_w, fill=selected)) +
    geom_col() +
    geom_hline(aes(yintercept=thr), color="black", linetype="dashed") +
    coord_flip() +
    labs(title=paste0("||peso|| con umbral (", view_sel, ", ", factor_sel, ")"),
         x=NULL, y="|weight|") +
    theme_minimal() +
    guides(fill="none")
   print(p)
   return(W_thr_df$feature[which(W_thr_df$selected == T & W_thr_df$view == view_sel & W_thr_df$factor ==factor_sel)])
}

# Factores más importantes por vista según var explicada
var_exp <- MOFA2::plot_variance_explained(modelo_final,plot_total =F)



transform_to_r2_matrix <- function(X) {
  # Transformar a formato wide: filas = view, columnas = factor
  r2_per_factor <- X %>%
    select(factor, view, value) %>%
    pivot_wider(
      names_from = factor,
      values_from = value
    )
  
  # Convertir a matriz y establecer los nombres de fila
  result_matrix <- as.matrix(r2_per_factor[, -1])
  rownames(result_matrix) <- r2_per_factor$view
  
  # Agregar fila de r2_total (suma por columna/factor)
  r2_total <- colSums(result_matrix, na.rm = TRUE)
  result_matrix <- rbind(result_matrix, r2_total = r2_total)
  
  return(result_matrix)
}

# Usar la función
var_exp_r2_matrix <- transform_to_r2_matrix(var_exp$data)
print(var_exp_r2_matrix)

# var_exp$r2_per_factor es matriz: filas = vistas (+ "r2_total"), cols = factores
r2_pf <- var_exp_r2_matrix
views  <- setdiff(rownames(r2_pf), "r2_total")
top_m  <- 2  # nº de factores top por vista
top_factors_by_view <- lapply(views, function(v) {
  ord <- order(r2_pf[v,], decreasing = TRUE)
  which(colnames(r2_pf) %in% colnames(r2_pf)[ord[seq_len(min(top_m, length(ord)))]])
}) |> setNames(views)

# Conjuntos de features:
# - union_sel: features que superan umbral en cualquier factor de su vista
# - inter_sel: intersección solo de los factores top por vista
union_sel <- W_thr %>%
  group_by(view, feature) %>%
  summarise(sel_any = any(selected), .groups="drop") %>%
  filter(sel_any)

inter_sel <- map_dfr(views, function(v) {
  v<- "transcriptomica"
  fset <- top_factors_by_view[[v]]
  dfv  <- W_thr %>% filter(view==v)
  dfv %>%
    group_by(feature) %>%
    summarise(sel_all = all(selected), .groups="drop") %>%
    filter(sel_all) %>%
    mutate(view=v)
})

# =======================================================
# PCA de la reconstrucción completa y de seleccionadas
# =======================================================
Xhat_list <- reconstruct_views(modelo_final)      # lista: vista -> matriz n x p
meta <- modelo_final@samples_metadata


# 1) PCA por vista sobre toda la reconstrucción
pca_all <- lapply(names(Xhat_list), function(v) {
  pc <- do_pca(Xhat_list[[v]])
  list(
    view = v,
    pca  = pc,
    score = plot_scores(pc, meta, "grupo", paste0("Scores PCA (recon) - ", v)),
    loading = plot_loadings(pc, 25, paste0("Loadings PC1 (recon) - ", v)),
    scree = plot_scree(pc, paste0("Scree (recon) - ", v))
  )
}) |> setNames(names(Xhat_list))

# 2) PCA por vista usando features seleccionadas (intersección de top factores)
pca_inter <- lapply(names(Xhat_list), function(v) {
  # Corregido: usar inter_sel directamente, no pipeado desde Xhat_list
  feats <- inter_sel %>% 
    filter(view == v) %>% 
    pull(feature)
  
  if (length(feats) < 2) return(NULL)
  
  Xv <- Xhat_list[[v]]
  feats <- intersect(feats, colnames(Xv))
  
  if (length(feats) < 2) return(NULL)
  
  pc <- do_pca(Xv[, feats, drop = FALSE])
  
  list(
    view = v,
    features = feats,
    pca = pc,
    score = plot_scores(pc, meta, "grupo", paste0("Scores PCA (intersección) - ", v)),
    loading = plot_loadings(pc, 25, paste0("Loadings PC1 (intersección) - ", v)),
    scree = plot_scree(pc, paste0("Scree (intersección) - ", v))
  )
}) |> setNames(names(Xhat_list))
# 3) Opcional: PCA usando la UNION de seleccionadas
pca_union <- lapply(names(Xhat_list), function(v) {
  feats <- union_sel %>% filter(view==v, sel_any) %>% pull(feature)
  if (length(feats) < 2) return(NULL)
  Xv <- Xhat_list[[v]]
  feats <- intersect(feats, colnames(Xv))
  if (length(feats) < 2) return(NULL)
  pc <- do_pca(Xv[, feats, drop=FALSE])
  list(
    view = v,
    features = feats,
    pca  = pc,
    score = plot_scores(pc, meta, "grupo", paste0("Scores PCA (unión) - ", v)),
    loading = plot_loadings(pc, 25, paste0("Loadings PC1 (unión) - ", v)),
    scree = plot_scree(pc, paste0("Scree (unión) - ", v))
  )
}) |> setNames(names(Xhat_list))

# =======================================================
# Gráficos de pesos absolutos con umbral (ejemplos)
# =======================================================
# Cambia "view" y "factor" según corresponda
feats_tr_f1 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[1], factor_sel = "Factor1",n_max = 30))
feats_tr_f2 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[1], factor_sel = "Factor2",n_max = 30))


feats_tr <- unique(c(feats_tr_f1,feats_tr_f2))

feats_pr_f1 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[2], factor_sel = "Factor1",n_max = 30))
feats_pr_f2 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[2], factor_sel = "Factor2",n_max = 30))

feats_pr <- unique(c(feats_pr_f1,feats_pr_f2))

feats_me_f1 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[3], factor_sel = "Factor1",n_max = 30))
feats_me_f2 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[3], factor_sel = "Factor2",n_max = 30))

feats_me <- unique(c(feats_me_f1,feats_me_f2))

feats_cl_f1 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[4], factor_sel = "Factor1",n_max = 30))
feats_cl_f2 <- (plot_abs_weights_with_thr(W_thr, view_sel = views[4], factor_sel = "Factor2",n_max = 30))

feats_cl <- unique(c(feats_cl_f1,feats_cl_f2))

(feats_view <- unique(c(feats_tr,feats_pr,feats_me,feats_cl)))


# --- helper: renombra columnas Entrez -> GeneSymbol solo en transcriptomica
rename_transcrip_if_needed <- function(mat, view, fm = rfm$features_metadata) {
  if (!identical(view, "transcriptomica")) return(mat)
  stopifnot(all(c("EntrezGeneID","GeneSymbol") %in% names(fm)))
  cols <- colnames(mat)
  idx  <- match(cols, fm$EntrezGeneID)
  new  <- cols
  hit  <- !is.na(idx)
  new[hit] <- ifelse(is.na(fm$GeneSymbol[idx[hit]]), cols[hit], fm$GeneSymbol[idx[hit]])
  # Desambigua duplicados: SYMBOL|ENTREZ para los que colisionan
  dup <- duplicated(new) | duplicated(new, fromLast = TRUE)
  new[dup] <- paste0(new[dup], "|", cols[dup])
  colnames(mat) <- make.unique(new, sep = "_")
  mat
}

# --- aplica renombrado por vista
Xhat_list_renamed <- Map(function(x, nm) rename_transcrip_if_needed(x, nm), 
                         Xhat_list, names(Xhat_list))

# --- filtra por feats_view y concatena
X_hat_view <- Reduce(
  "cbind",
  lapply(Xhat_list_renamed, function(x) {
    keep <- intersect(colnames(x), feats_view)
    x[, keep, drop = FALSE]
  })
)

# --- respeta el orden de feats_view
ord <- intersect(feats_view, colnames(X_hat_view))
X_hat_view <- X_hat_view[, ord, drop = FALSE]
# --- respeta el orden de feats_view
ord <- intersect(feats_view, colnames(X_hat_view))
X_hat_view <- X_hat_view[, ord, drop = FALSE]

# --- helper robusto: obtiene colnames del mismo objeto usado en el cbind
colnames_safe <- function(x) {
  if (is.list(x) && !(is.data.frame(x) || is.matrix(x))) x <- x[[1L]]
  colnames(x)
}
# 0) Partimos de X_hat_view ya filtrada y ordenada por 'feats_view'
feats_in_mat <- colnames(X_hat_view)  # sujeto x feature

# 1) Mapa feature -> ómica usando EXACTAMENTE los colnames de cada vista renombrada
feat2omic <- setNames(rep(NA_character_, length(feats_in_mat)), feats_in_mat)
for (nm in names(Xhat_list_renamed)) {
  cn <- intersect(feats_in_mat, colnames(Xhat_list_renamed[[nm]]))
  if (length(cn)) feat2omic[cn] <- nm
}

# Sanidad
unmapped <- names(feat2omic)[is.na(feat2omic)]
if (length(unmapped)) warning("Features sin mapeo de ómica: ",
                              paste(unmapped, collapse = ", "))

# 2) Anotación como factor 4 niveles en el orden de la lista
feat_annot <- data.frame(
  feature = names(feat2omic),
  omic    = factor(unname(feat2omic), levels = names(Xhat_list_renamed)),
  stringsAsFactors = FALSE
)

# 3) Matriz larga: filas = sujeto×feature; columnas finales incluyen 'omic'
#    Cumple “filas sujetos, columna final nombre de la ómica”
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

X_long <- as.data.frame(X_hat_view) %>%
  tibble::rownames_to_column("sujeto") %>%
  tidyr::pivot_longer(-sujeto, names_to = "feature", values_to = "valor") %>%
  dplyr::left_join(feat_annot, by = "feature") %>%
  dplyr::relocate(omic, .after = dplyr::last_col())

# 4) Si necesitas la versión “features en filas” con omic como columna (sin perder info):
X_feat_rows <- as.data.frame(t(X_hat_view))
X_feat_rows$feature <- rownames(X_feat_rows)
X_feat_rows$omic <- feat_annot$omic[match(X_feat_rows$feature, feat_annot$feature)]
# rownames se mantienen; no se pisa la matriz original.

# 5) Si quieres conservar el mapeo junto a la matriz wide sin mezclar tipos:
attr(X_hat_view, "feature_omic") <- feat_annot$omic  # vector factor alineado a colnames(X_hat_view)
