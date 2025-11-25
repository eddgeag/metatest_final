library(MOFA2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(factoextra)
library(tibble)
library(ggrepel)
library(Hmisc)   
library(limma)
library(Bolstad)
library(tidyr)
library(ggnewscale)

plot_pca_feats <- function(pcx, df_weights_unique, feats, axes = c(1,2)) {
  
  # Paleta fija
  omics_colors <- c(
    "Transcriptomics" = "red",
    "Proteomics"      = "black",
    "Metabolomics"    = "olivedrab",
    "Clinical"        = "blue"
  )
  
  # Extraer loadings de PCA
  loadings <- as.data.frame(pcx$rotation[, axes]) %>%
    tibble::rownames_to_column("Feature") %>%
    filter(Feature %in% feats) %>%
    left_join(df_weights_unique, by = "Feature") %>%
    filter(!is.na(view))   # quitamos las que no tienen ómica asignada
  
  # Renombrar ejes
  axis_names <- paste0("PC", axes)
  colnames(loadings)[2:3] <- axis_names
  
  varianza <-  round(100*(pcx$sdev^2)/sum(pcx$sdev^2),2)
  # --- Plot ---
  ggplot(loadings, aes_string(
    x = axis_names[1],
    y = axis_names[2],
    label = "Feature",
    color = "view")) +
    geom_point(size = 3) +
    geom_text_repel(max.overlaps = 50, size = 3) +
    scale_color_manual(values = omics_colors) +
    labs(
      title = "PCA loadings",
      x = paste0(axis_names[1], " (", varianza[as.numeric(gsub("PC","",axis_names[1]))], "%)"),
      y = paste0(axis_names[2], " (", varianza[as.numeric(gsub("PC","",axis_names[2]))], "%)"),
      color = "Ómic"
    ) +
    theme_minimal(base_size = 14)
}


get_corr <- function(df, feats) {
  cor(df[, c("Factor1", "Factor2")], df[, feats, drop = FALSE], method = "spearman")
}

make_corr_heatmap_global <- function(data, annotation_col, ann_colors) {
  sub <- data %>% select(-y)
  feats <- intersect(annotation_col$Label, colnames(sub))
  
  # orden por vista
  ann_sub <- annotation_col[feats, ]
  ann_sub <- ann_sub[order(ann_sub$View), , drop = FALSE]
  feats <- rownames(ann_sub)
  
  corr_sub <- get_corr(sub, feats)
  colnames(corr_sub) <- feats
  
  pheatmap(
    corr_sub,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 12,
    fontsize_col = 7,
    main = "Global Spearman correlation (Factors vs Features)",
    annotation_col = ann_sub["View", drop = FALSE],
    annotation_colors = ann_colors,
    border_color = NA
  )
}

make_corr_heatmap_group <- function(data,
                                    group_name,
                                    annotation_col,
                                    ann_colors) {
  sub <- data %>% filter(y == group_name) %>% select(-y)
  
  feats <- intersect(annotation_col$Label, colnames(sub))
  
  ann_sub <- annotation_col[feats, ]
  ann_sub <- ann_sub[order(ann_sub$View), , drop = FALSE]
  feats <- rownames(ann_sub)
  
  corr_sub <- get_corr(sub, feats)
  colnames(corr_sub) <- feats
  
  pheatmap(
    corr_sub,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 12,
    fontsize_col = 7,
    main = paste("Spearman correlation:", group_name),
    annotation_col = ann_sub["View", drop = FALSE],
    annotation_colors = ann_colors,
    border_color = NA
  )
}
dir_create <- function(x) {  dir.create(x, recursive = T)}
normalize_name <- function(x) {
  
  x %>%
    gsub("_", "", .) %>%
    gsub("-", "", .) %>%
    gsub("\\.", "", .) %>%
    toupper()
}

get_corr_sig <- function(df, feats) {
  mat <- as.matrix(df[, feats, drop = FALSE])
  res <- Hmisc::rcorr(mat, type = "spearman")
  corr <- res$r
  pval <- res$P
  sig <- ifelse(pval < 0.05, "*", "")
  diag(corr) <- 1
  diag(sig)  <- ""
  list(corr = corr, sig = sig)
  
}

# --- Heatmap sin dendrograma ---
# --- Heatmap sin dendrograma ---
make_feature_corr_heatmap <- function(data, annotation_col, ann_colors, group_name = "Global") {
  
  if (group_name != "Global") {
    sub <- data %>% filter(y == group_name)
  } else {
    sub <- data
  }
  
  # Eliminar columnas solo si existen
  drop_cols <- c("y", "Factor1", "Factor2")
  sub <- sub %>% select(-any_of(drop_cols))
  
  feats <- intersect(annotation_col$Label, colnames(sub))
  
  # ordenar por vista
  ann_sub <- annotation_col[feats, , drop = FALSE]
  ann_sub <- ann_sub[order(ann_sub$View), , drop = FALSE]
  feats <- rownames(ann_sub)
  
  corr_res <- get_corr_sig(sub, feats)
  
  pheatmap(
    corr_res$corr[feats, feats],
    color = colorRampPalette(c("blue","white","red"))(100),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 6,
    fontsize_col = 6,
    main = paste("Spearman correlations between features -", group_name),
    annotation_col = ann_sub["View", drop = FALSE],
    annotation_colors = ann_colors,
    display_numbers = corr_res$sig[feats, feats],
    number_color = "black",
    border_color = NA
  )
}


contrast_t_test <- function(datas, factor_, response, contrast_expr) {
  
 
  # Verificar que el factor y la respuesta están en el dataframe
  if (!all(c(factor_, response) %in% colnames(datas))) {
    stop("El factor y la variable de respuesta deben estar en las columnas del dataframe.")
  }
  
  # Extraer los niveles del factor y crear la matriz de diseño
  levels_factor <- levels(datas[[factor_]])
  design <- model.matrix( ~ 0 + datas[[factor_]])  # Diseño sin intercepto
  
  # Asignar nombres a las columnas de la matriz de diseño
  colnames(design) <- levels_factor
  
  # Crear la matriz de contrastes con makeContrasts
  contrasts <- makeContrasts(contrasts = contrast_expr, levels = levels(datas[[factor_]]))
  
  # Verificar que la matriz de contrastes sea válida
  if (nrow(contrasts) != ncol(design)) {
    stop("El contraste especificado no coincide con los niveles del factor.")
  }
  
  # Aplicar el contraste al diseño para obtener grupos
  contrast_values <- design %*% contrasts
  
  # Extraer los grupos contrastados para el t-test
  x <- datas[[response]][contrast_values[, 1] > 0]
  y <- datas[[response]][contrast_values[, 1] < 0]
  
  
  
  # Compute empirical estimates
  m_x <- mean(x)
  m_y <- mean(y)
  v_x <- var(x)
  v_y <- var(y)
  
  # Set prior parameters based on empirical estimates
  m <- c(m_x, m_y)            # Prior means
  n0 <- c(1 / v_x, 1 / v_y)       # Prior precisions (inverse of variances)
  sig.med <- median(c(sd(x), sd(y)))  # Prior median standard deviation
  kappa <- 1                  # Degrees of freedom for prior on sigma
  
  # Perform Bayesian t-test with empirical priors
  t_test_result <- bayes.t.test(
    x,
    y,
    var.equal = F,
    prior = "joint.conj",
    m = m,
    n0 = n0,
    sig.med = sig.med,
    kappa = kappa
  )
  
  
  return(t_test_result)
}


plot_heatmap_stats_ordered <- function(statistic, p_adj, df_weights, alpha = 0.05) {
  # Pasar a formato largo
  stat_long <- statistic %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Contrast", values_to = "Statistic")
  
  pval_long <- p_adj %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Contrast", values_to = "p_adj")
  
  df_long <- left_join(stat_long, pval_long, by = c("Feature", "Contrast"))
  df_long <- right_join(df_long, df_weights, by = "Feature")
  df_long <- df_long %>%
    mutate(
      sig = ifelse(p_adj < alpha, "*", ""),
      view = factor(view, levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Clinical")),
      Feature = factor(Feature, levels = df_weights$Feature)
    )
  
  # Colores de ómicas
  omics_colors <- c(
    "Transcriptomics" = "red",
    "Proteomics"      = "black",
    "Metabolomics"    = "olivedrab",
    "Clinical"        = "blue"
  )
  
  # Data frame auxiliar para anotación de ómica
  df_bar <- df_long %>% distinct(Feature, view)
  
  ggplot(df_long, aes(x = Contrast, y = Feature)) +
    # Heatmap principal
    geom_tile(aes(fill = Statistic), color = "white") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "Statistic"
    ) +
    
    # Estrellitas
    geom_text(aes(label = sig), color = "black", size = 3) +
    
    # Barra lateral con colores de ómica (otra escala de fill, independiente)
    new_scale_fill() +  # <- requiere library(ggnewscale)
    geom_tile(
      data = df_bar,
      aes(x = 0, y = Feature, fill = view),
      inherit.aes = FALSE,
      width = 0.3
    ) +
    scale_fill_manual(values = omics_colors, name = "Omic") +
    
    facet_grid(view ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text.y = element_text(angle = 0, face = "bold"),
      panel.spacing.y = unit(0.5, "lines")
    ) +
    labs(
      title = "Empirical Bayes t-test",
      x = "Contrast", y = "Feature"
    )
}

modelo <- readRDS("./modelo.rds")
biomarcadores <- readRDS("./biomarcadores_7.rds")
trans_metadata <- readRDS("~/Documents/PhD/metatest_final/repo/ready_for_modeling.rds")$features_metadata
df_train <- readRDS("~/Documents/PhD/metatest_final/repo/df_train_final.rds")
df_test <- readRDS("~/Documents/PhD/metatest_final/repo/df_test_final.rds")

feats_included <- colnames(df_train)
all_data <- bind_rows(df_train, df_test)

folder <- "./27-09-25_7"

dir_create(folder)

##===plots_mofa2=====

###===variance explained====

folder_integracion <- "integracion"
dir_create(file.path(folder, folder_integracion))

ps <- plot_variance_explained(modelo, plot_total = T)

df_factors <- ps[[1]]$data %>%
  mutate(
    factor = recode(factor, "Factor1" = "Factor 1", "Factor2" = "Factor 2"),
    view = recode(
      view,
      "transcriptomica" = "Transcriptomics",
      "proteomica"      = "Proteomics",
      "metabolomica"    = "Metabolomics",
      "clinical"        = "Clinical"
    )
  )

df_total <- ps[[2]]$data %>%
  mutate(
    view = recode(
      view,
      "transcriptomica" = "Transcriptomics",
      "proteomica"      = "Proteomics",
      "metabolomica"    = "Metabolomics",
      "clinical"        = "Clinical"
    )
  )

# Definir paleta
omics_colors <- c(
  "Clinical"       = "blue",
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab"
)

# p1
p1 <- ggplot(df_factors, aes(x = factor, y = value, fill = view)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = omics_colors) +
  labs(x = "Factor", y = "Explained variance (%)", fill = "Omics view") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# p2
p2 <- ggplot(df_total, aes(x = view, y = R2, fill = view)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = omics_colors) +
  labs(x = "Omics view", y = "Total explained variance (%)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

fig <- (p1 + p2) + plot_annotation(tag_levels = "A")

fig


ggsave(
  filename = file.path(folder, folder_integracion, "./Variance_Explained.jpeg"),
  plot = fig,
  width = 13,
  height = 8
)


###====factores=======

Z <- as.data.frame(get_factors(modelo)[[1]])
Z$grupo <- factor(modelo@samples_metadata$grupo)
Z$sample <- modelo@samples_metadata$sample  # etiquetas

p_scatter <- ggplot(Z, aes(
  x = Factor1,
  y = Factor2,
  color = grupo,
  label = sample
)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3,
                  max.overlaps = 20,
                  show.legend = FALSE) +
  labs(
    title = "MOFA2 Latent Factors",
    x = "Factor 1",
    y = "Factor 2",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Mostrar
p_scatter

ggsave(
  filename = file.path(folder, folder_integracion, "./MOFA2_latent_factors.jpeg"),
  plot = fig,
  width = 8,
  height = 8
)

###====pesos=======


aux <- plot_top_weights(modelo,
                        view = 1:4,
                        factors = 1:2,
                        nfeatures = Inf)
aux$data$feature

# Definir colores consistentes
omics_colors <- c(
  "Clinical"        = "blue",
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab"
)

# --- Preparar df_weights ---
# --- Preparar df_weights ---

df_weights <- aux$data %>%
  # Normalizar nombres de vista
  mutate(
    view = recode(
      view,
      "transcriptomica" = "Transcriptomics",
      "proteomica"      = "Proteomics",
      "metabolomica"    = "Metabolomics",
      "clinical"        = "Clinical"
    ),
    weight_signed = ifelse(sign == "+", value, -value),
    feature_id = as.character(feature)
  ) %>%
  # Mapear transcriptómicas con GeneSymbol (usando metadata deduplicada)
  left_join(
    trans_metadata %>%
      distinct(EntrezGeneID, .keep_all = TRUE) %>%
      transmute(EntrezGeneID = as.character(EntrezGeneID), GeneSymbol),
    by = c("feature_id" = "EntrezGeneID")
  ) %>%
  mutate(label = case_when(
    view == "Transcriptomics" & !is.na(GeneSymbol) ~ GeneSymbol,
    TRUE ~ feature_id
  ))

# Mantener solo las que están en all_data
df_weights <- df_weights[df_weights$label %in% colnames(all_data), ]

# Factor 1
df_factor1 <- df_weights %>%
  filter(factor == "Factor1") %>%
  group_by(label) %>%
  slice_max(order_by = abs(weight_signed), n = 1) %>%  # quedarnos con el peso dominante
  ungroup() %>%
  arrange(weight_signed) %>%
  mutate(order = row_number())

p1 <- ggplot(df_factor1, aes(x = order, y = weight_signed, fill = view)) +
  geom_col() +
  scale_fill_manual(values = omics_colors) +
  coord_flip() +
  scale_x_continuous(breaks = df_factor1$order, labels = df_factor1$label) +
  labs(title = "MOFA2 Feature Weights (Factor 1)",
       x = "Feature (ordered)",
       y = "Signed weight",
       fill = "Omics view") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ===================== #
# Factor 2
df_factor2 <- df_weights %>%
  filter(factor == "Factor2") %>%
  group_by(label) %>%
  slice_max(order_by = abs(weight_signed), n = 1) %>%
  ungroup() %>%
  arrange(weight_signed) %>%
  mutate(order = row_number())

p2 <- ggplot(df_factor2, aes(x = order, y = weight_signed, fill = view)) +
  geom_col() +
  scale_fill_manual(values = omics_colors) +
  coord_flip() +
  scale_x_continuous(breaks = df_factor2$order, labels = df_factor2$label) +
  labs(title = "MOFA2 Feature Weights (Factor 2)",
       x = "Feature (ordered)",
       y = "Signed weight",
       fill = "Omics view") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ===================== #
# Figura compuesta
fig <- p1 + p2 + plot_annotation(tag_levels = "A")
fig

ggsave(
  filename = file.path(folder, folder_integracion, "./MOFA2_used_weights.jpeg"),
  plot = fig,
  width = 17,
  height = 15
)





##===bivariante=====
# Colores fijos
omics_colors <- c(
  "Clinical"        = "blue",
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab"
)

valid_feats <- intersect(df_weights$label, colnames(all_data))

annotation_col <- df_weights %>%
  filter(label %in% valid_feats) %>%
  distinct(label, .keep_all = TRUE) %>%   # quita duplicados por Factor1/2
  transmute(ID    = feature_id,
            Label = label,
            View  = factor(
              view,
              levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Clinical")
            ))

# rownames = Label para que coincida con los nombres de columnas del heatmap
rownames(annotation_col) <- annotation_col$Label

ann_colors <- list(View = omics_colors)


p1 <- make_corr_heatmap_global(all_data, annotation_col, ann_colors)
p2 <- make_corr_heatmap_group(all_data, "NP", annotation_col, ann_colors)
p3 <- make_corr_heatmap_group(all_data, "SO", annotation_col, ann_colors)
p4 <- make_corr_heatmap_group(all_data, "FD", annotation_col, ann_colors)

folder_correlaciones <- "correlaciones"
dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(
    folder,
    folder_correlaciones,
    "./correlaciones_globales.jpeg"
  ),
  plot = p1,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_NP.jpeg"),
  plot = p2,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_SO.jpeg"),
  plot = p3,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_FD.jpeg"),
  plot = p4,
  width = 10,
  height = 8
)




# --- Colores de ómicas ---
omics_colors <- c(
  "Transcriptomics" = "red",
  "Proteomics"      = "black",
  "Metabolomics"    = "olivedrab",
  "Clinical"        = "blue"
)

# --- Anotación de features ---
# --- Anotación de features ---
valid_feats <- setdiff(colnames(all_data), c("y","Factor1","Factor2"))
annotation_col <- df_weights %>%
  filter(label %in% valid_feats) %>%
  distinct(label, .keep_all = TRUE) %>%
  transmute(
    Label = label,
    View  = factor(view,
                   levels = c("Transcriptomics","Proteomics","Metabolomics","Clinical"))
  )
rownames(annotation_col) <- annotation_col$Label
ann_colors <- list(View = omics_colors)

# --- Función auxiliar: correlación + significancia ---


# --- Ejecutar ---
p_global <- make_feature_corr_heatmap(all_data, annotation_col, ann_colors, "Global")
p_NP     <- make_feature_corr_heatmap(all_data, annotation_col, ann_colors, "NP")
p_SO     <- make_feature_corr_heatmap(all_data, annotation_col, ann_colors, "SO")
p_FD     <- make_feature_corr_heatmap(all_data, annotation_col, ann_colors, "FD")



folder_correlaciones <- "correlaciones_bivariante"
dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(
    folder,
    folder_correlaciones,
    "./correlaciones_globales_sin_integrar.jpeg"
  ),
  plot = p_global,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_NP_sin_integrar.jpeg"),
  plot = p_NP,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_SO_sin_integrar.jpeg"),
  plot = p_SO,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_FD_sin_integrar.jpeg"),
  plot = p_FD,
  width = 10,
  height = 10
)


tst_Data   <- readRDS("./ready_for_modeling.rds")
tst_Data <- tst_Data[1:4]
tst_Data <- lapply(tst_Data, function(x) x$test)
tst_Data <- t(Reduce("rbind",tst_Data))
Z <- get_expectations(modelo,variable = "Z")
W <- get_expectations(modelo,variable = "W")

# --- Función para proyectar nuevos datos al espacio latente de MOFA2 ---
project_to_mofa <- function(tst_Data, W) {
  # tst_Data: matriz (features × muestras)
  # W: lista de pesos de MOFA2, uno por vista
  
  # asegurar que las features estén en filas
  if (nrow(tst_Data) < ncol(tst_Data)) {
    tst_Data <- t(tst_Data)
  }
  
  # Separar tst_Data en lista de vistas usando rownames(W)
  tst_list <- lapply(W, function(Wv) {
    feats <- intersect(rownames(Wv), rownames(tst_Data))
    tst_Data[feats, , drop = FALSE]
  })
  
  # Proyección vista por vista
  Z_estimates <- lapply(names(W), function(vista) {
    X <- tst_list[[vista]]
    Wv <- W[[vista]]
    
    feats <- intersect(rownames(Wv), rownames(X))
    X <- X[feats, , drop = FALSE]
    Wv <- Wv[feats, , drop = FALSE]
    
    # pseudoinversa de Wv
    pinv <- solve(t(Wv) %*% Wv) %*% t(Wv)
    Zv <- t(pinv %*% X)   # muestras × factores
    Zv
  })
  
  # Combinar proyecciones promediando entre vistas
  Z_combined <- Reduce("+", Z_estimates) / length(Z_estimates)
  return(Z_combined)
}

# --- Ejemplo de uso ---
# tst_Data ya cargado (10 × 504)
# W ya calculado con get_expectations(modelo, "W")

W_ <- Reduce("rbind",W)

Z_new <- project_to_mofa(tst_Data, W)

dim(Z_new)  # debería ser 10 muestras × 2 factores
head(Z_new)

R_new <- Z_new %*% t(W_)


R <- MOFA2::predict(modelo)
R <- lapply(R, function(x) x[[1]])
R <- t(Reduce("rbind",R))
R <- as.data.frame(R)


all_R <- as.data.frame(rbind(R,R_new))
# Crear diccionario EntrezID -> GeneSymbol (solo transcriptómica)
map_entrez2symbol <- trans_metadata %>%
  filter(!is.na(EntrezGeneID), !is.na(GeneSymbol)) %>%
  distinct(EntrezGeneID, GeneSymbol)

# Aseguramos que EntrezGeneID sea texto, igual que tus colnames
map_entrez2symbol <- map_entrez2symbol %>%
  mutate(EntrezGeneID = as.character(EntrezGeneID))

# Reemplazar los nombres de las columnas en R
new_names <- colnames(all_R)

# Solo cambiamos donde hay match
new_names <- ifelse(new_names %in% map_entrez2symbol$EntrezGeneID,
                    map_entrez2symbol$GeneSymbol[match(new_names, map_entrez2symbol$EntrezGeneID)],
                    new_names)

colnames(all_R) <- new_names


train_feats <- colnames(df_train)[-c(1:3)]

train_norm <- normalize_name(train_feats)
R_norm     <- normalize_name(colnames(all_R))

# Fuzzy match
matches <- stringdist::amatch(train_norm, R_norm, maxDist = 2)  # tolera pequeñas diferencias

# Crear diccionario train -> R
map_feats <- data.frame(
  df_train = train_feats,
  R_col    = colnames(all_R)[matches]
)

# Mostrar mapeo
print(map_feats)

# Seleccionar las columnas de R en el orden de df_train
R_sub <- all_R[, map_feats$R_col, drop = FALSE]


y_all <- sub("_.*", "", rownames(R_sub))
R_sub$y <-y_all
# --- Ejecutar ---
p_global <- make_feature_corr_heatmap(R_sub, annotation_col, ann_colors, "Global")
p_NP     <- make_feature_corr_heatmap(R_sub, annotation_col, ann_colors, "NP")
p_SO     <- make_feature_corr_heatmap(R_sub, annotation_col, ann_colors, "SO")
p_FD     <- make_feature_corr_heatmap(R_sub, annotation_col, ann_colors, "FD") ## must have more than >4 observation its has 3




folder_correlaciones <- "correlaciones_bivariante_integrado"
dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(
    folder,
    folder_correlaciones,
    "./correlaciones_globales_sin_integrar.jpeg"
  ),
  plot = p_global,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_NP_sin_integrar.jpeg"),
  plot = p_NP,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_SO_sin_integrar.jpeg"),
  plot = p_SO,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(folder, folder_correlaciones, "./correlaciones_FD_sin_integrar.jpeg"),
  plot = p_FD,
  width = 10,
  height = 10
)



##====multivariante ======
# PCA
# PCA
pcx <- prcomp(all_data[, -c(1:3)], scale. = TRUE)

# Contribuciones
axes <- 1:2
p1 <- factoextra::fviz_contrib(pcx, choice = "var", axes = axes)
dd <- facto_summarize(pcx,
                      element = "var",
                      result = "contrib",
                      axes = axes)
contrib <- dd$contrib
names(contrib) <- rownames(dd)

# expected Average contribution
theo_contrib <- 100 / length(contrib)

if (length(axes) > 1) {
  # Ajuste por eigenvalues
  eig <- get_eigenvalue(pcx)[axes, 1]
  theo_contrib <- sum(theo_contrib * eig) / sum(eig)
}

# Subset de features candidatos
final_feats <- biomarcadores$features %>%
  filter(Candidate == "Yes")

feats <- as.character(p1$data$name[p1$data$contrib>=theo_contrib])
final_feats$Feature <-  as.character(final_feats$Feature)
final_feats <- final_feats$Feature
# --- PCA loadings con ómica ---


df_weights_unique <- df_weights %>%
  distinct(label, view)

colnames(df_weights_unique)[1] <- "Feature"

p1 <- plot_pca_feats(pcx, df_weights_unique, feats, axes = c(1,2))

pcx2 <- prcomp(all_data[, final_feats], scale. = TRUE)





p2 <- plot_pca_feats(pcx2, df_weights_unique, final_feats, axes = c(1,2))





# Preparar datos
plt <- data.frame(pcx$x, Group = all_data$y, name = rownames(pcx$x))
varianza <- round(100 * (pcx$sdev^2) / sum(pcx$sdev^2), 2)
# PCA plot 1
p3 <- ggplot(plt, aes(PC1, PC2, color = Group, label = name)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(max.overlaps = 50, size = 3) +
  labs(
    title = "PCA score plot",
    x = paste0("PC1 (", varianza[1], "%)"),
    y = paste0("PC2 (", varianza[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c("NP" = "black", "SO" = "blue", "FD" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )
plt <- data.frame(pcx2$x, Group = all_data$y, name = rownames(pcx2$x))
varianza <- round(100 * (pcx2$sdev^2) / sum(pcx2$sdev^2), 2)
# PCA plot 2
p4 <- ggplot(plt, aes(PC1, PC2, color = Group, label = name)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(max.overlaps = 50, size = 3) +
  labs(
    title = "PCA score plot",
    x = paste0("PC1 (", varianza[1], "%)"),
    y = paste0("PC2 (", varianza[2], "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c("NP" = "black", "SO" = "blue", "FD" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )




folder_multivariante <- "multivariante"
dir_create(file.path(folder, folder_multivariante))

ggsave(
  filename = file.path(folder, folder_multivariante, "./contrib_feats_loadingplot.jpeg"),
  plot = p1,
  width = 10,
  height = 10
)


ggsave(
  filename = file.path(folder, folder_multivariante, "./contrib_final_feats_loadingplot.jpeg"),
  plot = p2,
  width = 10,
  height = 10
)


ggsave(
  filename = file.path(folder, folder_multivariante, "./contrib_score_plot.jpeg"),
  plot = p3,
  width = 10,
  height = 10
)



ggsave(
  filename = file.path(folder, folder_multivariante, "./contrib_score_plot_final.jpeg"),
  plot = p4,
  width = 10,
  height = 10
)


p1
p2
p3
p4
### arreglar para solo hacer el pca con las variables de interes. no seleccinoando


##===univariante=====

###====integrado=====


R_sub$y <- as.factor(R_sub$y)

FD_NP <- lapply(colnames(R_sub)[-ncol(R_sub)], function(col) {
  csa <- contrast_t_test(
    R_sub,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

SO_NP <- lapply(colnames(R_sub)[-ncol(R_sub)], function(col) {
  csa <- contrast_t_test(
    R_sub,
    factor_ = "y",
    response = col,
    contrast_expr = "SO-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})


FD_SO <- lapply(colnames(R_sub)[-ncol(R_sub)], function(col) {
  csa <- contrast_t_test(
    R_sub,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-SO"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

FD_NP <- Reduce("rbind",FD_NP)
rownames(FD_NP) <- colnames(R_sub)[-ncol(R_sub)]

SO_NP <- Reduce("rbind",SO_NP)
rownames(SO_NP) <- colnames(R_sub)[-ncol(R_sub)]

FD_SO <- Reduce("rbind",FD_SO)
rownames(FD_SO) <- colnames(R_sub)[-ncol(R_sub)]


statistic <-  data.frame(FD_NP = FD_NP[,1],
                         SO_NP = SO_NP[,1],
                         FD_SO = FD_SO[,1])

pval <-  data.frame(FD_NP = FD_NP[,2],
                         SO_NP = SO_NP[,2],
                         FD_SO = FD_SO[,2])

p_adj <- apply(pval, 2, function(col) p.adjust(col, method = "BH"))




p1 <- plot_heatmap_stats_ordered(statistic = statistic,
                           p_adj = p_adj ,
                         df_weights = df_weights_unique,
                         alpha = 0.05)

folder_correlaciones <- "univariante_integrado"

dir_create(file.path(folder, folder_correlaciones))

ggsave(
  filename = file.path(folder, folder_correlaciones, "./univariante_integrado.jpeg"),
  plot = p1,
  width = 15,
  height = 10
)


###====integrado=====

all_data_ <- all_data[,-c(2:3)]


FD_NP <- lapply(colnames(all_data_)[-1], function(col) {
  csa <- contrast_t_test(
    all_data_,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

SO_NP <- lapply(colnames(all_data_)[-1], function(col) {
  csa <- contrast_t_test(
    all_data_,
    factor_ = "y",
    response = col,
    contrast_expr = "SO-NP"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})


FD_SO <- lapply(colnames(all_data_)[-1], function(col) {
  csa <- contrast_t_test(
    all_data_,
    factor_ = "y",
    response = col,
    contrast_expr = "FD-SO"
  )
  
  return(c(statistic = csa$statistic,
           p_value = csa$p.value))
  
})

FD_NP <- Reduce("rbind",FD_NP)
rownames(FD_NP) <- colnames(all_data_)[-1]

SO_NP <- Reduce("rbind",SO_NP)
rownames(SO_NP) <- colnames(all_data_)[-1]

FD_SO <- Reduce("rbind",FD_SO)
rownames(FD_SO) <- colnames(all_data_)[-1]


statistic <-  data.frame(FD_NP = FD_NP[,1],
                         SO_NP = SO_NP[,1],
                         FD_SO = FD_SO[,1])

pval <-  data.frame(FD_NP = FD_NP[,2],
                    SO_NP = SO_NP[,2],
                    FD_SO = FD_SO[,2])

p_adj <- apply(pval, 2, function(col) p.adjust(col, method = "BH"))




p1 <- plot_heatmap_stats_ordered(statistic = statistic,
                           p_adj = p_adj ,
                           df_weights = df_weights_unique,
                           alpha = 0.05)


ggsave(
  filename = file.path(folder, folder_correlaciones, "./univariante_sin_integrado.jpeg"),
  plot = p1,
  width = 15,
  height = 10
)
