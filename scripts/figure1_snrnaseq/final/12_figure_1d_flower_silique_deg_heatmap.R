suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(pheatmap)
  library(scales)
})

set.seed(123)

input_rds <- "data/processed/integrated_data.rds"
output_dir <- "results/figure1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
stopifnot(all(c("orig.ident", "seurat_clusters") %in% colnames(obj@meta.data)))
DefaultAssay(obj) <- "RNA"

# Seurat v5 objects may retain separate Flower and Silique layers.
if (length(Layers(obj[["RNA"]])) > 1) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
if (!"data" %in% Layers(obj[["RNA"]])) {
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

# Standardize tissue labels without changing the original metadata.
origin <- as.character(obj$orig.ident)
tissue <- ifelse(
  grepl("flower", origin, ignore.case = TRUE), "Flower",
  ifelse(grepl("silique", origin, ignore.case = TRUE), "Silique", NA_character_)
)
stopifnot(!anyNA(tissue), all(c("Flower", "Silique") %in% tissue))
obj$dataset <- factor(tissue, levels = c("Flower", "Silique"))

# Reproduce the DEG definition used for the associated source-data table.
Idents(obj) <- "dataset"
deg <- FindMarkers(
  obj,
  ident.1 = "Flower",
  ident.2 = "Silique",
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
) |>
  as.data.frame() |>
  tibble::rownames_to_column("gene")

top100 <- deg |>
  filter(!is.na(p_val_adj), p_val_adj < 0.05) |>
  arrange(desc(abs(avg_log2FC)), p_val_adj) |>
  slice_head(n = 100)

stopifnot(nrow(top100) == 100)
write.csv(
  top100,
  file.path(output_dir, "Fig1d_top100_flower_vs_silique_DEGs.csv"),
  row.names = FALSE
)

# Log-normalized expression, displayed as a row Z-score.
expr <- LayerData(obj, assay = "RNA", layer = "data")
genes <- intersect(top100$gene, rownames(expr))
stopifnot(length(genes) == 100)

# Order nuclei first by Seurat cluster and then by tissue. All nuclei are used.
cluster_chr <- as.character(obj$seurat_clusters)
cluster_num <- suppressWarnings(as.numeric(cluster_chr))
cluster_order <- if (all(!is.na(cluster_num))) cluster_num else cluster_chr
cell_order <- order(cluster_order, obj$dataset)
cells <- colnames(obj)[cell_order]

mat <- as.matrix(expr[genes, cells, drop = FALSE])
mat_z <- t(scale(t(mat)))
mat_z[!is.finite(mat_z)] <- 0
mat_z[mat_z < -2] <- -2
mat_z[mat_z > 2] <- 2

cluster_levels <- sort(unique(cluster_chr))
if (all(!is.na(suppressWarnings(as.numeric(cluster_levels))))) {
  cluster_levels <- cluster_levels[order(as.numeric(cluster_levels))]
}
annotation_col <- data.frame(
  Cluster = factor(cluster_chr[cell_order], levels = cluster_levels),
  row.names = cells
)
cluster_colors <- hue_pal()(length(cluster_levels))
names(cluster_colors) <- cluster_levels

pheatmap(
  mat_z,
  filename = file.path(output_dir, "Fig1d_flower_silique_top100_DEG_heatmap.pdf"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = list(Cluster = cluster_colors),
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  color = colorRampPalette(c("#313695", "#F7F7F7", "#A50026"))(201),
  breaks = seq(-2, 2, length.out = 202),
  width = 12,
  height = 8
)

