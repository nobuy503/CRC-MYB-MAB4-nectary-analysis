suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)
})

set.seed(123)

input_rds <- "data/processed/GSE226097_flower_230221.rds"
output_dir <- "results/figure1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

flower <- readRDS(input_rds)
stopifnot("seurat_clusters" %in% colnames(flower@meta.data))
stopifnot("umap" %in% names(flower@reductions))
DefaultAssay(flower) <- "RNA"

seurat_cluster <- as.character(flower$seurat_clusters)
cluster_levels <- sort(unique(seurat_cluster))
if (all(!is.na(suppressWarnings(as.numeric(cluster_levels))))) {
  cluster_levels <- cluster_levels[order(as.numeric(cluster_levels))]
}
if (length(cluster_levels) != 5) {
  stop(
    "Fig. 1p requires exactly five Seurat clusters in the Flower object; found: ",
    paste(cluster_levels, collapse = ", ")
  )
}
flower$Fig1p_cluster <- factor(seurat_cluster, levels = cluster_levels)

# Fig. 1p: the five Seurat clusters in the Flower-only dataset.
p_cluster <- DimPlot(
  flower,
  reduction = "umap",
  group.by = "Fig1p_cluster",
  pt.size = 0.5
) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Seurat cluster") +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, "Fig1p_flower_five_Seurat_clusters.pdf"), p_cluster,
       width = 5.5, height = 4.5)

counts_flower <- GetAssayData(flower, assay = "RNA", layer = "counts")
cell_metadata <- flower@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(counts_flower),
  row.names = rownames(counts_flower)
)

cds <- new_cell_data_set(
  expression_data = counts_flower,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
cds <- preprocess_cds(cds, num_dim = 30)

# Use the same UMAP coordinates as Fig. 1p, then learn the Monocle3 graph.
seurat_umap <- Embeddings(flower, reduction = "umap")
seurat_umap <- seurat_umap[colnames(cds), , drop = FALSE]
reducedDims(cds)$UMAP <- seurat_umap
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# Root selection: choose the principal-graph vertex most frequently nearest to
# nuclei with detectable CRC (AT1G69180) counts.
crc_gene <- "AT1G69180"
stopifnot(crc_gene %in% rownames(cds))
crc_cells <- colnames(cds)[counts(cds)[crc_gene, ] > 0]
if (length(crc_cells) == 0) stop("No CRC-expressing nuclei were detected.")
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
crc_vertex <- as.character(closest_vertex[crc_cells, 1])
root_node <- names(which.max(table(crc_vertex)))
message("CRC-defined root principal node: ", root_node)

cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_node)
pt <- pseudotime(cds)

# Fig. 1q: pseudotime on the same UMAP coordinates as Fig. 1p.
p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.5,
  trajectory_graph_segment_size = 0.7
) + theme_classic(base_size = 12)
ggsave(file.path(output_dir, "Fig1q_flower_Monocle3_pseudotime.pdf"), p_pt,
       width = 5.5, height = 4.5)

# Fig. 1r-s: expression relative to the maximum observed value for each gene.
make_expression_umap <- function(gene_id, gene_label, filename) {
  x <- log1p(as.numeric(counts(cds)[gene_id, ]))
  pct_max <- if (max(x) > 0) 100 * x / max(x) else x
  df <- data.frame(
    UMAP1 = reducedDims(cds)$UMAP[, 1],
    UMAP2 = reducedDims(cds)$UMAP[, 2],
    expression = pct_max
  )
  p <- ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(data = filter(df, expression == 0), color = "grey85", size = 0.5) +
    geom_point(
      data = filter(df, expression > 0),
      aes(color = expression), size = 0.7
    ) +
    scale_color_viridis_c(name = "% Max") +
    coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2", title = gene_label) +
    theme_classic(base_size = 12)
  ggsave(file.path(output_dir, filename), p, width = 5.5, height = 4.5)
}

make_expression_umap("AT1G69180", "CRC", "Fig1r_CRC_expression_UMAP.pdf")
make_expression_umap("AT1G65970", "At1g65970", "Fig1s_AT1G65970_expression_UMAP.pdf")

# Fig. 1t: log-transformed expression with separate y-axis scaling by gene.
genes <- c(
  CRC = "AT1G69180",
  MAB4 = "AT4G31820",
  SWEET9 = "AT2G39060",
  MYB57 = "AT3G01530",
  MYB21 = "AT3G27810",
  At1g65970 = "AT1G65970"
)
stopifnot(all(unname(genes) %in% rownames(cds)))
finite_cells <- names(pt)[is.finite(pt)]
expr <- log1p(as.matrix(counts(cds)[unname(genes), finite_cells, drop = FALSE]))
rownames(expr) <- names(genes)

curve_df <- data.frame(pseudotime = pt[finite_cells], t(expr), check.names = FALSE) |>
  pivot_longer(-pseudotime, names_to = "gene", values_to = "expression") |>
  mutate(gene = factor(gene, levels = names(genes)))

p_curves <- ggplot(curve_df, aes(pseudotime, expression, color = gene)) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    linewidth = 1
  ) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  scale_color_manual(values = setNames(viridis(6, end = 0.9), names(genes))) +
  labs(x = "Pseudotime", y = "Log-transformed expression") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", hjust = 0)
  )
ggsave(file.path(output_dir, "Fig1t_gene_expression_along_pseudotime.pdf"),
       p_curves, width = 5, height = 10)

writeLines(root_node, file.path(output_dir, "Fig1p_t_root_principal_node.txt"))
saveRDS(cds, file.path(output_dir, "Fig1p_t_flower_monocle3_cds.rds"))

