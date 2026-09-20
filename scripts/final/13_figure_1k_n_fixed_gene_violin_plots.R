suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

input_rds <- "data/processed/integrated_data.rds"
output_dir <- "results/figure1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
stopifnot("seurat_clusters" %in% colnames(obj@meta.data))
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]])) > 1) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
if (!"data" %in% Layers(obj[["RNA"]])) {
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

genes <- c(
  AT1G65970 = "At1g65970",
  AT3G01530 = "MYB57",
  AT1G55670 = "At1g55670",
  AT1G55330 = "At1g55330"
)
stopifnot(all(names(genes) %in% rownames(obj)))

cluster_levels <- as.character(0:21)
observed <- unique(as.character(obj$seurat_clusters))
stopifnot(all(observed %in% cluster_levels))

expr <- FetchData(
  obj,
  vars = names(genes),
  assay = "RNA",
  layer = "data"
)
plot_df <- expr |>
  mutate(
    cell = rownames(expr),
    cluster = factor(as.character(obj$seurat_clusters[cell]), levels = cluster_levels)
  ) |>
  pivot_longer(
    cols = all_of(names(genes)),
    names_to = "gene_id",
    values_to = "expression"
  ) |>
  mutate(gene = factor(unname(genes[gene_id]), levels = unname(genes)))

cluster_colors <- hue_pal()(length(cluster_levels))
names(cluster_colors) <- cluster_levels

make_violin <- function(gene_label) {
  ggplot(filter(plot_df, gene == gene_label), aes(cluster, expression, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, color = "black") +
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    labs(x = "Seurat cluster", y = "Log-normalized expression", title = gene_label) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "italic"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

plots <- lapply(unname(genes), make_violin)
names(plots) <- c("Fig1k", "Fig1l", "Fig1m", "Fig1n")

for (panel in names(plots)) {
  ggsave(
    file.path(output_dir, paste0(panel, "_violin.pdf")),
    plots[[panel]], width = 4.2, height = 3.2
  )
}

combined <- wrap_plots(plots, nrow = 1)
ggsave(
  file.path(output_dir, "Fig1k_n_fixed_gene_violin_plots.pdf"),
  combined, width = 16.8, height = 3.2
)

