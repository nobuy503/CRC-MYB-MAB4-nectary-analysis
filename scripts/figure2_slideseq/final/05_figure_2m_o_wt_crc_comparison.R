# Fig. 2m and 2o: WT-versus-crc Slide-seq comparison
#
# The independently generated WT_3 and crc_2 Curio Seeker Seurat objects are
# used for sample-level integration, marker visualization and spot-level
# differential-expression analysis.
#
# Panel correspondence:
#   m: WT-versus-crc volcano plot
#   o: DotPlot of selected nectary-associated genes in WT and crc

library(Seurat)
library(dplyr)
library(ggplot2)
library(future)

wt_file <- "data/processed/WT_3_seurat.rds"
crc_file <- "data/processed/crc_2_seurat.rds"
output_dir <- "results/figure_2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

WT_3_seurat <- readRDS(wt_file)
crc_2_seurat <- readRDS(crc_file)

# Construct RNA-only Seurat objects from the vendor-generated count matrices.
WT <- CreateSeuratObject(
  counts = GetAssayData(WT_3_seurat, assay = "RNA", layer = "counts"),
  project = "WT"
)
crc <- CreateSeuratObject(
  counts = GetAssayData(crc_2_seurat, assay = "RNA", layer = "counts"),
  project = "crc"
)

WT <- NormalizeData(WT)
WT <- FindVariableFeatures(WT)
crc <- NormalizeData(crc)
crc <- FindVariableFeatures(crc)

WT <- RenameCells(WT, add.cell.id = "WT")
crc <- RenameCells(crc, add.cell.id = "crc")

if (length(intersect(Cells(WT), Cells(crc))) > 0) {
  stop("Duplicated spot identifiers remain after sample-prefix assignment.")
}

# Integrate WT and crc for joint visualization.
plan(multisession, workers = 4)
options(future.globals.maxSize = 50 * 1024^3)

object_list <- list(WT = WT, crc = crc)
anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:20)
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(integrated_data) <- "integrated"
integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:20, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Fig. 2o: expression of selected nectary-associated genes by sample.
DefaultAssay(integrated_data) <- "RNA"
integrated_data <- JoinLayers(integrated_data)

genes_to_plot <- c("CRC", "SWEET9", "TPS24", "AT1G62480")
available_genes <- intersect(genes_to_plot, rownames(integrated_data))
missing_genes <- setdiff(genes_to_plot, available_genes)
if (length(missing_genes) > 0) {
  stop("Genes absent from the integrated object: ",
       paste(missing_genes, collapse = ", "))
}

dot_plot <- DotPlot(
  integrated_data,
  features = available_genes,
  group.by = "orig.ident",
  dot.scale = 10
) +
  scale_colour_gradient(low = "white", high = "blue4") +
  RotatedAxis() +
  labs(
    x = "Genes",
    y = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(output_dir, "fig2o_WT_crc_dotplot.pdf"),
  dot_plot,
  width = 6,
  height = 4
)

# Fig. 2m: spot-level differential-expression analysis.
# FindMarkers(test.use = "wilcox") uses a two-sided Wilcoxon rank-sum test.
# Seurat reports Bonferroni-adjusted P values in p_val_adj.
Idents(integrated_data) <- "orig.ident"

markers <- FindMarkers(
  integrated_data,
  ident.1 = "WT",
  ident.2 = "crc",
  assay = "RNA",
  logfc.threshold = 0,
  test.use = "wilcox"
)

markers <- markers %>%
  mutate(
    gene = rownames(.),
    log10_adjusted_p = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 1 ~ "Higher in WT",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Higher in crc",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  markers,
  file.path(output_dir, "fig2m_WT_vs_crc_FindMarkers_results.csv"),
  row.names = FALSE
)

genes_to_label <- c("AT1G62480", "TPS24", "SWEET9", "CRC")

volcano_plot <- ggplot(
  markers,
  aes(x = avg_log2FC, y = log10_adjusted_p, color = significance)
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_text(
    data = markers %>% filter(gene %in% genes_to_label),
    aes(label = gene),
    color = "black",
    size = 3,
    vjust = -1
  ) +
  scale_color_manual(
    values = c(
      "Higher in WT" = "purple",
      "Higher in crc" = "#008080",
      "Not significant" = "grey"
    )
  ) +
  labs(
    x = "Average log2 fold change (WT versus crc)",
    y = "-log10 adjusted P value"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "top"
  )

ggsave(
  file.path(output_dir, "fig2m_WT_vs_crc_volcano.pdf"),
  volcano_plot,
  width = 6,
  height = 5
)

wt_high_genes <- markers %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

crc_high_genes <- markers %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC)

write.csv(
  wt_high_genes,
  file.path(output_dir, "fig2m_genes_higher_in_WT.csv"),
  row.names = FALSE
)
write.csv(
  crc_high_genes,
  file.path(output_dir, "fig2m_genes_higher_in_crc.csv"),
  row.names = FALSE
)
