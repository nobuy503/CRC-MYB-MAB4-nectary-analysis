# Supplementary Fig. 3: Slide-seq replicate from a wild-type nectary
#
# This script analyzes the independently prepared WT_5 Slide-seq library.
# In this library, Seurat cluster 6 was annotated as the nectary-enriched
# cluster. The Curio Seeker Pipeline-generated Seurat object is the input.
#
# Panel correspondence:
#   a: spatial marker maps (LHCB1.1, PI, PSAO and CRC); the tissue schematic
#      and anatomical outlines were added during figure assembly
#   b: UMAP colored by Seurat cluster
#   c-f: UMAP expression plots for CRC, SWEET9, TPS24 and AT1G62480
#   g: spatial map of Seurat clusters
#   h-k: spatial expression maps for CRC, SWEET9, TPS24 and AT1G62480

library(Seurat)
library(ggplot2)
library(viridis)

input_file <- "data/processed/WT_5_seurat.rds"
output_dir <- "results/supplementary_figure_3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

WT_5_seurat <- readRDS(input_file)
DefaultAssay(WT_5_seurat) <- "RNA"

# Apply the library-specific Seurat workflow.
WT_5_seurat <- NormalizeData(
  WT_5_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
WT_5_seurat <- FindVariableFeatures(WT_5_seurat)
WT_5_seurat <- ScaleData(WT_5_seurat)
WT_5_seurat <- RunPCA(WT_5_seurat, npcs = 30)
WT_5_seurat <- FindNeighbors(WT_5_seurat, dims = 1:10)
WT_5_seurat <- FindClusters(WT_5_seurat, resolution = 0.5)
WT_5_seurat <- RunUMAP(WT_5_seurat, dims = 1:10)

# Confirm cluster identities and the presence of plotted genes.
print(table(WT_5_seurat$seurat_clusters))

tissue_markers <- c(
  photosynthetic_tissue = "LHCB1.1",
  petal = "PI",
  sepal = "PSAO",
  nectary = "CRC"
)
nectary_markers <- c("CRC", "SWEET9", "TPS24", "AT1G62480")
required_genes <- unique(c(tissue_markers, nectary_markers))

missing_genes <- setdiff(required_genes, rownames(WT_5_seurat))
if (length(missing_genes) > 0) {
  stop("Genes absent from the WT_5 object: ", paste(missing_genes, collapse = ", "))
}

# Supplementary Fig. 3b: cluster-level UMAP.
p_umap_clusters <- DimPlot(
  WT_5_seurat,
  reduction = "umap",
  group.by = "seurat_clusters"
)
ggsave(
  file.path(output_dir, "supp_fig3b_cluster_umap.pdf"),
  p_umap_clusters,
  width = 6,
  height = 5
)

# Supplementary Fig. 3c-f: nectary-marker expression on the UMAP.
for (gene in nectary_markers) {
  p <- FeaturePlot(
    WT_5_seurat,
    features = gene,
    reduction = "umap",
    order = TRUE,
    cols = c("grey95", "blue4"),
    min.cutoff = 0,
    max.cutoff = 2,
    pt.size = 1
  )

  ggsave(
    file.path(output_dir, paste0("supp_fig3_umap_", gene, ".pdf")),
    p,
    width = 6,
    height = 5
  )
}

# Supplementary Fig. 3g: spatial map of Seurat clusters.
p_spatial_clusters <- SpatialDimPlot(
  WT_5_seurat,
  group.by = "seurat_clusters"
)
ggsave(
  file.path(output_dir, "supp_fig3g_spatial_clusters.pdf"),
  p_spatial_clusters,
  width = 6,
  height = 5
)

# Supplementary Fig. 3a: tissue-annotation marker maps.
for (gene in unname(tissue_markers)) {
  p <- SpatialFeaturePlot(
    WT_5_seurat,
    features = gene,
    min.cutoff = 0,
    max.cutoff = "q90",
    pt.size.factor = 2,
    alpha = c(0.1, 1)
  ) +
    scale_fill_viridis_c(option = "viridis") +
    theme(panel.grid = element_blank())

  ggsave(
    file.path(output_dir, paste0("supp_fig3a_spatial_", gene, ".pdf")),
    p,
    width = 6,
    height = 5
  )
}

# Supplementary Fig. 3h-k: nectary-marker spatial maps.
for (gene in nectary_markers) {
  p <- SpatialFeaturePlot(
    WT_5_seurat,
    features = gene,
    min.cutoff = 0,
    max.cutoff = "q90",
    pt.size.factor = 2,
    alpha = c(0.1, 1)
  ) +
    scale_fill_viridis_c(option = "viridis") +
    theme(panel.grid = element_blank())

  ggsave(
    file.path(output_dir, paste0("supp_fig3_spatial_", gene, ".pdf")),
    p,
    width = 6,
    height = 5
  )
}

# Cluster 6 was the nectary-enriched cluster in the WT_5 replicate.
Idents(WT_5_seurat) <- "seurat_clusters"
cluster6_markers <- FindMarkers(
  WT_5_seurat,
  ident.1 = "6",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.05,
  logfc.threshold = 0.10
)
cluster6_markers$gene <- rownames(cluster6_markers)

write.csv(
  cluster6_markers,
  file.path(output_dir, "WT_5_cluster6_positive_markers.csv"),
  row.names = FALSE
)
