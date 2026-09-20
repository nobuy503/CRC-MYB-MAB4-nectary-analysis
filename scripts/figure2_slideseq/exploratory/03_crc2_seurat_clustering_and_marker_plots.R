# Exploratory Slide-seq analysis of the crc_2 library
#
# This script performs Seurat normalization, clustering and visualization for
# the independently prepared crc_2 Slide-seq library. Spatial expression maps
# for TPS24, AT1G62480, LHCB1.1 and AT1G56610 correspond to the analyses shown
# in Fig. 2p-s. Additional genes were inspected during tissue annotation.
# The Curio Seeker Pipeline-generated Seurat object is the input.

library(Seurat)
library(ggplot2)
library(viridis)

input_file <- "data/processed/crc_2_seurat.rds"
crc_2_seurat <- readRDS(input_file)
DefaultAssay(crc_2_seurat) <- "RNA"

# Inspect the vendor-generated object and its spatial coordinates.
print(crc_2_seurat)
head(crc_2_seurat@meta.data)
SpatialDimPlot(crc_2_seurat)

# Genes inspected in the crc spatial data.
# The first four genes are displayed in Fig. 2p-s.
figure_markers <- c("TPS24", "AT1G62480", "LHCB1.1", "AT1G56610")
exploratory_markers <- c(
  "SWEET9", "LHCB1.3", "RBCS-1A", "AT3G41768"
)
spatial_markers <- c(figure_markers, exploratory_markers)

for (gene in spatial_markers) {
  print(
    SpatialFeaturePlot(crc_2_seurat, features = gene) +
      scale_fill_viridis(option = "viridis") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  )
}

# Seurat normalization, dimensionality reduction and clustering.
crc_2_seurat <- NormalizeData(crc_2_seurat)
crc_2_seurat <- FindVariableFeatures(crc_2_seurat)
crc_2_seurat <- ScaleData(crc_2_seurat)
crc_2_seurat <- RunPCA(crc_2_seurat, npcs = 30)
crc_2_seurat <- FindNeighbors(crc_2_seurat, dims = 1:10)
crc_2_seurat <- FindClusters(crc_2_seurat, resolution = 0.5)
crc_2_seurat <- RunUMAP(crc_2_seurat, dims = 1:10)

print(DimPlot(crc_2_seurat, reduction = "umap", label = TRUE, label.size = 5))
print(table(crc_2_seurat$seurat_clusters))
print(SpatialDimPlot(crc_2_seurat, group.by = "seurat_clusters"))

# UMAP expression plots used to inspect marker distribution across clusters.
umap_markers <- c("SWEET9", "TPS24", "AT1G62480", "LHCB1.1")

for (gene in umap_markers) {
  print(
    FeaturePlot(
      crc_2_seurat,
      features = gene,
      reduction = "umap",
      order = TRUE,
      cols = c("grey95", "blue4"),
      min.cutoff = 0,
      max.cutoff = 2,
      pt.size = 1
    )
  )
}

# Optional processed-object output for subsequent cross-library analysis.
# saveRDS(crc_2_seurat, "data/processed/crc_2_seurat_clustered.rds")
