# Exploratory Slide-seq analysis of the WT_5 library
#
# This script applies the same Seurat workflow used for WT_3 to the
# independently prepared WT_5 library. It was used to inspect the spatial
# distribution of nectary markers and the reproducibility of cluster structure.
# The Curio Seeker Pipeline-generated Seurat object is the input.

library(Seurat)
library(ggplot2)
library(viridis)

input_file <- "data/processed/WT_5_seurat.rds"
WT_5_seurat <- readRDS(input_file)

# Inspect the vendor-generated object and its spatial coordinates.
print(WT_5_seurat)
head(WT_5_seurat@meta.data)
SpatialDimPlot(WT_5_seurat)

# Marker genes inspected in the spatial data.
spatial_markers <- c(
  "SWEET9", "CRC", "TPS24", "AT1G62480", "PRXIIC",
  "AT2G27770", "NSP1", "SHP1", "KAT5"
)

for (gene in spatial_markers) {
  print(
    SpatialFeaturePlot(WT_5_seurat, features = gene) +
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
WT_5_seurat <- NormalizeData(WT_5_seurat)
WT_5_seurat <- FindVariableFeatures(WT_5_seurat)
WT_5_seurat <- ScaleData(WT_5_seurat)
WT_5_seurat <- RunPCA(WT_5_seurat, npcs = 30)
WT_5_seurat <- FindNeighbors(WT_5_seurat, dims = 1:10)
WT_5_seurat <- FindClusters(WT_5_seurat, resolution = 0.5)
WT_5_seurat <- RunUMAP(WT_5_seurat, dims = 1:10)

print(DimPlot(WT_5_seurat, reduction = "umap", label = TRUE, label.size = 5))
print(table(WT_5_seurat$seurat_clusters))
print(SpatialDimPlot(WT_5_seurat, group.by = "seurat_clusters"))

# UMAP expression plots for the principal nectary markers examined in WT_5.
umap_markers <- c("SWEET9", "TPS24", "CRC", "AT1G62480")

for (gene in umap_markers) {
  print(
    FeaturePlot(
      WT_5_seurat,
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
# saveRDS(WT_5_seurat, "data/processed/WT_5_seurat_clustered.rds")
