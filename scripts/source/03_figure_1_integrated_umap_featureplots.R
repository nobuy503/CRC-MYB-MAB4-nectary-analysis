# Integrated flower and silique snRNA-seq UMAP visualizations for Figure 1
#
# Panels generated or supported by this source script:
#   Fig. 1b: UMAP colored by tissue/sample of origin (orig.ident)
#   Fig. 1c: UMAP colored by Seurat cluster, with cluster 15 emphasized
#   Fig. 1e: CRC (AT1G69180) expression
#   Fig. 1f: MAB4 (AT4G31820) expression
#   Fig. 1h: MYB57 (AT3G01530) expression
#   Fig. 1i: MYB21 (AT3G27810) expression
#   Fig. 1j: At1g65970 expression
#
# This script loads a previously generated integrated Seurat object. It does
# not perform the flower/silique integration itself and does not generate
# Fig. 1d (heatmap) or Fig. 1g (SWEET9).
#
# The input path has been standardized to:
#   data/processed/integrated_data.rds

library(Seurat)
library(ggplot2)
library(patchwork)

# ---- データ読み込み ----
seurat_integrated <- readRDS("data/processed/integrated_data.rds")

# ---- DimPlot（確認用）----
# まずUMAPが正しく埋め込まれているか確認1
DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP overview (integrated data)")

# まずUMAPが正しく埋め込まれているか確認2
DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 0.6
) + ggtitle("UMAP by sample (orig.ident)")


#クラスター8を上に
library(Seurat)
library(ggplot2)

Idents(seurat_integrated) <- seurat_integrated$seurat_clusters

# 通常のDimPlot
p <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.6
)

# 実際に使われた色を取得
used_cols <- ggplot_build(p)$data[[1]]$colour
used_lvls <- ggplot_build(p)$data[[1]]$group
mapping <- data.frame(level = levels(seurat_integrated$seurat_clusters),
                      color = unique(used_cols[order(as.numeric(used_lvls))]))

# クラスタ8の色を抽出
col15 <- mapping$color[mapping$level == "15"]

# クラスタ8のUMAP座標を抽出
emb <- Embeddings(seurat_integrated, "umap")
meta <- seurat_integrated@meta.data
cells15 <- rownames(meta)[meta$seurat_clusters == "15"]
df15 <- data.frame(UMAP_1 = emb[cells15, 1], UMAP_2 = emb[cells15, 2])

# クラスタ8だけ上書き描画（元と同じ色）
p + geom_point(
  data = df15,
  aes(x = UMAP_1, y = UMAP_2),
  inherit.aes = FALSE,
  color = col15,
  size = 1.2,
  alpha = 1
) + ggtitle("Cluster 15 emphasized, Seurat colors preserved")






# ---- FeaturePlot ----
p <- FeaturePlot(
  seurat_integrated,
  features  = "AT1G69180",
  reduction = "umap",
  order     = TRUE,
  cols      = c("grey90", "blue4"),
  min.cutoff = 0,
  max.cutoff = 2,
  pt.size = 0.7                # ← デフォルトは0.5前後、ここを大きく
)
p + scale_color_gradientn(colours = c("grey90", "blue4"), na.value = "grey90")

p <- FeaturePlot(
  seurat_integrated,
  features  = "AT4G31820",
  reduction = "umap",
  order     = TRUE,
  cols      = c("grey90", "blue4"),
  min.cutoff = 0,
  max.cutoff = 2,
  pt.size = 0.7                # ← デフォルトは0.5前後、ここを大きく
)
p + scale_color_gradientn(colours = c("grey90", "blue4"), na.value = "grey90")


p <- FeaturePlot(
  seurat_integrated,
  features  = "AT3G01530",
  reduction = "umap",
  order     = TRUE,
  cols      = c("grey90", "blue4"),
  min.cutoff = 0,
  max.cutoff = 2,
  pt.size = 0.7                # ← デフォルトは0.5前後、ここを大きく
)
p + scale_color_gradientn(colours = c("grey90", "blue4"), na.value = "grey90")

p <- FeaturePlot(
  seurat_integrated,
  features  = "AT3G27810",
  reduction = "umap",
  order     = TRUE,
  cols      = c("grey90", "blue4"),
  min.cutoff = 0,
  max.cutoff = 2,
  pt.size = 0.7                # ← デフォルトは0.5前後、ここを大きく
)
p + scale_color_gradientn(colours = c("grey90", "blue4"), na.value = "grey90")

p <- FeaturePlot(
  seurat_integrated,
  features  = "AT1G65970",
  reduction = "umap",
  order     = TRUE,
  cols      = c("grey90", "blue4"),
  min.cutoff = 0,
  max.cutoff = 2,
  pt.size = 0.7                # ← デフォルトは0.5前後、ここを大きく
)
p + scale_color_gradientn(colours = c("grey90", "blue4"), na.value = "grey90")

