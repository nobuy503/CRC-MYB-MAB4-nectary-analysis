# Exploratory WT_3 Slide-seq analysis for Fig. 2
#
# Input:
#   data/processed/WT_3_seurat.rds
#
# Main analyses:
# - Inspection of the Curio-generated spatial Seurat object.
# - Spatial expression plots for nectary-associated genes.
# - Seurat normalization, variable-feature selection, scaling and PCA.
# - Neighbor graph construction using PCs 1-10.
# - Clustering at resolution 0.5.
# - UMAP visualization and spatial cluster mapping.
# - UMAP feature plots for SWEET9, TPS24, CRC and AT1G62480.
#
# Figure relationship:
# This source script provides the main WT_3 analysis underlying Fig. 2c, d and
# f-k. Fig. 2b quality mapping, panel e schematic, WT_5 analyses, crc analyses,
# differential-expression analysis, Venn diagram and final dot plot are not
# included.
#
# Provenance note:
# This is the exploratory source used during figure preparation. Repeated
# plotting blocks are retained. A streamlined final script will be prepared
# after all WT and crc source scripts have been identified.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
})

input_rds <- "data/processed/WT_3_seurat.rds"
WT_3_seurat <- readRDS(input_rds)


# オブジェクトの基本情報を確認
print(WT_3_seurat)

# メタデータの確認
head(WT_3_seurat@meta.data)

# 空間データが正しく含まれているか確認
SpatialDimPlot(WT_3_seurat)

# 空間的遺伝子発現のプロット
SpatialFeaturePlot(WT_3_seurat, features = c("SWEET9"))

# 空間的遺伝子発現のプロット
SpatialFeaturePlot(WT_3_seurat, features = c("CRC"))

# 空間的遺伝子発現のプロット
SpatialFeaturePlot(WT_3_seurat, features = c("TPS24"))

# 空間的遺伝子発現のプロット
SpatialFeaturePlot(WT_3_seurat, features = c("AT1G62480"))



# 必要なライブラリの読み込み
library(viridis)
library(ggplot2)

# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("SWEET9"))


# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )





# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("TPS24"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )




# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("CRC"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )




# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("AT1G62480"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )


# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("PRXIIC"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )


# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("AT2G27770"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )



# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("NSP1"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )

# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("SHP1"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )

# SpatialFeaturePlotのプロットを取得
plot <- SpatialFeaturePlot(WT_3_seurat, features = c("KAT5"))

# ggplot2でviridisカラースケールを適用し、グリッドを非表示
plot +
  scale_fill_viridis(option = "viridis") + # viridisカラースケールを適用
  theme_minimal() + # ミニマルなテーマ
  theme(
    panel.grid = element_blank(),         # グリッドを削除
    panel.grid.major = element_blank(),   # メジャーグリッドを削除
    panel.grid.minor = element_blank()    # マイナーグリッドを削除
  )




# データの正規化
WT_3_seurat <- NormalizeData(WT_3_seurat)

# 変動遺伝子の検出
WT_3_seurat <- FindVariableFeatures(WT_3_seurat)

# データのスケーリング
WT_3_seurat <- ScaleData(WT_3_seurat)

# PCAによる次元削減
WT_3_seurat <- RunPCA(WT_3_seurat, npcs = 30)

# 最近傍グラフの構築
WT_3_seurat <- FindNeighbors(WT_3_seurat, dims = 1:10)

# クラスタリングの実行
WT_3_seurat <- FindClusters(WT_3_seurat, resolution = 0.5) # 解像度は適宜調整

# UMAPの計算
WT_3_seurat <- RunUMAP(WT_3_seurat, dims = 1:10)

DimPlot(WT_3_seurat, reduction = "umap", label = TRUE, label.size = 5)  # オプション：テーマの変更

# クラスターラベルを確認
table(WT_3_seurat$seurat_clusters)

SpatialDimPlot(WT_3_seurat, group.by = "seurat_clusters")

FeaturePlot(
  WT_3_seurat,
  features = c("SWEET9"),
  reduction = "umap",
  order = TRUE,  # 発現値の高い順にプロット
  cols = c("grey95", "blue4"),  # カラースケールを淡い灰色から濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2,  # カットオフの最大値を指定
  pt.size = 1  # ポイントサイズを指定
  )

FeaturePlot(
  WT_3_seurat,
  features = c("TPS24"),
  reduction = "umap",
  order = TRUE,  # 発現値の高い順にプロット
  cols = c("grey95", "blue4"),  # カラースケールを淡い灰色から濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2,  # カットオフの最大値を指定
  pt.size = 1  # ポイントサイズを指定
)



FeaturePlot(
  WT_3_seurat,
  features = c("CRC"),
  reduction = "umap",
  order = TRUE,  # 発現値の高い順にプロット
  cols = c("grey95", "blue4"),  # カラースケールを淡い灰色から濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2,  # カットオフの最大値を指定
  pt.size = 1  # ポイントサイズを指定
)



FeaturePlot(
  WT_3_seurat,
  features = c("AT1G62480"),
  reduction = "umap",
  order = TRUE,  # 発現値の高い順にプロット
  cols = c("grey95", "blue4"),  # カラースケールを淡い灰色から濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2,  # カットオフの最大値を指定
  pt.size = 1  # ポイントサイズを指定
)

