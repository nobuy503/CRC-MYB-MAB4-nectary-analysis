# Legacy flower-and-silique integration and visualization workflow
#
# This earlier analysis script:
#   - loads the flower and silique Seurat objects from NCBI GEO GSE226097;
#   - performs Seurat anchor-based integration;
#   - runs scaling, PCA, UMAP, neighbor finding and clustering;
#   - generates integrated UMAP, FeaturePlot, DotPlot and violin plots;
#   - includes SWEET9 (AT2G39060) expression visualization.
#
# The original file also contained an obsolete combined Flower+Silique
# Monocle3 pseudotime analysis. That section has been intentionally excluded
# because the final manuscript reports a Flower-only pseudotime analysis.
#
# Local Desktop paths were replaced with repository-relative paths.
# This file documents the earlier workflow and will be streamlined against the
# final integrated object before public release.

# ライブラリの読み込み
library(Seurat)
library(dplyr)


# データのロード
GSE226097_flower_230221 <- readRDS("data/processed/GSE226097_flower_230221.rds")
GSE226097_silique_230221 <- readRDS("data/processed/GSE226097_silique_230221.rds")

# 新しい Seurat オブジェクトを作成
silique_new <- CreateSeuratObject(
  counts = GSE226097_silique_230221@assays$RNA@counts, 
  project = "Silique"
)
flower_new <- CreateSeuratObject(
  counts = GSE226097_flower_230221@assays$RNA@counts, 
  project = "Flower"
)

# 正規化
silique_new <- NormalizeData(silique_new)
flower_new <- NormalizeData(flower_new)

# 変動遺伝子を特定
silique_new <- FindVariableFeatures(silique_new)
flower_new <- FindVariableFeatures(flower_new)

# 変動遺伝子の数を確認
cat("Silique variable features: ", length(VariableFeatures(silique_new)), "\n")
cat("Flower variable features: ", length(VariableFeatures(flower_new)), "\n")

# データリストを作成
object_list <- list(silique_new, flower_new)


# silique_newに接頭辞 "Silique_" を追加
silique_new <- RenameCells(silique_new, add.cell.id = "Silique")

# flower_newに接頭辞 "Flower_" を追加
flower_new <- RenameCells(flower_new, add.cell.id = "Flower")



# オブジェクトの細胞名を取得
silique_cells <- Cells(silique_new)
flower_cells <- Cells(flower_new)

# 重複している細胞名を確認
duplicate_cells <- intersect(silique_cells, flower_cells)
print(duplicate_cells)

# データリストを作成
object_list <- list(silique_new, flower_new)


dim(silique_new)  # silique_newの行列の次元数を確認
dim(flower_new)   # flower_newの行列の次元数を確認

# object_listの中身を確認
object_list

library(future)
plan(multisession, workers = 4)
options(future.globals.maxSize = 8 * 1024^3)

# アンカー検索の実行
anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:20, verbose = TRUE)

# 結果を確認
print(anchors)

# future.globals.maxSizeを50GBに設定
options(future.globals.maxSize = 50 * 1024^3)

# データ統合
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20)

print(integrated_data)

# スケーリング
integrated_data <- ScaleData(integrated_data, verbose = TRUE)
# 次元削減（PCA)
integrated_data <- RunPCA(integrated_data, verbose = TRUE)
# 結果確認
print(integrated_data[["pca"]])  # PCAの詳細を確認
ElbowPlot(integrated_data)      # 主成分の寄与を可視化
# 可視化（UMAPまたはt-SNE）
integrated_data <- RunUMAP(integrated_data, dims = 1:20)
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident")  # 元データごとに可視化
# クラスタリング
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
DimPlot(integrated_data, reduction = "umap", label = TRUE, group.by = "seurat_clusters")  # クラスターを可視化
saveRDS(integrated_data, file = "data/processed/integrated_data.rds")

# 統合オブジェクトが作成されたか確認
seurat_integrated <- IntegrateData(anchorset = anchors)  # anchors で統合

# 次元削減のリストを確認
Reductions(seurat_integrated)


# 正規化を実行
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)




# アクティブなアッセイをRNAに設定
DefaultAssay(seurat_integrated) <- "RNA"

# スケーリングの実行
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)






# データの構造を確認
str(seurat_integrated@meta.data)

# 行名の確認（細胞名）
head(colnames(seurat_integrated))

# 遺伝子発現データの確認
head(seurat_integrated[["integrated"]]@data)



# 明示的に遺伝子を指定してスケーリング
seurat_integrated <- ScaleData(
  seurat_integrated,
  features = rownames(seurat_integrated),
  verbose = FALSE
)

# スケールデータの確認
head(seurat_integrated[["integrated"]]@scale.data)





# PCAを実行
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = FALSE)

# PCA結果の確認
ElbowPlot(seurat_integrated)  # 主成分の重要性を確認

# UMAPを実行
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30, verbose = FALSE)

# 次元削減の確認
Reductions(seurat_integrated)

# 元データごとに可視化
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident")  
# DimPlotの凡例を消す
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident") +
  NoLegend()  # 凡例を非表示

# クラスターを可視化
DimPlot(integrated_data, reduction = "umap", label = TRUE, group.by = "seurat_clusters")  +
  NoLegend()  # 凡例を非表示



# 凝ったFeaturePlotを作成 1
FeaturePlot(
  seurat_integrated, 
  features = c("AT1G69180", "AT4G31820", "AT2G39060"), 
  ncol = 3,  
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)



# 凝ったFeaturePlotを作成 1
FeaturePlot(
  seurat_integrated, 
  features = c("AT3G01530", "AT3G25810", "AT1G65970"), 
  ncol = 3,  
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)



# 凝ったFeaturePlotを作成 1
FeaturePlot(
  seurat_integrated, 
  features = c("AT1G69180", "AT4G31820", "AT2G39060"), 
  ncol = 3,  
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)



# 凝ったFeaturePlotを作成 2
FeaturePlot(
  seurat_integrated, 
  features = c("AT3G01530", "AT3G25810", "AT1G62480"), 
  ncol = 3, 
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)

# 凝ったFeaturePlotを作成 3
FeaturePlot(
  seurat_integrated,
  features = c("AT5G20830", "AT2G26580", "AT1G77110"), 
  ncol = 3, 
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)

# 凝ったFeaturePlotを作成 4
FeaturePlot(
  seurat_integrated, 
  features = c("AT2G36190", "AT1G65970", "AT1G62480"), 
  ncol = 3, 
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)






#ここからヒートマップ

# 現在のメタデータ列名を確認
colnames(seurat_integrated@meta.data)

# 近傍探索を実行
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30)

# クラスター解析を実行
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)

#メタデータの確認
colnames(seurat_integrated@meta.data)
head(seurat_integrated@meta.data$seurat_clusters)


#レイヤーの確認
Layers(seurat_integrated[["RNA"]])

# スケールデータを取得
scale_data <- seurat_integrated[["RNA"]]@layers[["scale.data"]]

# データを確認
head(scale_data)

DoHeatmap(
  seurat_integrated, 
  features = c("AT1G69180", "AT4G31820", "AT2G39060", 
               "AT3G01530", "AT3G27810", "AT5G40350", 
               "AT5G20830", "AT2G26580", "AT1G77110", 
               "AT2G36190"), 
  group.by = "seurat_clusters", 
  slot = "scale.data"
)





# クラスターごとのヒートマップを作成
library(ggplot2)
DoHeatmap(
  seurat_integrated, 
  features = c("AT3G62230", "AT4G11543", "AT2G45403", "AT1G75870", "AT1G34095"), # 遺伝子リスト
  group.by = "seurat_clusters",  # クラスターごとに表示
  slot = "scale.data"            # スケール済みデータを使用
) +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0
  ) +
  theme(
    axis.text.y = element_text(face = "italic") # Y軸ラベルを斜体に設定
  )


# 現在のアクティブアッセイを確認
DefaultAssay(seurat_integrated)

# 必要であればRNAアッセイに切り替え
DefaultAssay(seurat_integrated) <- "RNA"



# 選んだ遺伝子をRNAアッセイで

all(c("AT1G69180", "AT4G31820", "AT2G39060", 
      "AT3G01530", "AT3G27810", "AT5G40350", 
      "AT5G20830", "AT2G26580", "AT1G77110", 
      "AT2G36190") %in% rownames(seurat_integrated))

head(seurat_integrated[["RNA"]]@scale.data)


DoHeatmap(
  seurat_integrated, 
  features = c("AT1G69180", "AT4G31820", "AT2G39060", 
               "AT3G01530", "AT3G27810", "AT5G40350", 
               "AT5G20830", "AT2G26580", "AT1G77110", 
               "AT2G36190"), 
  group.by = "seurat_clusters",  # クラスターごとに表示
  slot = "scale.data"
)








# 各クラスターで可変遺伝子から上位100を選択してヒートマップ

colnames(seurat_integrated@meta.data)

# 利用可能なグラフを確認
Reductions(seurat_integrated)
# グラフ名を指定して作成
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30, graph.name = "custom_snn", verbose = FALSE)

# 指定したグラフ名を使用してクラスター計算
seurat_integrated <- FindClusters(seurat_integrated, graph.name = "custom_snn", resolution = 0.5, verbose = FALSE)

head(seurat_integrated@meta.data$seurat_clusters)

# クラスター数の確認
table(seurat_integrated@meta.data$seurat_clusters)

# クラスターごとの細胞数
summary(seurat_integrated@meta.data$seurat_clusters)

#  ヒートマップの作成
DoHeatmap(
  seurat_integrated, 
  features = c("AT1G69180", "AT4G31820", "AT2G39060", 
               "AT3G01530", "AT3G27810", "AT5G40350", 
               "AT5G20830", "AT2G26580", "AT1G77110", 
               "AT2G36190"), 
  group.by = "seurat_clusters", 
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(face = "italic"))  # 遺伝子名を斜体に

# クラスター0の細胞
cluster_0_cells <- WhichCells(seurat_integrated, idents = 0)

# クラスター0の細胞数を確認
length(cluster_0_cells)

VlnPlot(seurat_integrated, features = c("AT1G62480", "AT3G25810"), group.by = "seurat_clusters")









# 各クラスターで可変遺伝子から上位100を選択
# 例: 最初のクラスターのデータを取得して確認
# 例: 最初のクラスターに対して処理

# デバッグ用に1つのクラスターで確認
cluster <- unique(seurat_integrated$seurat_clusters)[1]  # 最初のクラスター
cells_in_cluster <- WhichCells(seurat_integrated, idents = cluster)  # 細胞を取得
subset_obj <- subset(seurat_integrated, cells = cells_in_cluster)  # サブセット作成

# スケールデータが存在しなければ作成
if (!"scale.data" %in% Layers(subset_obj[["RNA"]])) {
  subset_obj <- ScaleData(subset_obj, verbose = FALSE)
}

# 可変遺伝子の取得
subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 2000)
top_genes <- VariableFeatures(subset_obj)[1:min(100, length(VariableFeatures(subset_obj)))]

print(top_genes)  # 可変遺伝子を確認


# クラスターごとに可変遺伝子を取得
top_genes_per_cluster <- lapply(unique(seurat_integrated$seurat_clusters), function(cluster) {
  # クラスターごとの細胞を抽出
  cells_in_cluster <- WhichCells(seurat_integrated, idents = cluster)
  
  # クラスター内の細胞だけでサブセット
  subset_obj <- subset(seurat_integrated, cells = cells_in_cluster)
  
  # スケールデータが存在しない場合は計算
  if (!"scale.data" %in% Layers(subset_obj[["RNA"]])) {
    subset_obj <- ScaleData(subset_obj, verbose = FALSE)
  }
  
  # 可変遺伝子を取得
  subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 2000)
  top_genes <- VariableFeatures(subset_obj)
  
  # 上位100遺伝子を返す（不足している場合はすべて返す）
  return(top_genes[1:min(100, length(top_genes))])
})

# エラーの原因となる部分をデバッグ
print(lengths(top_genes_per_cluster))  # 各クラスターの遺伝子数を確認










library(ggplot2)

# ヒートマップを作成
DoHeatmap(
  seurat_integrated, 
  features = unique_genes, # 各クラスターで選ばれた遺伝子を使用
  group.by = "seurat_clusters",  # クラスターごとに表示
  slot = "scale.data"
) +
  labs(title = "Heatmap of top 100 variable genes per cluster") +
  theme(
    axis.text.y = element_text(face = "italic"), # 遺伝子名を斜体に
    plot.title = element_text(hjust = 0.5)      # タイトルを中央揃え
  )




# DotPlotを作成
# 現在のアクティブなアッセイを確認
DefaultAssay(seurat_integrated)

# 必要ならRNAアッセイに変更
DefaultAssay(seurat_integrated) <- "RNA"

# 遺伝子名をマッピング
gene_mapping <- c(
  "AT1G69180" = "CRC",
  "AT4G31820" = "MAB4",
  "AT2G39060" = "SWEET9",
  "AT3G01530" = "MYB57",
  "AT3G27810" = "MYB21",
  "AT5G40350" = "MYB24",
  "AT5G20830" = "SUS1",
  "AT2G26580" = "YAB5",
  "AT1G77110" = "PIN6",
  "AT2G36190" = "cwINV4"
)

# DotPlotの作成
DotPlot(
  seurat_integrated, 
  features = names(gene_mapping),  # 遺伝子IDを指定
  group.by = "orig.ident",         # flowerとsiliqueでグループ化
  dot.scale = 12                    # 円の大きさを設定
) +
  scale_color_gradientn(colors = c("grey95", "blue4")) +  # 青をblue4に設定
  scale_y_discrete(labels = gene_mapping) +              # 遺伝子名を表示
  labs(title = "DotPlot of Gene Expression in Flower and Silique") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # X軸ラベルを45度傾ける
    axis.text.y = element_text(face = "italic"),        # Y軸ラベルを斜体に
    plot.title = element_text(hjust = 0.5)             # タイトルを中央揃え
  )

# DotPlotオブジェクトを保存
p <- DotPlot(
  seurat_integrated, 
  features = names(gene_mapping),  # 遺伝子IDを指定
  group.by = "orig.ident",         # flowerとsiliqueでグループ化
  dot.scale = 12                   # 円の大きさを設定
)
# データのY軸ラベルをマッピングに置き換え
p$data$features.plot <- gene_mapping[as.character(p$data$features.plot)]
# プロットを描画
p +
  scale_color_gradientn(colors = c("grey95", "blue4")) +  # 青をblue4に設定
  labs(title = "DotPlot of Gene Expression in Flower and Silique") + 
  theme(
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),  # X軸ラベルを斜体に
    axis.text.y = element_text(face = "plain"),                         # Y軸ラベルを通常に
    plot.title = element_text(hjust = 0.5)                              # タイトルを中央揃え
  )

# スタックされたバイオリンプロットを作成
VlnPlot(
  object = seurat_integrated, 
  features = c("AT3G01530", "AT3G27810", "AT5G40350"),  # 表示する遺伝子
  group.by = "orig.ident",                             # flowerとsiliqueでグループ化
  pt.size = 0,                                         # ポイントサイズ（必要に応じて調整）
  cols = c("blue", "red")                              # カラースケールを指定
) +
  labs(title = "Stacked Violin Plot of Gene Expression in Flower and Silique") +
  theme(
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),  # X軸ラベルを斜体に
    plot.title = element_text(hjust = 0.5)                              # タイトルを中央揃え
  )
# Seuratオブジェクトを保存
saveRDS(seurat_integrated, file = "seurat_integrated.rds")

# メタデータの保存
write.csv(seurat_integrated@meta.data, file = "metadata.csv")

