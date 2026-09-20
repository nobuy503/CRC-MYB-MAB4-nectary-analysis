# Source script for silique snRNA-seq visualizations in Supplementary Figure 1
#
# Dataset: silique portion of NCBI GEO GSE226097.
#
# This file preserves the analysis code supplied by the authors. The original
# local input path was replaced with the repository-relative path
# data/GSE226097_silique_230221.rds, and full-width spaces that prevent R
# parsing were converted to standard spaces. The analysis commands themselves
# were not changed.
#
# This source file includes exploratory plotting commands in addition to code
# used for the figure. A streamlined reproducible version will be prepared
# after the final panel(s) and object-preparation steps are confirmed.

#パッケージのインストール
install.packages("Seurat")
install.packages("dplyr")
install.packages("patchwork")


#必要なライブラリーの準備
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(viridis)

#UMAPを描くだけ
GSE226097_silique_230221 <- readRDS("data/GSE226097_silique_230221.rds")
class(GSE226097_silique_230221) 
p1 <- DimPlot(GSE226097_silique_230221)
print(p1)

#クラスターの順にUMAPを描くだけ
# オブジェクトのクラスを確認
object_class <- class(GSE226097_silique_230221)
cat("Class of the object:", object_class, "\n")
# クラスターの順序を確認
current_levels <- levels(Idents(GSE226097_silique_230221))
cat("Current cluster levels:", current_levels, "\n")
# 正しい順序に並べ替える（例: 数値順）
desired_order <- sort(as.numeric(as.character(current_levels)))
Idents(GSE226097_silique_230221) <- factor(Idents(GSE226097_silique_230221), levels = as.character(desired_order))
# DimPlotを作成して表示
p1 <- DimPlot(GSE226097_silique_230221)
print(p1)



#普通にFeaturePlotを描く
DefaultAssay(GSE226097_silique_230221) <- "RNA" #RNAで指定が必要
FeaturePlot(GSE226097_silique_230221, features = c("AT1G69180", "AT4G31820", "AT2G39060"), ncol = 3) 

#陽性細胞を全面に出してにFeaturePlotを描く
DefaultAssay(GSE226097_silique_230221) <- "RNA" #RNAで指定が必要
FeaturePlot(GSE226097_silique_230221, features = c("AT1G69180", "AT4G31820", "AT2G39060"), ncol = 3, order =T) 

# 凝ったFeaturePlotを作成 1
FeaturePlot(
  GSE226097_silique_230221, 
  features = c("AT1G69180", "AT4G31820", "AT2G39060"), 
  ncol = 3,  
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)

# 凝ったFeaturePlotを作成 2
FeaturePlot(
  GSE226097_silique_230221, 
  features = c("AT3G01530", "AT3G27810", "AT5G40350"), 
  ncol = 3, 
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)

# 凝ったFeaturePlotを作成 3
FeaturePlot(
  GSE226097_silique_230221, 
  features = c("AT5G20830", "AT2G26580", "AT1G77110"), 
  ncol = 3, 
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)

# 凝ったFeaturePlotを作成 4
FeaturePlot(
  GSE226097_silique_230221, 
  features = c("AT2G36190", "AT1G65970", "AT1G62480"), 
  ncol = 3, 
  order = TRUE, 
  cols = c("grey95", "blue4"),  # カラースケールを濃い青に変更
  min.cutoff = 0,  # カットオフの最小値を指定
  max.cutoff = 2  # カットオフの最大値を指定（必要に応じて変更）
)








#普通にdot plotを描く
DefaultAssay(GSE226097_silique_230221) <- "RNA" #RNAで指定が必要
DotPlot(GSE226097_silique_230221, features = c("AT1G69180", "AT4G31820", "AT2G39060", "AT3G01530", "AT3G27810", "AT5G40350", "AT5G20830", "AT2G26580", "AT1G77110", "AT2G36190"), dot.scale = 20) 

#dot plotを描く、gene mappingとラベルを斜めに、Y軸ラベルを直接編集
DefaultAssay(GSE226097_silique_230221) <- "RNA" #RNAで指定が必要
# DotPlotを作成してオブジェクトとして保存
p <- DotPlot(GSE226097_silique_230221, 
             features = c("AT1G69180", "AT4G31820", "AT2G39060", 
                          "AT3G01530", "AT3G27810", "AT5G40350", 
                          "AT5G20830", "AT2G26580", "AT1G77110", 
                          "AT2G36190"), 
             dot.scale = 20)
# 現在のY軸ラベルを確認
p$data$features.plot
# gene_mappingを調整（例: 大文字・小文字やフォーマットを揃える）
adjusted_gene_mapping <- c(
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
# Y軸ラベルを直接編集
p$data$features.plot <- adjusted_gene_mapping[as.character(p$data$features.plot)]
# プロットを表示（X軸ラベルを斜体に設定）
p + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), # X軸ラベルを斜体に
    axis.text.y = element_text()                                      # Y軸はデフォルト
  )
# プロットを表示（X軸ラベルを斜体にし、Viridisカラーマップを適用）
p + 
  scale_color_viridis() +  # Viridisカラーマップを追加
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), # X軸ラベルを斜体に
    axis.text.y = element_text()                                      # Y軸はデフォルト
  )


# 普通にVlnplotを表示
VlnPlot(GSE226097_silique_230221,
        features = c("AT1G69180", "AT4G31820", "AT2G39060", 
                     "AT3G01530", "AT3G27810", "AT5G40350", 
                     "AT5G20830", "AT2G26580", "AT1G77110", 
                     "AT2G36190"),
        stack = T) 

# 普通にRidgeplotを表示
RidgePlot(
  GSE226097_silique_230221,
  features = c("AT1G69180"),
  ncol = 2  # プロットを3列で配置
  )

# 普通にHeatmapを表示
# スケーリングを実行
GSE226097_silique_230221 <- ScaleData(GSE226097_silique_230221, features = rownames(GSE226097_silique_230221))
DoHeatmap(GSE226097_silique_230221, features = c("AT1G69180", "AT4G31820", "AT2G39060", 
                                                "AT3G01530", "AT3G27810", "AT5G40350", 
                                                "AT5G20830", "AT2G26580", "AT1G77110", 
                                                "AT2G36190"))
# スロット値を表示, 表示域を変更 
DoHeatmap(GSE226097_silique_230221, features = c("AT1G69180", "AT4G31820", "AT2G39060", 
                                                "AT3G01530", "AT3G27810", "AT5G40350", 
                                                "AT5G20830", "AT2G26580", "AT1G77110", 
                                                "AT2G36190"), slot = "data", disp.max = 1.5)
# 遺伝子IDと遺伝子名の対応表を作成
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

# DoHeatmapを作成
p <- DoHeatmap(
  GSE226097_silique_230221, 
  features = names(gene_mapping), # 遺伝子IDを指定
  slot = "data", 
  disp.max = 1.5
)

# Y軸ラベルを遺伝子名に置き換える
p + scale_y_discrete(labels = gene_mapping) + theme(
  axis.text.y = element_text(face = "italic") # Y軸ラベルを斜体に設定
)

# DoHeatmapを作成 色違い
# 遺伝子IDと遺伝子名の対応表を作成
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

# DoHeatmapを作成
p <- DoHeatmap(
  GSE226097_silique_230221, 
  features = names(gene_mapping), # 遺伝子IDを指定
  slot = "data", 
  disp.max = 1.5
)

# カスタムカラースケールを追加
p + 
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0.2
  ) +
  scale_y_discrete(labels = gene_mapping) + # Y軸ラベルを遺伝子名に変更
  theme(
    axis.text.y = element_text(face = "italic") # Y軸ラベルを斜体に設定
  )


# 遺伝子を増やしたヒートマップ
# 遺伝子IDと遺伝子名の対応表を作成
gene_mapping <- c(
  "AT1G69180" = "CRC", "AT4G31820" = "MAB4", "AT2G39060" = "SWEET9", 
  "AT3G01530" = "MYB57", "AT3G27810" = "MYB21", "AT5G40350" = "MYB24", 
  "AT5G20830" = "SUS1", "AT2G26580" = "YAB5", "AT1G77110" = "PIN6", 
  "AT2G36190" = "cwINV4"
)

# 指定された遺伝子リスト
gene_list <- c(
  "AT4G36220", "AT1G51680", "AT2G26750", "AT2G27770", "AT1G01490",
  "AT2G26580", "AT5G22500", "AT5G54160", "AT4G31820", "AT2G36090",
  "AT5G24270", "AT1G32540", "AT2G39060", "AT1G18720", "AT4G12530",
  "AT1G23010", "AT4G00950", "AT1G23300", "AT3G10040", "AT5G06720",
  "AT1G51340", "AT5G20830", "AT1G43800", "AT5G65110", "AT5G11090",
  "AT4G29010", "AT1G23200", "AT1G77110", "AT3G60780", "AT1G12805",
  "AT1G65970", "AT1G19640", "AT5G40350", "AT3G21270", "AT5G24580",
  "AT3G10340", "AT1G09155", "AT2G35920", "AT2G38940", "AT1G03620",
  "AT4G04020", "AT5G39890", "AT1G62480", "AT1G68825", "AT2G32510",
  "AT5G26310", "AT1G77590", "AT1G70680", "AT3G25640", "AT2G28840",
  "AT1G78780", "AT2G28110", "AT3G05640", "AT2G22800", "AT4G24120",
  "AT3G11430", "AT4G34135", "AT1G17960", "AT1G70670", "AT5G15120",
  "AT5G24150", "AT3G14990"
)

# 重複をチェックし、統一リストを作成
unique_genes <- unique(c(names(gene_mapping), gene_list))

# 統一リストの確認
cat("統一された遺伝子リスト:\n", unique_genes, "\n")

# DoHeatmapを作成
p <- DoHeatmap(
  GSE226097_silique_230221,
  features = unique_genes, # 統一リストを使用
  slot = "data",           # 表示するスロット
  disp.max = 1.5           # カラースケールの最大値
)

# 遺伝子名マッピングを適用してヒートマップを表示
p +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0.2
  ) +
  scale_y_discrete(labels = function(x) ifelse(!is.na(gene_mapping[x]), gene_mapping[x], x)) + # マッピングがない場合AGIコードのまま表示
  theme(
    axis.text.y = element_text(face = "italic") # Y軸ラベルを斜体に設定
  )

# 遺伝子名マッピングを適用してヒートマップを表示 ver2
p +
  scale_y_discrete(labels = function(x) ifelse(!is.na(gene_mapping[x]), gene_mapping[x], x)) + # マッピングがない場合AGIコードのまま表示
  theme(
    axis.text.y = element_text(face = "italic") # Y軸ラベルを斜体に設定
  )
