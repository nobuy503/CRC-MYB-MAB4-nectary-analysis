# GO enrichment dot plots for integrated snRNA-seq clusters
#
# Input:
#   data/processed/GO_for_R.xlsx
#
# This script visualizes previously generated GO-enrichment results for
# clusters 0, 15 and 19. It does not perform the enrichment test itself.
# The enrichment results were obtained separately using agriGO v2.0
# Singular Enrichment Analysis (SEA).
#
# Figure relationship:
#   The cluster 15 plot is associated with Fig. 1o.
#
# Statistical-display note:
#   This source script maps the raw column named "p-value" to the x-axis
#   and colour scale. If the final figure/legend reports FDR-adjusted
#   P values, the corresponding FDR column from GO_for_R.xlsx must be
#   selected explicitly in the final plotting code.

## ================================
##  GO dotplot for clusters (0,15,19)
##  - X軸: -log10(p-value)（上限 20）
##  - Y軸: group + Description（p=0 でもラベルは残す）
##  - 色: -log10(p-value)（viridis）
##  - 丸サイズ: number_in_list
## ================================

library(readxl)
library(dplyr)
library(forcats)
library(ggplot2)

## データ読み込み -----------------------
input_dir <- "data/processed"
output_dir <- "results/figure1_go"
go_path  <- file.path(input_dir, "GO_for_R.xlsx")

go_all <- read_excel(go_path)

## プロット用関数 -----------------------
plot_go_cluster <- function(df, cl){
  
  # そのクラスタの全行（p=0も含む）
  df_base <- df %>%
    filter(cluster == cl) %>%
    mutate(
      # p>0 のときだけ -log10(p) を計算（p=0 は NA にしておく）
      neglog10_p_raw = if_else(`p-value` > 0,
                               -log10(`p-value`),
                               NA_real_),
      
      # group の順番を固定したければここを調整
      group = factor(
        group,
        levels = c("development", "pollen", "glucose", "stomata", "response")
      ),
      
      # Y軸ラベル用：group + Description
      label_raw = paste(group, Description, sep = " : ")
    ) %>%
    # 並び順：group → 有意度（NAは最後）
    arrange(group, dplyr::desc(neglog10_p_raw)) %>%
    mutate(
      # 並べた順をそのままy軸に反映（上から下）
      label = factor(label_raw, levels = rev(unique(label_raw)))
    )
  
  # 実際に点を描くのは p>0 の行だけ
  df_points <- df_base %>%
    filter(!is.na(neglog10_p_raw)) %>%
    mutate(
      # プロット用に上限20でクリップ
      neglog10_p = pmin(neglog10_p_raw, 20)
    )
  
  ggplot() +
    geom_point(
      data  = df_points,
      aes(x = neglog10_p,
          y = label,
          size  = number_in_list,
          color = neglog10_p)
    ) +
    # 色スケール（viridis、0〜20）
    scale_color_viridis_c(
      limits = c(0, 20),
      name   = expression(-log[10]("p-value"))
    ) +
    # サイズスケール
    scale_size(
      name  = "Number in list",
      range = c(2, 8)
    ) +
    # Y軸は df_base のラベル全部を使う（p=0だけの行も残す）
    scale_y_discrete(
      limits = levels(df_base$label),
      drop   = FALSE
    ) +
    labs(
      title = paste("GO enrichment – cluster", cl),
      x     = expression(-log[10]("p-value")),
      y     = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid       = element_blank(),  # 灰色の補助線を消す
      axis.text.y      = element_text(size = 8)
    )
}

## クラスタごとの図を作成 ----------------
# クラスタ番号は必要に応じて変更
p0  <- plot_go_cluster(go_all, 0)
p15 <- plot_go_cluster(go_all, 15)
p19 <- plot_go_cluster(go_all, 19)

p0
p15
p19

## 保存したい場合 ------------------------
ggsave(file.path(output_dir, "GO_cluster0.pdf"),  p0,  width = 7, height = 6)
ggsave(file.path(output_dir, "GO_cluster15.pdf"), p15, width = 7, height = 6)
ggsave(file.path(output_dir, "GO_cluster19.pdf"), p19, width = 7, height = 6)
