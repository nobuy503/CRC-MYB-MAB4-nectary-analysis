# Flower-only pseudotime analysis for Fig. 1p-t
#
# Input:
#   data/processed/GSE226097_flower_230221.rds
#
# Main outputs represented in this source script:
#   Fig. 1q: Monocle3 pseudotime trajectory
#   Fig. 1r-s: CRC and AT1G65970 expression projected on the trajectory
#   Fig. 1t: GAM-smoothed expression of six genes along pseudotime
#
# Important provenance notes:
# - This is the analysis source used during figure preparation and contains
#   exploratory plotting blocks in addition to the final plotting code.
# - The script calculates the principal-graph vertex most frequently associated
#   with CRC-expressing nuclei, but then explicitly sets root_node to "Y_111".
# - The five Seurat clusters described for Fig. 1p are not explicitly selected
#   in this script; cluster/cell-type annotations are inherited from the input
#   Seurat object's metadata.
# - Pseudotime represents a putative ordering of transcriptional states and
#   should not be interpreted as a directly observed developmental timeline.

flower <- readRDS("data/processed/GSE226097_flower_230221.rds")
class(flower)
dim(flower)

Assays(flower)
colnames(flower@meta.data)
head(flower@meta.data)

DefaultAssay(flower) <- "RNA"
Layers(flower[["RNA"]])

counts_flower <- GetAssayData(
  flower,
  assay = "RNA",
  layer = "counts"
)

dim(counts_flower)
class(counts_flower)

cell_metadata <- flower@meta.data

gene_metadata <- data.frame(
  gene_short_name = rownames(counts_flower),
  row.names = rownames(counts_flower)
)

identical(colnames(counts_flower), rownames(cell_metadata))


library(monocle3)

cds_flower <- new_cell_data_set(
  expression_data = counts_flower,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds_flower

cds_flower <- preprocess_cds(
  cds_flower,
  num_dim = 30
)

cds_flower

cds_flower <- reduce_dimension(
  cds_flower,
  reduction_method = "UMAP"
)

cds_flower

cds_flower <- cluster_cells(
  cds_flower,
  reduction_method = "UMAP"
)

plot_cells(
  cds_flower,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  show_trajectory_graph = FALSE
)

cds_flower <- learn_graph(cds_flower)

plot_cells(
  cds_flower,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  label_branch_points = TRUE,
  label_leaves = TRUE,
  show_trajectory_graph = TRUE
)

plot_cells(
  cds_flower,
  color_cells_by = "CellType",
  label_cell_groups = TRUE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  label_branch_points = FALSE,
  label_leaves = FALSE
)

genes <- c(
  CRC   = "AT1G69180",
  MYB57 = "AT3G01530",
  MYB21 = "AT3G27810",
  MYB24 = "AT5G40350",
  MAB4  = "AT4G31820"
)

genes %in% rownames(cds_flower)

plot_cells(
  cds_flower,
  genes = unname(genes),
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  min_expr = 0.5
)

crc_counts <- counts(cds_flower)["AT1G69180", ]

crc_cells <- colnames(cds_flower)[crc_counts > 0]

length(crc_cells)


closest_vertex <- cds_flower@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

crc_vertex <- closest_vertex[crc_cells, 1]

root_node <- names(sort(table(crc_vertex), decreasing = TRUE))[1]

root_node

root_node <- "Y_111"

cds_flower <- order_cells(
  cds_flower,
  root_pr_nodes = root_node
)

plot_cells(
  cds_flower,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = TRUE,
  label_leaves = TRUE
)

pt <- pseudotime(cds_flower)

expr <- log1p(as.matrix(
  counts(cds_flower)[
    c("AT1G69180", "AT3G01530", "AT3G27810", "AT5G40350", "AT4G31820"),
    names(pt)
  ]
))

plot_df <- data.frame(
  pseudotime = pt,
  CRC   = expr["AT1G69180", ],
  MYB57 = expr["AT3G01530", ],
  MYB21 = expr["AT3G27810", ],
  MYB24 = expr["AT5G40350", ],
  MAB4  = expr["AT4G31820", ]
)

head(plot_df)

plot_df2 <- plot_df[is.finite(plot_df$pseudotime), ]

nrow(plot_df2)
range(plot_df2$pseudotime)


library(dplyr)
library(tidyr)
library(ggplot2)

plot_df2$pt_bin <- cut(
  plot_df2$pseudotime,
  breaks = 10,
  include.lowest = TRUE,
  labels = FALSE
)

mean_df <- plot_df2 %>%
  group_by(pt_bin) %>%
  summarise(
    pseudotime = mean(pseudotime),
    CRC   = mean(CRC),
    MYB57 = mean(MYB57),
    MYB21 = mean(MYB21),
    MYB24 = mean(MYB24),
    MAB4  = mean(MAB4),
    .groups = "drop"
  )

mean_df


mean_long <- mean_df %>%
  pivot_longer(
    cols = c(CRC, MYB57, MYB21, MYB24, MAB4),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(gene) %>%
  mutate(
    scaled_expression = as.numeric(scale(expression))
  ) %>%
  ungroup()

ggplot(
  mean_long,
  aes(
    x = pseudotime,
    y = scaled_expression,
    group = gene,
    linetype = gene
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Pseudotime",
    y = "Scaled mean expression"
  ) +
  theme_classic()


library(ggplot2)
library(tidyr)
library(dplyr)

plot_long <- plot_df2 %>%
  pivot_longer(
    cols = c(CRC, MYB57, MYB21, MYB24, MAB4),
    names_to = "gene",
    values_to = "expression"
  )

ggplot(
  plot_long,
  aes(
    x = pseudotime,
    y = expression
  )
) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_classic()

crc_info <- data.frame(
  cluster = clusters(cds_flower)[crc_cells],
  CellType = colData(cds_flower)[crc_cells, "CellType"]
)

table(crc_info$cluster)
table(crc_info$CellType)

colData(cds_flower)$CRC_positive <- ifelse(
  colnames(cds_flower) %in% crc_cells,
  "CRC-positive",
  "Other"
)

plot_cells(
  cds_flower,
  color_cells_by = "CRC_positive",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE
)

plot_cells(
  cds_flower,
  color_cells_by = "CRC_positive",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.4
) +
  scale_color_manual(
    values = c(
      "Other" = "grey90",
      "CRC-positive" = "red"
    )
  )

flower$CRC_positive <- ifelse(
  colnames(flower) %in% crc_cells,
  "CRC-positive",
  "Other"
)

library(Seurat)
DimPlot(
  flower,
  reduction = "umap",
  group.by = "CRC_positive",
  order = c("Other", "CRC-positive"),
  pt.size = 0.2
) +
  scale_color_manual(
    values = c(
      "Other" = "grey90",
      "CRC-positive" = "red"
    )
  )


FeaturePlot(
  flower,
  features = "AT1G69180",
  reduction = "umap",
  order = TRUE,
  pt.size = 0.3
)

plot_cells(
  cds_flower,
  genes = "AT1G69180",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

plot_cells(
  cds_flower,
  genes = "AT2G39060",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

plot_cells(
  cds_flower,
  genes = "AT1G65970",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

plot_cells(
  cds_flower,
  genes = "AT1G62480",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)


plot_cells(
  cds_flower,
  genes = "AT3G01530",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

plot_cells(
  cds_flower,
  genes = "AT3G27810",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

plot_cells(
  cds_flower,
  genes = "AT5G40350",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

plot_cells(
  cds_flower,
  genes = "AT4G31820",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
)

p_pt <- plot_cells(
  cds_flower,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.5,
  trajectory_graph_segment_size = 0.7
) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 8, 5, 5)
  )

p_pt

p_pt <- p_pt +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

p_pt


p_type <- plot_cells(
  cds_flower,
  color_cells_by = "CellType",
  show_trajectory_graph = TRUE,
  label_cell_groups = TRUE,
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.5,
  trajectory_graph_segment_size = 0.7
) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

p_type


p_type <- plot_cells(
  cds_flower,
  color_cells_by = "CellType",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.5,
  trajectory_graph_segment_size = 0.7
) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(),
    axis.line = element_line(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

p_type


p_pt <- plot_cells(
  cds_flower,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.5,
  trajectory_graph_segment_size = 0.7
) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(),
    axis.line = element_line(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

p_pt




crc_expr <- log1p(as.numeric(
  SingleCellExperiment::counts(cds_flower)["AT1G69180", ]
))

SummarizedExperiment::colData(cds_flower)$CRC_percent_max <-
  100 * crc_expr / max(crc_expr)

"CRC_percent_max" %in% colnames(
  SummarizedExperiment::colData(cds_flower)
)


p_crc <- plot_cells(
  cds_flower,
  genes = "AT1G69180",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  cell_size = 0.8,
  min_expr = 0.1
) +
  scale_color_gradientn(
    colours = viridisLite::viridis(3),
    values = c(0, 0.25, 1),
    name = "% Max"
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 12)

p_crc



library(dplyr)
library(tidyr)
library(ggplot2)

plot_long <- plot_df2 %>%
  pivot_longer(
    cols = c(CRC, MYB57, MYB21, MYB24, MAB4),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(gene) %>%
  mutate(
    expression_scaled = expression / max(expression, na.rm = TRUE)
  ) %>%
  ungroup()


p_gene_pt <- ggplot(
  plot_long,
  aes(
    x = pseudotime,
    y = expression_scaled,
    color = gene
  )
) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    linewidth = 1.2
  ) +
  scale_color_viridis_d(
    option = "viridis",
    end = 0.9
  ) +
  labs(
    x = "Pseudotime",
    y = "Scaled expression",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

p_gene_pt




plot_long <- plot_df2 %>%
  pivot_longer(
    cols = c(CRC, MYB57, MYB21, MYB24, MAB4),
    names_to = "gene",
    values_to = "expression"
  )

gene_cols <- setNames(
  viridisLite::viridis(5, option = "viridis", end = 0.9),
  c("CRC", "MYB57", "MYB21", "MYB24", "MAB4")
)

p_gene_pt <- ggplot(
  plot_long,
  aes(
    x = pseudotime,
    y = expression,
    color = gene
  )
) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    linewidth = 1.1
  ) +
  facet_wrap(
    ~ gene,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(values = gene_cols) +
  labs(
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic",
      size = 11,
      hjust = 0
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.spacing = grid::unit(0.35, "lines")
  )

p_gene_pt



genes_use <- c(
  "AT1G69180",  # CRC
  "AT4G31820",  # MAB4
  "AT2G39060",  # SWEET9
  "AT3G01530",  # MYB57
  "AT3G27810",  # MYB21
  "AT1G65970"
)

gene_names <- c(
  "CRC",
  "MAB4",
  "SWEET9",
  "MYB57",
  "MYB21",
  "AT1G65970"
)

expr6 <- log1p(as.matrix(
  counts(cds_flower)[genes_use, names(pt)]
))

plot_df6 <- data.frame(
  pseudotime = pt,
  t(expr6)
)

colnames(plot_df6)[-1] <- gene_names

plot_df6 <- plot_df6[
  is.finite(plot_df6$pseudotime),
]

plot_long6 <- plot_df6 %>%
  pivot_longer(
    cols = all_of(gene_names),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = c(
        "CRC",
        "MAB4",
        "SWEET9",
        "MYB57",
        "MYB21",
        "AT1G65970"
      )
    )
  )


gene_cols <- setNames(
  viridisLite::viridis(
    6,
    option = "viridis",
    end = 0.9
  ),
  levels(plot_long6$gene)
)

p_gene_pt6 <- ggplot(
  plot_long6,
  aes(
    x = pseudotime,
    y = expression,
    color = gene
  )
) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    linewidth = 1.1
  ) +
  facet_wrap(
    ~ gene,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(values = gene_cols) +
  labs(
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic",
      size = 11,
      hjust = 0
    ),
    panel.spacing = grid::unit(0.35, "lines")
  )

p_gene_pt6
