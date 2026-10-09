# -------------------------------
# Step 0: 参数设置
# 说明：所有可修改参数集中放在开头，便于手动调整和复现
# -------------------------------
input_rdata <- "spa .rdata"
spa_object_name <- "spa"

output_root <- file.path("fig", "Tcell_gene_gradient_ring_pooled_95CI")
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_root, paste0("run_", run_id))

sample_indices <- NULL
spatial_assay <- "Spatial"
expression_layer <- "data"
coord_cols <- c("x", "y")
cell_type_col <- "type"
niche_col <- "nich4"

target_t_cells <- c("CTL", "Tex_T8", "TRM_T8")
ring_breaks <- c(0, 10, 30, 50, 70, 90, 110)
ring_lower <- head(ring_breaks, -1)
ring_upper <- tail(ring_breaks, -1)
ring_labels <- paste0(ring_lower, "-", ring_upper)
include_self <- FALSE

target_genes <- c(
  "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
  "Ccl11", "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21a",
  "Cxcl9", "Cxcl10", "Cxcl12", "Cxcl13", "Cxcl16",
  "Tgfb1", "Tnf", "Ifng", "Il2", "Il6", "Il7", "Il12", "Il17a"
)

# 只控制出图的基因；设为 NULL 时画 target_genes 中所有基因
plot_genes <- c("Ccl19", "Ccl21a", "Cxcl9", "Cxcl10", "Cxcl13", "Tgfb1")

plot_ncol <- 4
plot_width <- 14
plot_height_per_row <- 2.2

# -------------------------------
# Step 1: 加载依赖包
# 说明：缺包时直接停止，避免长任务运行到中途才报错
# -------------------------------
required_packages <- c("Seurat", "dplyr", "ggplot2", "Matrix")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "缺少以下 R 包：", paste(missing_packages, collapse = ", "),
    "。请先手动安装后再运行本脚本。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

# -------------------------------
# Step 2: 工具函数
# 说明：这些短函数用于基因匹配、稳健汇总和 95% 置信区间计算
# -------------------------------
standardize_gene_name <- function(x) {
  # 统一基因名格式，减少大小写或符号差异导致的匹配失败
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

safe_mean <- function(x) {
  # 没有有效值时返回 NA，避免 NaN 混入输出表
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_se <- function(x) {
  # 标准误用于目标细胞层面的 95% 置信区间
  valid_x <- x[!is.na(x)]
  if (length(valid_x) <= 1) {
    return(NA_real_)
  }
  sd(valid_x) / sqrt(length(valid_x))
}

safe_sd <- function(x) {
  # 标准差用于输出表记录细胞层面离散程度
  valid_x <- x[!is.na(x)]
  if (length(valid_x) <= 1) {
    return(NA_real_)
  }
  sd(valid_x)
}

scale_for_plot <- function(x) {
  # 热图按每个基因内部做 z-score，只用于展示相对趋势
  out <- rep(NA_real_, length(x))
  valid_idx <- which(!is.na(x))
  if (length(valid_idx) <= 1) {
    return(out)
  }
  x_sd <- sd(x[valid_idx])
  if (!is.finite(x_sd) || x_sd == 0) {
    out[valid_idx] <- 0
    return(out)
  }
  out[valid_idx] <- as.numeric(scale(x[valid_idx]))
  out
}

match_target_genes <- function(target_genes, feature_names) {
  # 将目标基因匹配到空间表达矩阵的 feature 名，并保留原始基因顺序
  feature_map <- data.frame(
    matched_feature = feature_names,
    feature_key = standardize_gene_name(feature_names),
    stringsAsFactors = FALSE
  )
  feature_map <- feature_map[!duplicated(feature_map$feature_key), , drop = FALSE]

  gene_key <- standardize_gene_name(target_genes)
  matched_idx <- match(gene_key, feature_map$feature_key)

  data.frame(
    requested_gene = target_genes,
    matched_feature = feature_map$matched_feature[matched_idx],
    matched = !is.na(matched_idx),
    stringsAsFactors = FALSE
  )
}

calculate_gene_ring_one_sample <- function(obj, sample_label) {
  # 对单个样本计算：每个目标 T 细胞为中心，每个环带内的基因平均表达
  if (!inherits(obj, "Seurat")) {
    stop(sample_label, " 不是 Seurat 对象。")
  }
  if (!spatial_assay %in% Assays(obj)) {
    stop(sample_label, " 中没有 assay = '", spatial_assay, "'。")
  }
  if (!expression_layer %in% Layers(obj[[spatial_assay]])) {
    stop(sample_label, " 的 assay = '", spatial_assay, "' 中没有 layer = '", expression_layer, "'。")
  }

  meta <- obj@meta.data[colnames(obj), , drop = FALSE]
  if (!all(coord_cols %in% colnames(meta))) {
    stop(sample_label, " 的 meta.data 中缺少坐标列：", paste(coord_cols, collapse = ", "))
  }
  if (!cell_type_col %in% colnames(meta)) {
    stop(sample_label, " 的 meta.data 中缺少细胞类型列：", cell_type_col)
  }

  gene_match <- match_target_genes(target_genes, rownames(obj[[spatial_assay]]))
  matched_gene_match <- gene_match %>%
    filter(matched)

  if (nrow(matched_gene_match) == 0) {
    stop(sample_label, " 中没有匹配到任何目标基因，请检查 assay 或基因名。")
  }
  if (any(!gene_match$matched)) {
    warning(
      sample_label, " 中以下基因没有匹配到，已跳过：",
      paste(gene_match$requested_gene[!gene_match$matched], collapse = ", ")
    )
  }

  expr_mat <- GetAssayData(obj, assay = spatial_assay, layer = expression_layer)[
    matched_gene_match$matched_feature,
    ,
    drop = FALSE
  ]
  rownames(expr_mat) <- matched_gene_match$requested_gene

  global_gene_mean <- Matrix::rowMeans(expr_mat)
  global_gene_mean[!is.finite(global_gene_mean) | global_gene_mean == 0] <- NA_real_

  spot_df <- data.frame(
    cell_id = colnames(obj),
    x = as.numeric(meta[[coord_cols[1]]]),
    y = as.numeric(meta[[coord_cols[2]]]),
    cell_type = as.character(meta[[cell_type_col]]),
    stringsAsFactors = FALSE
  )
  if (niche_col %in% colnames(meta)) {
    spot_df$niche <- as.character(meta[[niche_col]])
  } else {
    spot_df$niche <- NA_character_
  }

  valid_coord <- is.finite(spot_df$x) & is.finite(spot_df$y)
  coords <- as.matrix(spot_df[, c("x", "y")])
  result_list <- list()
  row_id <- 1

  for (cell_type_i in target_t_cells) {
    target_idx <- which(spot_df$cell_type == cell_type_i & valid_coord)
    if (length(target_idx) == 0) {
      warning(sample_label, " 中没有找到 ", cell_type_i, "。")
      next
    }

    for (target_i in target_idx) {
      distances <- sqrt((coords[, 1] - coords[target_i, 1])^2 + (coords[, 2] - coords[target_i, 2])^2)
      distances[!valid_coord] <- Inf
      if (!include_self) {
        distances[target_i] <- Inf
      }

      for (ring_i in seq_along(ring_labels)) {
        in_ring <- distances > ring_lower[[ring_i]] & distances <= ring_upper[[ring_i]]
        n_bins <- sum(in_ring)

        if (n_bins > 0) {
          ring_mean <- Matrix::rowMeans(expr_mat[, in_ring, drop = FALSE])
        } else {
          ring_mean <- rep(NA_real_, nrow(expr_mat))
          names(ring_mean) <- rownames(expr_mat)
        }

        scaled_expression <- as.numeric(ring_mean / global_gene_mean[rownames(expr_mat)])
        scaled_expression[!is.finite(scaled_expression)] <- NA_real_

        result_list[[row_id]] <- data.frame(
          sample = sample_label,
          target_cell_type = cell_type_i,
          target_cell_id = spot_df$cell_id[[target_i]],
          target_x = spot_df$x[[target_i]],
          target_y = spot_df$y[[target_i]],
          target_niche = spot_df$niche[[target_i]],
          ring_label = ring_labels[[ring_i]],
          ring_lower = ring_lower[[ring_i]],
          ring_upper = ring_upper[[ring_i]],
          gene = rownames(expr_mat),
          mean_expression = as.numeric(ring_mean),
          global_gene_mean = as.numeric(global_gene_mean[rownames(expr_mat)]),
          scaled_expression = scaled_expression,
          n_bins_in_ring = n_bins,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
      }
    }
  }

  cell_level <- bind_rows(result_list)
  sample_info <- data.frame(
    sample = sample_label,
    object_size = format(object.size(obj), units = "auto"),
    n_cells = ncol(obj),
    n_features_spatial = nrow(obj[[spatial_assay]]),
    default_assay = DefaultAssay(obj),
    assays = paste(Assays(obj), collapse = ";"),
    n_CTL = sum(meta[[cell_type_col]] == "CTL", na.rm = TRUE),
    n_Tex_T8 = sum(meta[[cell_type_col]] == "Tex_T8", na.rm = TRUE),
    n_TRM_T8 = sum(meta[[cell_type_col]] == "TRM_T8", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    cell_level = cell_level,
    gene_match = mutate(gene_match, sample = sample_label, .before = 1),
    sample_info = sample_info
  )
}

# -------------------------------
# Step 3: 创建输出目录并加载对象
# 说明：每次运行使用带时间戳的新目录，避免覆盖旧 checkpoint
# -------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load(input_rdata)
if (!exists(spa_object_name)) {
  stop("输入文件中没有找到对象：", spa_object_name)
}

spa_input <- get(spa_object_name)
if (inherits(spa_input, "Seurat")) {
  spa_list <- list(sample1 = spa_input)
} else if (is.list(spa_input)) {
  spa_list <- spa_input
} else {
  stop(spa_object_name, " 既不是 Seurat 对象，也不是 Seurat 对象列表。")
}

if (is.null(names(spa_list))) {
  names(spa_list) <- paste0("sample", seq_along(spa_list))
}

if (is.null(sample_indices)) {
  selected_spa <- spa_list
} else {
  selected_spa <- spa_list[sample_indices]
}

# -------------------------------
# Step 4: 逐样本计算 ring 基因梯度
# 说明：每个目标 T 细胞为中心计算环带表达，后续直接混合所有样本的中心细胞
# -------------------------------
result_list <- list()
for (sample_i in names(selected_spa)) {
  result_list[[sample_i]] <- calculate_gene_ring_one_sample(
    obj = selected_spa[[sample_i]],
    sample_label = sample_i
  )
  message("完成 ring 基因梯度计算：", sample_i)
}

gene_ring_cell_level <- bind_rows(lapply(result_list, `[[`, "cell_level")) %>%
  mutate(
    target_cell_type = factor(target_cell_type, levels = target_t_cells),
    ring_label = factor(ring_label, levels = ring_labels),
    gene = factor(gene, levels = target_genes)
  )
gene_match_summary <- bind_rows(lapply(result_list, `[[`, "gene_match"))
sample_info <- bind_rows(lapply(result_list, `[[`, "sample_info"))

if (nrow(gene_ring_cell_level) == 0) {
  stop("没有得到任何 ring 基因梯度结果，请检查目标 T 细胞命名、坐标列和 ring 设置。")
}

saveRDS(
  list(
    gene_ring_cell_level = gene_ring_cell_level,
    gene_match_summary = gene_match_summary,
    sample_info = sample_info
  ),
  file = file.path(output_dir, "checkpoint_01_gene_ring_cell_level.rds")
)

write.csv(
  gene_match_summary,
  file = file.path(output_dir, "gene_matching_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  sample_info,
  file = file.path(output_dir, "sample_object_info.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  gene_ring_cell_level,
  file = file.path(output_dir, "Tcell_gene_ring_cell_level.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------------
# Step 5: 直接混合汇总并计算 95% 置信区间
# 说明：不区分样本来源，直接把所有目标 T 细胞中心点混合后计算均值和 95% CI
# -------------------------------
pooled_summary <- gene_ring_cell_level %>%
  group_by(target_cell_type, ring_label, ring_lower, ring_upper, gene) %>%
  summarise(
    n_target_cells = n_distinct(target_cell_id),
    n_samples = n_distinct(sample),
    pooled_mean_expression = safe_mean(mean_expression),
    pooled_scaled_expression = safe_mean(scaled_expression),
    pooled_sd_scaled_expression = safe_sd(scaled_expression),
    pooled_se_scaled_expression = safe_se(scaled_expression),
    ci95_low = pooled_scaled_expression - 1.96 * pooled_se_scaled_expression,
    ci95_high = pooled_scaled_expression + 1.96 * pooled_se_scaled_expression,
    mean_bins_in_ring = safe_mean(n_bins_in_ring),
    .groups = "drop"
  ) %>%
  mutate(
    target_cell_type = factor(target_cell_type, levels = target_t_cells),
    ring_label = factor(ring_label, levels = ring_labels),
    gene = factor(gene, levels = target_genes)
  )

saveRDS(
  list(
    pooled_summary = pooled_summary
  ),
  file = file.path(output_dir, "checkpoint_02_gene_ring_pooled_95CI_summary.rds")
)
write.csv(
  pooled_summary,
  file = file.path(output_dir, "Tcell_gene_ring_pooled_95CI_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------------
# Step 6: 生成直接混合图
# 说明：总图使用所有样本目标 T 细胞直接混合后的 scaled expression，bar 为 95% 置信区间
# -------------------------------
cell_colors <- c(
  CTL = "#3C5488",
  Tex_T8 = "#E64B35",
  TRM_T8 = "#00A087"
)

if (is.null(plot_genes)) {
  plot_genes_use <- target_genes
} else {
  plot_genes_use <- intersect(plot_genes, as.character(unique(pooled_summary$gene)))
  missing_plot_genes <- setdiff(plot_genes, plot_genes_use)
  if (length(missing_plot_genes) > 0) {
    warning("以下 plot_genes 没有可绘制结果，已跳过：", paste(missing_plot_genes, collapse = ", "))
  }
}

if (length(plot_genes_use) == 0) {
  stop("plot_genes 没有匹配到任何可绘制基因，请检查基因名。")
}

plot_summary <- pooled_summary %>%
  filter(as.character(gene) %in% plot_genes_use) %>%
  mutate(gene = factor(as.character(gene), levels = plot_genes_use))

n_gene <- length(unique(na.omit(as.character(plot_summary$gene))))
plot_height <- max(8, ceiling(n_gene / plot_ncol) * plot_height_per_row)

p_line <- ggplot(
  plot_summary,
  aes(
    x = ring_label,
    y = pooled_scaled_expression,
    color = target_cell_type,
    group = target_cell_type
  )
) +
  geom_line(linewidth = 0.85, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  geom_errorbar(
    aes(
      ymin = ci95_low,
      ymax = ci95_high
    ),
    width = 0.15,
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  facet_wrap(~ gene, scales = "free_y", ncol = plot_ncol) +
  scale_color_manual(values = cell_colors) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "半径环带",
    y = "直接混合表达 / 样本全局均值",
    color = "中心 T 细胞",
    title = "三类 T 细胞邻域基因梯度（ring 模式，直接混合，95% CI）"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey70")
  )

ggsave(
  filename = file.path(output_dir, "Tcell_gene_ring_pooled_95CI_lineplot_selected_genes.pdf"),
  plot = p_line,
  width = plot_width,
  height = plot_height,
  limitsize = FALSE,
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Tcell_gene_ring_pooled_95CI_lineplot_selected_genes.png"),
  plot = p_line,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

heatmap_df <- plot_summary %>%
  group_by(gene) %>%
  mutate(gene_z = scale_for_plot(pooled_scaled_expression)) %>%
  ungroup()

p_heatmap <- ggplot(
  heatmap_df,
  aes(x = ring_label, y = target_cell_type, fill = gene_z)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_wrap(~ gene, ncol = plot_ncol) +
  scale_fill_gradient2(
    low = "#3C5488",
    mid = "white",
    high = "#E64B35",
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    x = "半径环带",
    y = "中心 T 细胞",
    fill = "基因内 z-score",
    title = "三类 T 细胞邻域基因梯度热图（ring 模式）"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey70")
  )

ggsave(
  filename = file.path(output_dir, "Tcell_gene_ring_pooled_95CI_heatmap_selected_genes.pdf"),
  plot = p_heatmap,
  width = plot_width,
  height = plot_height,
  limitsize = FALSE,
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Tcell_gene_ring_pooled_95CI_heatmap_selected_genes.png"),
  plot = p_heatmap,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

# -------------------------------
# Step 7: 保存运行信息
# 说明：记录对象规模、参数和 sessionInfo，方便之后复现或排查差异
# -------------------------------
run_info <- c(
  paste0("input_rdata: ", input_rdata),
  paste0("spa_object_name: ", spa_object_name),
  paste0("output_dir: ", output_dir),
  paste0("selected_samples: ", paste(names(selected_spa), collapse = ", ")),
  paste0("target_t_cells: ", paste(target_t_cells, collapse = ", ")),
  paste0("ring_labels: ", paste(ring_labels, collapse = ", ")),
  paste0("include_self: ", include_self),
  paste0("spatial_assay: ", spatial_assay),
  paste0("expression_layer: ", expression_layer),
  paste0("target_genes: ", paste(target_genes, collapse = ", ")),
  paste0("plot_genes: ", paste(plot_genes_use, collapse = ", ")),
  "",
  "sample_info:",
  capture.output(print(sample_info))
)
writeLines(run_info, con = file.path(output_dir, "run_info.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("完成：T 细胞基因梯度 ring 直接混合 95% CI 分析已输出到：", output_dir)
