# -------------------------------
# Step 0: 参数设置
# 说明：所有可修改参数集中放在开头，便于手动调整和复现
# -------------------------------
input_rdata <- "spa .rdata"
spa_object_name <- "spa"

output_root <- file.path("fig", "Tcell_Hallmark_radius_signal")
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_root, paste0("run_", run_id))

sample_indices <- NULL
spatial_assay <- "Spatial"
expression_layer <- "data"
coord_cols <- c("x", "y")
cell_type_col <- "type"
niche_col <- "nich4"

target_t_cells <- c("CTL", "Tex_T8", "TRM_T8")
radii <- c(10, 30, 50, 70, 90, 110)
radius_mode <- "ring"
include_self <- FALSE
radius_label_levels <- if (radius_mode == "ring") {
  paste0(c(0, head(radii, -1)), "-", radii)
} else {
  paste0("0-", radii)
}

species_name <- "Mus musculus"
target_pathways <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_FATTY_ACID_METABOLISM"
)

min_genes_per_pathway <- 5
scale_genes_within_sample <- TRUE
zscore_clip <- 10

plot_width <- 10
plot_height <- 6

# -------------------------------
# Step 1: 加载依赖包
# 说明：缺包时直接停止，避免中途运行到一半才报错
# -------------------------------
required_packages <- c("Seurat", "dplyr", "ggplot2", "Matrix", "msigdbr")
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
  library(msigdbr)
})

# -------------------------------
# Step 2: 工具函数
# 说明：只保留会反复使用的短函数，减少主流程中的重复代码
# -------------------------------
standardize_gene_name <- function(x) {
  # 统一基因名格式，用于处理大小写或符号差异导致的匹配失败
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

safe_mean <- function(x) {
  # 在邻域没有有效 spot/bin 时返回 NA，避免 NaN 混入后续汇总
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  # 中位数用于辅助判断是否被少数高值 spot/bin 拉动
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

safe_se <- function(x) {
  # 标准误用于画汇总曲线的误差线
  valid_x <- x[!is.na(x)]
  if (length(valid_x) <= 1) {
    return(NA_real_)
  }
  sd(valid_x) / sqrt(length(valid_x))
}

scale_for_plot <- function(x) {
  # 热图内部按 pathway 做 z-score，只用于显示相对模式
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

get_hallmark_gene_sets <- function(species_name, target_pathways) {
  # 从 msigdbr 提取 Hallmark 基因集，并兼容新版和旧版参数名
  msigdbr_args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% msigdbr_args) {
    msig_tbl <- msigdbr::msigdbr(species = species_name, collection = "H")
  } else {
    msig_tbl <- msigdbr::msigdbr(species = species_name, category = "H")
  }

  if (!all(c("gs_name", "gene_symbol") %in% colnames(msig_tbl))) {
    stop("msigdbr 返回结果中缺少 gs_name 或 gene_symbol 列，请检查 msigdbr 版本。")
  }

  msig_tbl <- msig_tbl %>%
    filter(gs_name %in% target_pathways) %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  missing_pathways <- setdiff(target_pathways, unique(msig_tbl$gs_name))
  if (length(missing_pathways) > 0) {
    stop("MSigDB 中没有找到以下 pathway：", paste(missing_pathways, collapse = ", "))
  }

  gene_sets <- split(msig_tbl$gene_symbol, msig_tbl$gs_name)
  gene_sets[target_pathways]
}

match_gene_sets_to_features <- function(gene_sets, feature_names) {
  # 将 MSigDB 基因名匹配到当前空间表达矩阵的 feature 名，避免 assay 命名差异造成空集
  feature_map <- data.frame(
    feature = feature_names,
    feature_key = standardize_gene_name(feature_names),
    stringsAsFactors = FALSE
  )
  feature_map <- feature_map[!duplicated(feature_map$feature_key), , drop = FALSE]

  matched_sets <- list()
  summary_list <- list()

  for (pathway_i in names(gene_sets)) {
    requested_genes <- unique(gene_sets[[pathway_i]])
    requested_key <- standardize_gene_name(requested_genes)
    matched_idx <- match(requested_key, feature_map$feature_key)
    matched_features <- unique(feature_map$feature[matched_idx[!is.na(matched_idx)]])
    missing_genes <- requested_genes[is.na(matched_idx)]

    matched_sets[[pathway_i]] <- matched_features
    summary_list[[pathway_i]] <- data.frame(
      pathway = pathway_i,
      requested_genes = length(requested_genes),
      matched_genes = length(matched_features),
      missing_genes = length(missing_genes),
      matched_feature_list = paste(matched_features, collapse = ";"),
      missing_gene_list = paste(head(missing_genes, 100), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  list(
    matched_sets = matched_sets,
    match_summary = bind_rows(summary_list)
  )
}

score_spot_pathways <- function(obj, sample_label, gene_sets_matched, score_col_map) {
  # 在 spot/bin 层面计算 pathway 模块分数，后续半径分析直接汇总这些空间分数
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

  spot_df <- data.frame(
    sample = sample_label,
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

  expr_mat <- GetAssayData(obj, assay = spatial_assay, layer = expression_layer)

  for (pathway_i in names(gene_sets_matched)) {
    score_col <- score_col_map[[pathway_i]]
    matched_genes <- gene_sets_matched[[pathway_i]]

    if (length(matched_genes) < min_genes_per_pathway) {
      warning(
        sample_label, " 的 ", pathway_i, " 只匹配到 ", length(matched_genes),
        " 个基因，低于阈值，分数记为 NA。"
      )
      spot_df[[score_col]] <- NA_real_
      next
    }

    pathway_mat <- expr_mat[matched_genes, , drop = FALSE]
    pathway_dense <- as.matrix(pathway_mat)

    if (scale_genes_within_sample) {
      # 每个基因先在样本内标准化，再取均值，减少高表达基因对模块分数的支配
      gene_mean <- rowMeans(pathway_dense, na.rm = TRUE)
      gene_sd <- apply(pathway_dense, 1, sd, na.rm = TRUE)
      gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- NA_real_
      z_mat <- sweep(pathway_dense, 1, gene_mean, "-")
      z_mat <- sweep(z_mat, 1, gene_sd, "/")
      z_mat[!is.finite(z_mat)] <- 0
      z_mat[z_mat > zscore_clip] <- zscore_clip
      z_mat[z_mat < -zscore_clip] <- -zscore_clip
      pathway_score <- colMeans(z_mat, na.rm = TRUE)
    } else {
      # 如需绝对表达强度，可关闭标准化，直接取基因集平均表达
      pathway_score <- Matrix::colMeans(pathway_mat, na.rm = TRUE)
    }

    spot_df[[score_col]] <- as.numeric(pathway_score)
  }

  spot_df
}

summarise_radius_signal <- function(spot_df, score_col_map) {
  # 以每个目标 T 细胞为中心，汇总指定半径内所有 spot/bin 的 pathway 分数
  if (!radius_mode %in% c("cumulative", "ring")) {
    stop("radius_mode 只能是 'cumulative' 或 'ring'。")
  }

  valid_coord <- is.finite(spot_df$x) & is.finite(spot_df$y)
  coords <- as.matrix(spot_df[, c("x", "y")])
  n_target_total <- sum(spot_df$cell_type %in% target_t_cells & valid_coord)

  if (n_target_total == 0) {
    warning(unique(spot_df$sample), " 中没有找到目标 T 细胞。")
    return(data.frame())
  }

  result_list <- vector("list", n_target_total * length(radii) * length(score_col_map))
  row_id <- 1

  for (cell_type_i in target_t_cells) {
    target_idx <- which(spot_df$cell_type == cell_type_i & valid_coord)
    if (length(target_idx) == 0) {
      warning(unique(spot_df$sample), " 中没有找到 ", cell_type_i, "。")
      next
    }

    for (target_i in target_idx) {
      distances <- sqrt((coords[, 1] - coords[target_i, 1])^2 + (coords[, 2] - coords[target_i, 2])^2)
      distances[!valid_coord] <- Inf
      if (!include_self) {
        distances[target_i] <- Inf
      }

      for (radius_i in seq_along(radii)) {
        radius_upper <- radii[[radius_i]]
        radius_lower <- 0
        if (radius_mode == "ring" && radius_i > 1) {
          radius_lower <- radii[[radius_i - 1]]
        }

        if (radius_mode == "ring") {
          in_radius <- distances > radius_lower & distances <= radius_upper
        } else {
          in_radius <- distances <= radius_upper
        }

        for (pathway_i in names(score_col_map)) {
          score_col <- score_col_map[[pathway_i]]
          score_values <- spot_df[[score_col]][in_radius]
          valid_scores <- score_values[!is.na(score_values)]

          result_list[[row_id]] <- data.frame(
            sample = spot_df$sample[[target_i]],
            target_cell_type = cell_type_i,
            target_cell_id = spot_df$cell_id[[target_i]],
            target_x = spot_df$x[[target_i]],
            target_y = spot_df$y[[target_i]],
            target_niche = spot_df$niche[[target_i]],
            radius = radius_upper,
            radius_lower = radius_lower,
            radius_label = paste0(radius_lower, "-", radius_upper),
            radius_mode = radius_mode,
            pathway = pathway_i,
            pathway_score_col = score_col,
            mean_score = safe_mean(score_values),
            median_score = safe_median(score_values),
            sum_score = ifelse(length(valid_scores) == 0, NA_real_, sum(valid_scores)),
            n_bins_in_radius = sum(in_radius),
            n_bins_with_score = length(valid_scores),
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1
        }
      }
    }
  }

  bind_rows(result_list[seq_len(row_id - 1)])
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
# Step 4: 提取 MSigDB Hallmark 基因集
# 说明：只保留缺氧、糖酵解和脂肪酸代谢三个 pathway
# -------------------------------
hallmark_gene_sets <- get_hallmark_gene_sets(species_name, target_pathways)
score_col_map <- setNames(
  paste0(tolower(gsub("^HALLMARK_", "", target_pathways)), "_score"),
  target_pathways
)

saveRDS(
  hallmark_gene_sets,
  file = file.path(output_dir, "checkpoint_01_hallmark_gene_sets.rds")
)

# -------------------------------
# Step 5: 逐样本计算 spot/bin 层面 pathway score
# 说明：每个样本单独标准化，避免不同样本测序深度或表达尺度直接混合
# -------------------------------
spot_score_list <- list()
match_summary_list <- list()
sample_info_list <- list()

for (sample_i in names(selected_spa)) {
  obj <- selected_spa[[sample_i]]

  if (!inherits(obj, "Seurat")) {
    warning(sample_i, " 不是 Seurat 对象，已跳过。")
    next
  }

  feature_names <- rownames(obj[[spatial_assay]])
  match_result <- match_gene_sets_to_features(hallmark_gene_sets, feature_names)
  matched_sets <- match_result$matched_sets
  match_summary <- match_result$match_summary %>%
    mutate(sample = sample_i, .before = 1)
  match_summary_list[[sample_i]] <- match_summary

  spot_score_list[[sample_i]] <- score_spot_pathways(
    obj = obj,
    sample_label = sample_i,
    gene_sets_matched = matched_sets,
    score_col_map = score_col_map
  )

  meta <- obj@meta.data
  sample_info_list[[sample_i]] <- data.frame(
    sample = sample_i,
    object_size = format(object.size(obj), units = "auto"),
    n_cells = ncol(obj),
    n_features_spatial = length(feature_names),
    default_assay = DefaultAssay(obj),
    assays = paste(Assays(obj), collapse = ";"),
    n_CTL = sum(meta[[cell_type_col]] == "CTL", na.rm = TRUE),
    n_Tex_T8 = sum(meta[[cell_type_col]] == "Tex_T8", na.rm = TRUE),
    n_TRM_T8 = sum(meta[[cell_type_col]] == "TRM_T8", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  message("完成 spot/bin pathway 打分：", sample_i)
}

spot_scores <- bind_rows(spot_score_list)
gene_match_summary <- bind_rows(match_summary_list)
sample_info <- bind_rows(sample_info_list)

if (nrow(spot_scores) == 0) {
  stop("没有生成任何 spot/bin pathway score，请检查 sample_indices、输入对象和 Spatial assay。")
}

saveRDS(
  list(
    spot_scores = spot_scores,
    gene_match_summary = gene_match_summary,
    sample_info = sample_info
  ),
  file = file.path(output_dir, "checkpoint_02_spot_pathway_scores.rds")
)

write.csv(
  gene_match_summary,
  file = file.path(output_dir, "Tcell_Hallmark_gene_matching.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  sample_info,
  file = file.path(output_dir, "sample_object_info.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------------
# Step 6: 计算三类 T 细胞半径邻域信号
# 说明：以每个 CTL、Tex_T8、TRM_T8 为中心，汇总半径内 spot/bin 的 pathway score
# -------------------------------
radius_result_list <- lapply(
  split(spot_scores, spot_scores$sample),
  summarise_radius_signal,
  score_col_map = score_col_map
)

radius_cell_level <- bind_rows(radius_result_list) %>%
  mutate(
    target_cell_type = factor(target_cell_type, levels = target_t_cells),
    radius_label = factor(radius_label, levels = radius_label_levels),
    pathway = factor(pathway, levels = target_pathways),
    pathway_label = gsub("^HALLMARK_", "", as.character(pathway))
  )

if (nrow(radius_cell_level) == 0) {
  stop("没有得到任何半径邻域结果，请检查目标 T 细胞命名、坐标列和半径设置。")
}

saveRDS(
  radius_cell_level,
  file = file.path(output_dir, "checkpoint_03_radius_cell_level_results.rds")
)

write.csv(
  radius_cell_level,
  file = file.path(output_dir, "Tcell_radius_pathway_cell_level.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------------
# Step 7: 汇总样本内和跨样本结果
# 说明：先对每个样本内目标 T 细胞取均值，再跨样本汇总，避免样本细胞数差异主导整体结果
# -------------------------------
sample_summary <- radius_cell_level %>%
  group_by(sample, target_cell_type, radius, radius_lower, radius_label, radius_mode, pathway, pathway_label) %>%
  summarise(
    n_target_cells = n_distinct(target_cell_id),
    mean_score = safe_mean(mean_score),
    median_target_score = safe_median(mean_score),
    se_target_score = safe_se(mean_score),
    mean_bins_in_radius = safe_mean(n_bins_in_radius),
    mean_bins_with_score = safe_mean(n_bins_with_score),
    .groups = "drop"
  )

overall_summary <- sample_summary %>%
  group_by(target_cell_type, radius, radius_lower, radius_label, radius_mode, pathway, pathway_label) %>%
  summarise(
    n_samples = sum(!is.na(mean_score)),
    mean_score = safe_mean(mean_score),
    median_sample_score = safe_median(mean_score),
    se_across_samples = safe_se(mean_score),
    mean_target_cells = safe_mean(n_target_cells),
    mean_bins_in_radius = safe_mean(mean_bins_in_radius),
    .groups = "drop"
  )

saveRDS(
  list(
    sample_summary = sample_summary,
    overall_summary = overall_summary
  ),
  file = file.path(output_dir, "checkpoint_04_radius_summary.rds")
)

write.csv(
  sample_summary,
  file = file.path(output_dir, "Tcell_radius_pathway_sample_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  overall_summary,
  file = file.path(output_dir, "Tcell_radius_pathway_overall_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -------------------------------
# Step 8: 生成汇总图
# 说明：曲线展示半径趋势，热图展示三类 T 细胞在不同 pathway 和半径下的相对强弱
# -------------------------------
cell_colors <- c(
  CTL = "#3C5488",
  Tex_T8 = "#E64B35",
  TRM_T8 = "#00A087"
)

p_line <- ggplot(
  overall_summary,
  aes(x = radius_label, y = mean_score, color = target_cell_type, group = target_cell_type)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = mean_score - se_across_samples, ymax = mean_score + se_across_samples),
    width = 0.15,
    linewidth = 0.5,
    na.rm = TRUE
  ) +
  facet_wrap(~ pathway_label, scales = "free_y") +
  scale_color_manual(values = cell_colors) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "半径环带",
    y = "平均空间 pathway score",
    color = "中心 T 细胞",
    title = "三类 T 细胞周围 Hallmark pathway 信号随环带变化"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey70")
  )

ggsave(
  filename = file.path(output_dir, "Tcell_radius_pathway_lineplot.pdf"),
  plot = p_line,
  width = plot_width,
  height = plot_height,
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Tcell_radius_pathway_lineplot.png"),
  plot = p_line,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

heatmap_df <- overall_summary %>%
  group_by(pathway) %>%
  mutate(pathway_z = scale_for_plot(mean_score)) %>%
  ungroup()

p_heatmap <- ggplot(
  heatmap_df,
  aes(x = radius_label, y = target_cell_type, fill = pathway_z)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~ pathway_label) +
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
    fill = "pathway 内 z-score",
    title = "三类 T 细胞邻域 Hallmark pathway 相对信号"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey70")
  )

ggsave(
  filename = file.path(output_dir, "Tcell_radius_pathway_heatmap.pdf"),
  plot = p_heatmap,
  width = plot_width,
  height = plot_height,
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Tcell_radius_pathway_heatmap.png"),
  plot = p_heatmap,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

# -------------------------------
# Step 9: 保存运行信息
# 说明：记录对象规模、参数和 sessionInfo，方便之后复现或排查差异
# -------------------------------
run_info <- c(
  paste0("input_rdata: ", input_rdata),
  paste0("spa_object_name: ", spa_object_name),
  paste0("output_dir: ", output_dir),
  paste0("selected_samples: ", paste(names(selected_spa), collapse = ", ")),
  paste0("target_t_cells: ", paste(target_t_cells, collapse = ", ")),
  paste0("radii: ", paste(radii, collapse = ", ")),
  paste0("radius_labels: ", paste(radius_label_levels, collapse = ", ")),
  paste0("radius_mode: ", radius_mode),
  paste0("include_self: ", include_self),
  paste0("spatial_assay: ", spatial_assay),
  paste0("expression_layer: ", expression_layer),
  paste0("scale_genes_within_sample: ", scale_genes_within_sample),
  paste0("target_pathways: ", paste(target_pathways, collapse = ", ")),
  "",
  "sample_info:",
  capture.output(print(sample_info))
)
writeLines(run_info, con = file.path(output_dir, "run_info.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("完成：三类 T 细胞 Hallmark 半径信号分析已输出到：", output_dir)
