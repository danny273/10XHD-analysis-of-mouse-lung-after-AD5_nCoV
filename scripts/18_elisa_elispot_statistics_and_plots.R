# -------------------------------
# Step 1: 参数与依赖
# 说明：这里只读取已经整理好的 ELISA/ELISpot CSV，用于复现 2 行 8 图柱状散点图
# -------------------------------
library(tidyverse)
library(patchwork)
library(ggpubr)

project_dir <- if (file.exists("elisa.csv")) "." else "单免13m"
input_elisa <- file.path(project_dir, "elisa.csv")
input_elispot <- file.path(project_dir, "elisplot.csv")
output_dir <- file.path(project_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Con", "i.m.", "i.t.")
group_colors <- c("Con" = "#4D4D4D", "i.m." = "#D55E00", "i.t." = "#CC79A7")


# -------------------------------
# Step 2: 数据读取与长表整理
# 说明：统一 ELISA 和 ELISpot 的格式，便于用同一个绘图语法画四个 panel
# -------------------------------
elisa_long <- read_csv(input_elisa, show_col_types = FALSE) %>%
  rename(group = `组别`) %>%
  pivot_longer(-group, names_to = "metric", values_to = "value") %>%
  mutate(
    assay = "ELISA",
    group = factor(group, levels = group_levels),
    value = as.numeric(value)
  ) %>%
  filter(!is.na(group), !is.na(value))

elispot_long <- read_csv(input_elispot, show_col_types = FALSE) %>%
  select(`组别`, `BM-IgG`, `BM-IgA`, `Lung-IgG`, `Lung-IgA`) %>%
  rename(group = `组别`) %>%
  pivot_longer(-group, names_to = "metric", values_to = "value") %>%
  mutate(
    assay = "ELISpot",
    group = factor(group, levels = group_levels),
    value = as.numeric(str_trim(value))
  ) %>%
  filter(!is.na(group), !is.na(value))

plot_df <- bind_rows(elisa_long, elispot_long)
write_csv(plot_df, file.path(output_dir, "elisa_elispot_long.csv"))


# -------------------------------
# Step 3: Welch t tests on assay-appropriate values
# ELISA uses log10(value); ELISpot tests use log10(value + 0.5) to handle zeros.
# -------------------------------
run_welch_test <- function(metric_name, reference_group) {
  metric_df <- plot_df %>%
    filter(metric == metric_name, group %in% c(reference_group, "i.t.")) %>%
    mutate(
      test_value = if_else(
        assay == "ELISA",
        log10(value),
        log10(value + 0.5)
      )
    )

  values_it <- metric_df %>% filter(group == "i.t.") %>% pull(test_value)
  values_reference <- metric_df %>% filter(group == reference_group) %>% pull(test_value)
  raw_it <- metric_df %>% filter(group == "i.t.") %>% pull(value)
  raw_reference <- metric_df %>% filter(group == reference_group) %>% pull(value)
  test_result <- t.test(values_it, values_reference, var.equal = FALSE)

  tibble(
    metric = metric_name,
    assay = first(metric_df$assay),
    transformation = if_else(first(metric_df$assay) == "ELISA", "log10(value)", "log10(value + 0.5)"),
    test_method = "two-sided Welch t-test",
    comparison = paste("i.t. vs", reference_group),
    group1 = reference_group,
    group2 = "i.t.",
    n_1 = length(values_it),
    n_2 = length(values_reference),
    mean_1 = mean(raw_it),
    mean_2 = mean(raw_reference),
    mean_test_scale_1 = mean(values_it),
    mean_test_scale_2 = mean(values_reference),
    t_stat = unname(test_result$statistic),
    p_value = test_result$p.value
  )
}

test_stats <- map_dfr(unique(plot_df$metric), function(metric_name) {
  map_dfr(c("Con", "i.m."), function(reference_group) {
    run_welch_test(metric_name, reference_group)
  })
}) %>%
  mutate(
    p_label = map_chr(p_value, function(p) format.pval(p, digits = 2, eps = 0.0001))
  )

write_csv(
  test_stats %>% select(-group1, -group2, -p_label),
  file.path(output_dir, "ordinary_t_test_p_values.csv")
)


# -------------------------------
# Step 4: 单 panel 绘图
# 说明：柱子显示均值，误差线为 SEM，叠加每只小鼠的原始点；空心柱突出原始点分布
# -------------------------------
make_panel <- function(df, metric_name, y_label, y_max = NULL) {
  panel_df <- df %>%
    filter(metric == metric_name) %>%
    mutate(
      # Keep the plotted values on/above the biological zero baseline.
      plot_value = pmax(if_else(assay == "ELISA", log10(value), value), 0)
    )
  is_elisa <- first(panel_df$assay) == "ELISA"
  value_max <- max(panel_df$plot_value, na.rm = TRUE)
  lower_bracket_offset <- if (is_elisa) 0.25 else value_max * 0.12
  upper_bracket_offset <- if (is_elisa) 0.55 else value_max * 0.28
  panel_stats <- test_stats %>%
    filter(metric == metric_name) %>%
    mutate(
      y_position = value_max + if_else(
        group1 == "Con",
        lower_bracket_offset,
        upper_bracket_offset
      )
    )

  if (!is.null(y_max)) {
    y_max <- max(y_max, max(panel_stats$y_position, na.rm = TRUE) * 1.08)
  }

  p <- ggplot(panel_df, aes(x = group, y = plot_value, color = group)) +
    stat_summary(fun = mean, geom = "col", width = 0.66, fill = NA, linewidth = 0.85) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.24,
      linewidth = 0.75
    ) +
    # Horizontal jitter only; keep the measured y values unchanged.
    geom_jitter(width = 0.12, height = 0, size = 2.7, alpha = 0.95, stroke = 0) +
    stat_pvalue_manual(
      panel_stats,
      label = "p_label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y_position",
      size = 3.2,
      color = "black",
      bracket.size = 0.45,
      tip.length = 0.015,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(title = metric_name, x = NULL, y = y_label) +
    theme_classic(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.55),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = unit(2.2, "mm"),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.title.y = element_text(color = "black", size = 12, face = "bold", margin = margin(r = 6)),
      plot.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 10, 8, 10)
    )

  if (!is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max), clip = "on")
  }

  p
}


# -------------------------------
# Step 5: 复现参考图布局
# 说明：8 个指标全部绘制；上排为 ELISA，下排为 ELISpot
# -------------------------------
p1 <- make_panel(plot_df, "Serum-IgG", expression(log[10]~"Serum IgG"), y_max = 6.3)
p2 <- make_panel(plot_df, "Serum-IgA", expression(log[10]~"Serum IgA"), y_max = 3.8)
p3 <- make_panel(plot_df, "BALF-IgG", expression(log[10]~"BALF IgG"), y_max = 3.0)
p4 <- make_panel(plot_df, "BALF-IgA", expression(log[10]~"BALF IgA"), y_max = 5.2)

p5 <- make_panel(plot_df, "BM-IgG", "BM IgG ASC", y_max = 18)
p6 <- make_panel(plot_df, "BM-IgA", "BM IgA ASC", y_max = 24)
p7 <- make_panel(plot_df, "Lung-IgG", "Lung IgG ASC", y_max = 18)
p8 <- make_panel(plot_df, "Lung-IgA", "Lung IgA ASC", y_max = 8)

final_plot <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8) +
  plot_layout(widths = c(1, 1, 1, 1), heights = c(1, 1)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(output_dir, "elisa_elispot_8panel.png"), final_plot, width = 14.4, height = 7.2, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "elisa_elispot_8panel.pdf"), final_plot, width = 14.4, height = 7.2, bg = "white")
