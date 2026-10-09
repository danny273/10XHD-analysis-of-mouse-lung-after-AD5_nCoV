# 1. 环境准备与数据过滤
library(Seurat); library(tidyverse); library(ggpubr); library(eulerr); library(cowplot)
setwd("D:/atlas/")
load("balt/3月ibalt.rdata")
  load("balt/bulk_new.rdata")
head(pbmc)
  DimPlot(pbmc,label=T)
# 提取浆细胞子集并要求 Jchain 表达
pbmc3 <- subset(pbmc, group == "pla")
pbmc3 <- subset(pbmc3, subset = Jchain > 0)
exp <- FetchData(pbmc3, vars = c("Igha", "Ighm", "Ighg2c", "Ighg2b", "Ighg1", "Ighg3", "Jchain"))

# 2. 自动阈值计算与分类逻辑 (利用密度图低谷)
get_threshold <- function(x, valley_rank = 1) {
  if (length(unique(x)) < 3) return(Inf)
  d <- density(x, na.rm = TRUE)
  valleys <- which(diff(sign(diff(d$y))) == 2) + 1
  peaks   <- which(diff(sign(diff(d$y))) == -2) + 1
  if (length(valleys) == 0 || length(peaks) == 0) return(1)
  candidate_valleys <- valleys[valleys > peaks[1]]
  if (length(candidate_valleys) == 0) return(1)
  res <- d$x[candidate_valleys[valley_rank]]
  return(ifelse(is.na(res) || res < 1, 1, res))
}

genes_to_check <- intersect(c("Igha", "Ighm", "Ighg2c", "Ighg2b", "Ighg1", "Ighg3"), colnames(exp))

# Igha 的第一个谷偏保守，因此单独改为第二个谷；其余 Ig 基因仍使用第一个谷
threshold_valley_rank <- setNames(rep(1, length(genes_to_check)), genes_to_check)
threshold_valley_rank[intersect("Igha", names(threshold_valley_rank))] <- 2

ths <- vapply(
  genes_to_check,
  function(gene) get_threshold(exp[[gene]], valley_rank = threshold_valley_rank[[gene]]),
  numeric(1)
)

# -------------------------------
# Step 2.1: 绘制 Ig 阈值确定图
# 说明：这里只画密度曲线和最终阈值线，方便人工检查 Ig 分型是否合理
# -------------------------------
ig_threshold_dir <- file.path("fig", "三月pla_ig_threshold")
dir.create(ig_threshold_dir, recursive = TRUE, showWarnings = FALSE)

ig_density_df <- bind_rows(lapply(names(ths), function(gene) {
  gene_density <- density(exp[[gene]], na.rm = TRUE)
  tibble(gene = gene, x = gene_density$x, y = gene_density$y)
}))

ig_threshold_df <- tibble(
  gene = names(ths),
  threshold = as.numeric(ths),
  valley_rank = as.integer(threshold_valley_rank[names(ths)]),
  threshold_label = ifelse(is.finite(threshold), round(threshold, 2), Inf),
  positive_n = sapply(names(ths), function(gene) sum(exp[[gene]] > ths[[gene]], na.rm = TRUE)),
  total_n = nrow(exp)
) %>%
  mutate(
    positive_pct = positive_n / total_n,
    plot_label = paste0("第", valley_rank, "个谷\n阈值=", threshold_label, "\n阳性=", positive_n, " (", scales::percent(positive_pct, accuracy = 0.1), ")")
  )
write_csv(ig_threshold_df, file.path(ig_threshold_dir, "ig_threshold_values.csv"))

ig_threshold_plot <- ggplot(ig_density_df, aes(x = x, y = y)) +
  geom_area(fill = "#E8EEF7", alpha = 0.75) +
  geom_line(color = "#2F5C7C", size = 0.55) +
  geom_vline(
    data = ig_threshold_df %>% filter(is.finite(threshold)),
    aes(xintercept = threshold),
    color = "#D55E00",
    linetype = "dashed",
    size = 0.65
  ) +
  geom_text(
    data = ig_threshold_df %>% filter(is.finite(threshold)),
    aes(x = threshold, y = Inf, label = plot_label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 3.1,
    color = "#8B3A14"
  ) +
  facet_wrap(~ gene, scales = "free", ncol = 3) +
  labs(
    title = "Ig 基因表达阈值确定图",
    x = "表达量",
    y = "密度"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey92", color = "grey75"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(ig_threshold_dir, "ig_threshold_density.png"), ig_threshold_plot, width = 10, height = 6, dpi = 300)
ggsave(file.path(ig_threshold_dir, "ig_threshold_density.pdf"), ig_threshold_plot, width = 10, height = 6)

pbmc3$ig <- apply(exp, 1, function(row) {
  labels <- c()
  if (row["Igha"] > ths["Igha"]) labels <- c(labels, "IGA")
  if (any(row[intersect(c("Ighm", "Ighd"), names(ths))] > ths[intersect(c("Ighm", "Ighd"), names(ths))])) labels <- c(labels, "IGM")
  if (any(row[intersect(c("Ighg2c", "Ighg2b", "Ighg1", "Ighg3"), names(ths))] > ths[intersect(c("Ighg2c", "Ighg2b", "Ighg1", "Ighg3"), names(ths))])) labels <- c(labels, "IGG")
  if (length(labels) > 0) return(paste(labels, collapse = "_")) else return("Other")
})


###########################################################################################################################
# 3. 欧拉图绘制 (按时间点比例缩放)############################################################
my_fills <- c("IGA" = "#FF9999", "IGG" = "#99CCFF", "IGM" = "#99FF99")
df_euler <- as.data.frame(table(pbmc3$time, pbmc3$ig)) %>% 
  filter(Freq > 0) %>% mutate(Var2 = str_replace_all(Var2, "_", "&"))

plot_list <- list(); counts <- c()
for(tp in unique(df_euler$Var1)) {
  dat <- setNames(df_euler$Freq[df_euler$Var1==tp], df_euler$Var2[df_euler$Var1==tp])
  counts[tp] <- sum(dat)
  plot_list[[tp]] <- plot(euler(dat, shape = "ellipse"), quantities = TRUE, fills = my_fills, alpha = 0.6, edges = list(col = "white"))
}

final_euler <- plot_grid(plotlist = plot_list, labels = names(counts), scale = sqrt(counts/max(counts)), ncol = 3)
print(final_euler)


############################################################################################################
# 4. 空间相关性散点图 (DNBS)#############################################################################
# 颜色定义：单一色 + 高饱和中间色
my_palette <- c("IGA"="#FF9999", "IGG"="#99CCFF", "IGM"="#99FF99", 
                "IGA_IGG"="#D127D1", "IGA_IGM"="#FF9500", "IGM_IGG"="#00C5CD", "IGA_IGM_IGG"="#404040")

# 提取目标基因（以Igha为例，请根据需要修改变量名）
target_gene <- "Cxcl12" 
meta=pbmc3@meta.data
meta$plot_gene <- FetchData(pbmc3, vars = target_gene)[[1]]

ggplot(meta, aes(x = dist, y = plot_gene)) +
  geom_point(aes(color = ig), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", fill = "grey80", size = 0.8) +
  scale_color_manual(values = my_palette) +
  stat_cor(method = "pearson", label.x.npc = "right", label.y.npc = "top", 
           hjust = 1, color = "red3", fontface = "bold.italic", size = 5) +
  labs(title = bquote("Correlation: " ~ italic(.(target_gene)) ~ " vs DNBS"),
       x = "Distance to Nearest Bronchovascular Structure (μm)",
       y = bquote(italic(.(target_gene)) ~ " Expression"), color = "Ig Phenotype") +
  theme_classic() + theme(plot.title = element_text(face = "bold"), axis.title = element_text(face = "bold"))

####数量与结构关系
# 1. 确保数据准备好
meta$num2 <- log10(meta$num + 1)
meta$dist2 <- log10(meta$dist + 1)
# 2. 绘图
library(ggplot2)
library(ggpubr)
library(scales) # 用于处理坐标轴刻度

library(ggpubr)

library(ggplot2)
library(ggpubr)
library(patchwork)

num_plot_dir <- file.path("fig", "三月pla_num_dist")
dir.create(num_plot_dir, recursive = TRUE, showWarnings = FALSE)
num_label_y <- max(meta$num, na.rm = TRUE) * 1.05

# 1. 统计比较图：展示复合型 niche 细胞数量(num)显著更多
ggplot(meta, aes(x = ig, y = num, fill = ig)) +
  # 用小提琴图看分布，箱线图看中位数
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", size = 0.6) +
  # 关键：由于 num 差异可能很大，使用对数轴展示更清晰
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000)) +
  scale_fill_manual(values = my_palette) +
  # 自动添加显著性检验（例如：复合型 vs 其他）
  stat_compare_means(label = "p.signif", method = "t.test", 
                     ref.group = ".all.", label.y = num_label_y) + 
  labs(x = "Niche Phenotype", y = "Cell Count (Log-scale)") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

# 2. 空间相关性主图 (沿用之前的优化代码)
p_main <- ggplot(meta, aes(x = dist, y = num)) +
  geom_point(aes(color = ig), size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black", fill = "grey90") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000)) +
  scale_color_manual(values = my_palette) +
  stat_cor(method = "pearson", color = "red3", fontface = 4) +
  labs(x = "Distance to DNBS (μm)", y = "Cell Count (Log-scale)") +
  theme_classic()
print(p_main)
ggsave(file.path(num_plot_dir, "num_vs_dist_by_ig.png"), p_main, width = 7, height = 5, dpi = 300)

table(meta$time)
my_time_palette <- setNames(
  scales::hue_pal()(length(unique(meta$time))),
  sort(unique(meta$time))
)

p_time <- ggplot(meta, aes(x = time, y = num, fill = time)) +
  # 用小提琴图看分布，箱线图看中位数
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", size = 0.6) +
  # 关键：由于 num 差异可能很大，使用对数轴展示更清晰
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000)) +
  scale_fill_manual(values = my_time_palette) +
  # 自动添加显著性检验（例如：复合型 vs 其他）
  stat_compare_means(label = "p.signif", method = "t.test", 
                     ref.group = ".all.", label.y = num_label_y) + 
  labs(x = "Time", y = "Cell Count (Log-scale)") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
print(p_time)
ggsave(file.path(num_plot_dir, "num_by_time.png"), p_time, width = 7, height = 5, dpi = 300)
#################################################################################################################################
