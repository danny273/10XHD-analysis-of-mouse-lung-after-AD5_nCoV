
library(data.table)
library(pheatmap)

setwd("D:/atlas/")
path="combined_plots60_30_2_0.9_cosin"
# 高效读取数据 (比read.csv快得多)
data <- fread(paste0(path,"/leiden_labels_with_index.csv"))
sor=data$source
sor[sor=="0"] <- "d001"  
sor[sor=="1"] <- "d101"
sor[sor=="2"] <- "d102"
sor[sor=="3"] <- "d011"
sor[sor=="4"] <- "d211"
sor[sor=="5"] <- "d212"
sor[sor=="6"] <- "d451"
sor[sor=="7"] <- "d452"
data$source <- sor
# 移除索引列（假设第一列是自动生成的索引）
data[, V1 := NULL]  # 根据实际列名调整，例如如果列名是"X"则用 data[, X := NULL]
# 提取关键列并转换为data.table
heat=function(data,path){
  dt <- data[, .(type, merged_cluster)]
  
  # 使用data.table的dcast直接生成计数矩阵 (比tidyr/dplyr快得多)
  d1 <- dcast(dt, merged_cluster ~ type, 
              fun.aggregate = length, 
              drop = FALSE, 
              value.var = "type", 
              fill = 0)
  
  # 转换为数据框并设置行名
  d1_df <- as.data.frame(d1)
  rownames(d1_df) <- d1_df$merged_cluster
  d1_df$merged_cluster <- NULL
  
  # 筛选列总和>10的列
  col_sums <- colSums(d1_df)
  d2 <- d1_df[, col_sums > 10, drop = FALSE]
  
  # 绘制热图（如果数据过大建议移除display_numbers参数）
  png(path, width = 1200, height = 1200)
  pheatmap(d2,cluster_rows = T,cluster_cols = T,show_rownames = T,show_colnames = T,
           color = colorRampPalette(c("blue","white", "red"))(50),scale = "column",
           fontsize_row = 10,fontsize_col = 10,cellwidth = 20,cellheight = 20,display_numbers = d2,fontsize=12,
           main = "Cell type distribution in each niche",border_color = NA)
  dev.off()
}
heat(data,paste0(path,"/heatmap.png"))

write.csv(data,paste0(path,"/leiden.csv"),row.names = F)

for (sam in unique(data$source)) {
  dat=data[data$source==sam,]
  heat(dat,paste0(path,"/",sam,".png"))
}
########################################################################################################
library(ggplot2)
library(dplyr)

# 计算每个source中各niche的比例
data_summary <- data %>%
  group_by(source, merged_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(source) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 绘制堆叠柱状图
ggplot(data_summary, 
       aes(x = factor(source), 
           y = proportion, 
           fill = merged_cluster)) +
  geom_col(width = 0.7) +
  labs(x = "Source", 
       y = "Proportion", 
       fill = "Niche",
       title = "Niche Composition by Source") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())
###############################################################################################################
##计算各niche的比例随source变化的柱状图，#####################################################
library(tidyverse)

# 创建保存目录
if (!dir.exists(paste0(path,"/","niche_ration"))) {
  dir.create(paste0(path,"/","niche_ration"))
}

# 补全所有source和niche组合，缺失比例设为0
complete_data <- data_summary %>%
  complete(source, merged_cluster, fill = list(proportion = 0))

# 获取所有唯一的niche
niche_list <- unique(complete_data$merged_cluster)

# 为每个niche绘制图形
for (niche in niche_list) {
  p <- complete_data %>%
    filter(merged_cluster == niche) %>%
    ggplot(aes(x = source, y = proportion)) +
    geom_col(fill = "#0072B2", width = 0.7) +  # 使用蓝色系填充
    labs(
      title = paste("Proportion Distribution:", niche),
      subtitle = "Across Different Sources",
      x = "Data Source",
      y = "Proportion"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # 优化y轴显示范围
  
  # 保存图像
  ggsave(
    filename = paste0(path,"/","niche_ration/", niche, ".png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}
