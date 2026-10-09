
source("D:/atlas/配色.r")
spa211=spa[[5]]
# spa211=spa[[5]]
# meta=spam@meta.data
# table(meta$nich4)
# spam@meta.data$nich4 <- recode(spam@meta.data$nich4, "CN3_iBALT" = "CN3_BALCs")
# spa211@meta.data$nich4 <- recode(spa211@meta.data$nich4, "CN3_iBALT" = "CN3_BALCs")
# table(spa211$nich4)




# lev=levels(spa211$nich4)
# nb_colors <- length(lev)
# my_colors <- colorRampPalette(brewer.pal(12, "Set2"))(nb_colors)
# lev <- c(
#  "CN1_AM",
#  "CN2_ISG",
#  "CN3_BALCs",
#  "CN4_Lymphatic",
#  "CN5_Plasma",
#  "CN6_AdvFib",
# "CN7_Airway",
#  "CN8_Artery",
#  "CN9_SMC",
#  "CN10_Capillary",
#  "CN11_Alveolar"
# )
names(my_colors)=lev
p=DimPlot(spa[[6]],group.by = "nich4",raster=FALSE,cols = niche_col)
p
ggsave(p,height = 28,width = 30,file="d212.pdf")
#############################################################################
meta211=spa211@meta.data
nich4=meta211$nich4
table(nich4)
nich4=as.character(nich4)
nich4[nich4%in%c("CN7_Airway","CN8_Artery")]="struct"
meta211$nich4=nich4
knn_nich4=function(data,k=5,typeid=NULL,nihc4="blood"){
  library(FNN)
  # 提取坐标
  if (!is.null(typeid)) {
    cell1_data <-filter(data, nich4==typeid)
  } else {
    cell1_data <- data
  }
  cell2_data <- filter(data, nich4==nihc4)
  # 确保cell2细胞存在
  if (nrow(cell2_data) == 0) {
    stop("数据中没有", cell2, "细胞，无法计算距离。")
  }
  # 提取坐标
  cell1_coords <- cell1_data[, c("x", "y")]
  cell2_coords <- cell2_data[, c("x", "y")]
  # 计算每个细胞到最近的k个邻居的距离
  nearest_dist <- knnx.dist(
    query = cell1_coords,  # 需要查询的点
    data = cell2_coords,  # 参考点
    k = k # 找最近的k个邻居
  )
  # 将结果添加到inf数据中
  nearest_dist=as.data.frame(nearest_dist)
  cell1_data$nearest_dist <-rowSums(nearest_dist)
  return(cell1_data)
}
a=knn_nich4(meta211,nihc4="struct")

########叠加密度图###########################################################3
#########变量含义示意图
#major_colors <- c(AT="#D4CCD5", B="#A6DCC2", Mo="#C7D267", DC="#E8B374", Mac="#DED3B3", T="#BFBBD9", Gran="#EBD2DF", Airway="#C8A5C8", Endo="#BC82BD", Stroma="#D2EBBA", Fib="#FFED6F")
DimPlot(spa211,group.by = "type2",cols =maj_col,raster=FALSE)
ggsave(filename = "d211_type2.pdf",width = 40,height = 40)
spa211$dist_tra=a$nearest_dist
DimPlot(spa211,reduction ="spatial" )
FeaturePlot(spa211,reduction ="spatial",features = "dist_tra",raster=FALSE)
ggsave(filename = "fig/空间建模.tif",dpi = 600)
ggsave(filename = "fig/空间建模.pdf")
####################################################################
##新的密度图#####################################################
meta=spam@meta.data
nich4=meta$nich4
nich4=as.character(nich4)
nich4[nich4%in%c("CN7_Airway","CN8_Artery")]="struct"
source=c("d101","d102","d211","d212","d451","d452")
meta$nich4=nich4
#metat=subset(meta,type2=="T")
cellt=c("CN1_AM","CN2_ISG","CN3_BALCs","CN4_Lymphatic",
        "CN5_Plasma")
for (source_id in source) {
  print(source_id)
  
  # 【关键修改 1】：在每个 source 开始前初始化/清空列表
  density <- list() 
  
  meta211 = filter(meta, source == source_id)
  a = knn_nich4(meta211, nihc4 = "struct")
  
  for (cell_t in cellt) {
    dis = subset(a, nich4 == cell_t)
    
    # 检查是否有数据，防止 rbind 出错
    if(nrow(dis) > 0) {
      df = data.frame(Group = cell_t, Distance = dis$nearest_dist)
      density[[cell_t]] = df
    }
  }
  
  # 【关键修改 2】：检查 density 是否为空再合并
  if (length(density) == 0) {
    message(paste("Warning: No data for source", source_id))
    next
  }
  
  # 创建包含分组的数据框
  plot_data <- do.call(rbind, density)
  plot_data = filter(plot_data, Distance < 10000)
  
  # 绘制叠加密度图
  # 注意：确保 niche_col 已经提前定义好了
  p <- ggplot(plot_data, aes(x = Distance, color = Group, fill = Group)) +
    geom_density(alpha = 0, linewidth = 0.8) + 
    labs(title = paste0(source_id, ""),
         x = " Distance to Airway or Artery ", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = niche_col) + 
    scale_color_manual(values = niche_col) + 
    theme(legend.position = "bottom")
  
  print(p) # 在循环中建议显式 print 图片
  
  # 保存
  ggsave(filename = paste0(source_id, "细胞比例随气管距离变化", ".pdf"), 
         path = "fig/叠加密度图", plot = p)
}






















###########################################################################################################
###T细胞亚群的密度图#####################################
knn_cell4=function(data,k=5,typeid=NULL,nihc4="blood"){
  library(FNN)
  # 提取坐标
  if (!is.null(typeid)) {
    cell1_data <-filter(data, type==typeid)
  } else {
    cell1_data <- data
  }
  cell2_data <- filter(data, nich4==nihc4)
  # 确保cell2细胞存在
  if (nrow(cell2_data) == 0) {
    stop("数据中没有", cell2, "细胞，无法计算距离。")
  }
  # 提取坐标
  cell1_coords <- cell1_data[, c("x", "y")]
  cell2_coords <- cell2_data[, c("x", "y")]
  # 计算每个细胞到最近的k个邻居的距离
  nearest_dist <- knnx.dist(
    query = cell1_coords,  # 需要查询的点
    data = cell2_coords,  # 参考点
    k = k # 找最近的k个邻居
  )
  # 将结果添加到inf数据中
  nearest_dist=as.data.frame(nearest_dist)
  cell1_data$nearest_dist <-rowSums(nearest_dist)
  return(cell1_data)
}

source=c("d101","d102","d211","d212","d451","d452")
#metat=subset(meta,type2=="T")
cellt=c("CTL2","TRM8","pro_T","T4",
        "CTL1")
#metat=subset(meta,type2=="T")
cellt=c("DC1","DC2","DC3","IM1","IM2","inf_Mo","MoDC",
        "pro_DC")
for ( source_id in source) {
  print(source_id)
  meta211=filter(meta,source==source_id)
  density=list()
  for (cell_t in cellt) {
    dis=knn_cell4(meta211,nihc4="struct",typeid =cell_t )
    df=data.frame(Group = cell_t, Distance = dis$nearest_dist)
    density[[cell_t]]=df
  }
  # 创建包含分组的数据框
  plot_data <- do.call(rbind, density)
  plot_data=filter(plot_data,Distance<15000)
  # 绘制叠加密度图
  ggplot(plot_data, aes(x = Distance, color = Group, fill = Group)) +
    geom_density(alpha = 0, linewidth = 0.8) +  # 设置透明度和线宽
    labs(title = paste0(source_id,"细胞比例随距离变化"),
         x = "Nearest Distance", y = "Density") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +  # 填充色使用ColorBrewer配色
    scale_color_brewer(palette = "Set2") + # 边框色使用相同配色
    theme(legend.position = "bottom")  # 将图例放在底部
  ggsave(filename =paste0(source_id,"细胞比例随气管距离变3",".png"),path = "合并气管动脉/叠加密度图——T细胞" )
  ####结果明显提升
}
a=knn_nich4(meta211,nihc4="struct",typeid = celltype)
#######################################################################################################################################
###总T细胞的密度图#####################################
knn_cell4=function(data,k=5,typeid=NULL,nihc4="blood"){
  library(FNN)
  # 提取坐标
  if (!is.null(typeid)) {
    cell1_data <-filter(data, type2==typeid)
  } else {
    cell1_data <- data
  }
  cell2_data <- filter(data, nich4==nihc4)
  # 确保cell2细胞存在
  if (nrow(cell2_data) == 0) {
    stop("数据中没有", cell2, "细胞，无法计算距离。")
  }
  # 提取坐标
  cell1_coords <- cell1_data[, c("x", "y")]
  cell2_coords <- cell2_data[, c("x", "y")]
  # 计算每个细胞到最近的k个邻居的距离
  nearest_dist <- knnx.dist(
    query = cell1_coords,  # 需要查询的点
    data = cell2_coords,  # 参考点
    k = k # 找最近的k个邻居
  )
  # 将结果添加到inf数据中
  nearest_dist=as.data.frame(nearest_dist)
  cell1_data$nearest_dist <-rowSums(nearest_dist)
  return(cell1_data)
}

source=c("d101","d102","d211","d212","d451","d452")
#metat=subset(meta,type2=="T")
cellt=c("T")
table(meta$type2)
#metat=subset(meta,type2=="T")
# cellt=c("DC1","DC2","DC3","IM1","IM2","inf_Mo","MoDC",
#         "pro_DC")
for ( source_id in source) {
  print(source_id)
  meta211=filter(meta,source==source_id)
  density=list()
  for (cell_t in cellt) {
    dis=knn_cell4(meta211,nihc4="struct",typeid =cell_t )
    df=data.frame(Group = cell_t, Distance = dis$nearest_dist)
    density[[cell_t]]=df
  }
  # 创建包含分组的数据框
  plot_data <- do.call(rbind, density)
  plot_data=filter(plot_data,Distance<15000)
  # 绘制叠加密度图
  ggplot(plot_data, aes(x = Distance, color = Group, fill = Group)) +
    geom_density(alpha = 0, linewidth = 0.8) +  # 设置透明度和线宽
    labs(title = paste0(source_id,"细胞比例随距离变化"),
         x = "Nearest Distance", y = "Density") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +  # 填充色使用ColorBrewer配色
    scale_color_brewer(palette = "Set2") + # 边框色使用相同配色
    theme(legend.position = "bottom")  # 将图例放在底部
  ggsave(filename =paste0(source_id,"细胞比例随气管距离变3",".png"),path = "合并气管动脉/叠加密度图——T细胞" )
}
############################################################################################3
##基因的空间定位#########################################################
meta211=spa211@meta.data
nich4=meta211$nich4
table(nich4)
nich4=as.character(nich4)
nich4[nich4%in%c("CN7_Airway","CN8_Artery")]="struct"
meta211$nich4=nich4
exp=FetchData(spa211,vars = c("Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9", "Ccl11", "Ccl12", "Ccl17",
                              "Ccl19", "Ccl20", "Ccl21a", "Cxcl9", "Cxcl10", "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl16",
                              "Ifng", "Ifna", "Ifnb", "Il1b", "Il2", "Il4", "Il6", "Il7", "Il10","Il12a", "Il15", "Il17a", "Tgfb1", "Tnf"))
a=knn_nich4(meta211,nihc4="struct")

# 创建结果数据框
result_df <- data.frame(
  cell_id = a$object_id,
  avg_distance = a$nearest_dist
)
head(exp)
head(result_df)
result_df=cbind(result_df,exp)
##############################
#ggplot(result_df,aes(x=avg_distances,y=Cxcl10))+geom_point()
# df=data.frame(dis=result_df$avg_distance,expresion=result_df$Ccl19)
# df=filter(df,dis<10000)
# source=c("d101","d102","d211","d212","d451","d452")
# #metat=subset(meta,type2=="T")
# cellt=c("T")
# table(meta$type2)
#metat=subset(meta,type2=="T")
# cellt=c("DC1","DC2","DC3","IM1","IM2","inf_Mo","MoDC",
#         "pro_DC")
####################################################################################3
#####每个基因循环出图
# 计算每个分箱的统计量
# 设置分箱长度（例如每100单位一个分箱）
bin_width <- 200
for(gene in colnames(exp)){
  df=data.frame(dis=result_df$avg_distance,expresion=result_df[[gene]])
  df=filter(df,dis<10000)
  # 计算最小和最大距离
  min_dist <- min(df$dis)
  max_dist <- max(df$dis)
  
  # 创建分箱序列
  breaks <- seq(from = min_dist, to = max_dist + bin_width, by = bin_width)
  
  # 创建距离分箱
  df$distance_bin <- cut(df$dis, breaks = breaks, include.lowest = TRUE)
  df_binned <- df %>%
    group_by(distance_bin) %>%
    summarise(
      n_cells = n(),  # 细胞数量
      mean_expression = mean(expresion),  # 平均表达量
      expression_rate = mean(expresion > 0),  # 表达细胞比例
      median_expression = median(expresion),  # 中位数表达量
      expression_positive = mean(expresion[expresion > 0])  # 表达细胞的平均表达量
    ) %>%
    # 计算分箱中点距离
    mutate(
      mid_distance = as.numeric(sub("\\((.+),.*", "\\1", distance_bin)) + 
        (as.numeric(sub(".*,(.+)\\]", "\\1", distance_bin)) - 
           as.numeric(sub("\\((.+),.*", "\\1", distance_bin))) / 2
    )
  # 计算最大细胞数量和最大平均表达量
  max_n_cells <- max(df_binned$n_cells, na.rm = TRUE)
  max_mean_expression <- max(df_binned$mean_expression, na.rm = TRUE)
  
  # 创建组合图
  ggplot(df_binned, aes(x = mid_distance)) +
    geom_bar(aes(y = n_cells), 
             stat = "identity", 
             fill = "steelblue", 
             alpha = 0.5) +
    labs(title = "细胞数量与平均表达量随距离的变化",
         x = "到目标结构的距离")+
    geom_line(aes(y = mean_expression * max_n_cells / max_mean_expression), 
              color = "red", 
              size = 1) +
    scale_y_continuous(
      name = "细胞数量",
      sec.axis = sec_axis(~ . * max_mean_expression / max_n_cells, 
                          name = "平均表达量")
    ) +
    labs(title = paste0(gene,"细胞数量与平均表达量随距离的变化"),
         x = "到目标结构的距离") +
    theme_minimal() +
    theme(
      axis.title.y.left = element_text(color = "steelblue", size = 12),
      axis.text.y.left = element_text(color = "steelblue", size = 10),
      axis.title.y.right = element_text(color = "red", size = 12),
      axis.text.y.right = element_text(color = "red", size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12)
    )
  ggsave(filename =paste0(gene,"平均表达量距离脉管区",".pdf"),path = "fig/细胞因子d211" )
}



###################

