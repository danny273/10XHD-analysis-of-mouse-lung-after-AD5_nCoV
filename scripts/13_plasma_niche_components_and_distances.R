library(sf)
library(dbscan)
library(igraph)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(rstatix)
library(tidyverse)
library(rstatix)

find_nei=function(data,k=20,eps = 200){
  # 步骤1: 准备数据
  points_sf <- st_as_sf(data, coords = c("x", "y"))
  coords <- st_coordinates(points_sf)
  # 步骤2: 查找半径100内的邻居
  nn <- frNN(coords, eps = eps, sort = FALSE)
  # 步骤3: 构建边列表（修复孤立点问题）
  edges <- data.frame(
    from = rep(1:length(nn$id), lengths(nn$id)),
    to = unlist(nn$id)
  )
  edges <- edges[edges$from < edges$to, ]  # 去重
  # 步骤4: 创建图时显式包含所有顶点
  all_nodes <- data.frame(name = 1:nrow(data))  # 确保所有点作为顶点
  g <- graph_from_data_frame(
    edges, 
    directed = FALSE,
    vertices = all_nodes  # 关键修复：强制包含孤立点
  )
  # 步骤5: 连通分量分析
  clusters <- igraph::components(g)
  data$nichn <- clusters$membership
  print(table(data$nichn))
  # 步骤6: 过滤小于K的小群
  ind=table(data$nichn)
  names(ind)[ind>k]
  data=data[data$nichn%in%names(ind)[ind>k],]
  return(data)
}
add_feature=function(pbmc,fea){
  cell=colnames(pbmc)
  fea2=ifelse(cell%in%names(fea),fea[cell],"other")
  table(fea2)
  return(fea2)
}
##计算cell1到距离最近的k个cell2的距离
knn_cell=function(data,k=5,cell1="TLS",cell2="blood"){
  # 安装并加载FNN包（如果未安装）
  library(FNN)
  # 提取坐标
  cell1_data <- subset(data, cell == cell1)
  cell2_data <- subset(data, cell == cell2)
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
#对id排序，把pla154放到后面
sort_id<- function(data, id_column = "id") {
  data %>%
    mutate(
      time_point = str_extract(.data[[id_column]], "^d\\d+"),  # 提取时间点
      internal_num = as.integer(str_extract(.data[[id_column]], "\\d+$"))  # 提取末尾数字
    ) %>%
    arrange(time_point, internal_num) %>%  # 按时间点和内部数字排序
    select(-time_point, -internal_num)  # 移除临时列
}
# 定义分类函数
classify_immunoglobulin <- function(dat) {
  # 确保必要的列存在
  required_cols <- c("Ighg2c", "Ighg2b", "Ighg1", "Ighg3", "Igha", "Ighm", "Ighd")
  if (!all(required_cols %in% colnames(dat))) {
    stop("数据缺少必要的列名")
  }
  
  # 创建新列存储分类结果
  dat$Type <- NA
  
  # 遍历每一行数据
  for (i in 1:nrow(dat)) {
    # 提取当前行的值
    row <- dat[i, ]
    
    # 计算三类免疫球蛋白的值
    IGg_val <- max(row$Ighg2c, row$Ighg2b, row$Ighg1, row$Ighg3, na.rm = TRUE)
    IGa_val <- row$Igha
    IgM_D_val <- max(row$Ighm, row$Ighd, na.rm = TRUE)
    
    # 获取三类值并命名
    class_values <- c(IgG = IGg_val, IgA = IGa_val, `IgM/D` = IgM_D_val)
    
    # 降序排序并获取前两个值
    sorted_values <- sort(class_values, decreasing = TRUE)
    top1_value <- sorted_values[1]
    top2_value <- sorted_values[2]
    top1_class <- names(sorted_values)[1]
    top2_class <- names(sorted_values)[2]
    
    # 分类逻辑
    if (top1_value == 0) {
      dat$Type[i] <- "未分类"  # 所有表达为0
    } else if (top2_value >= (0.5) * top1_value) {
      dat$Type[i] <- paste(top1_class, top2_class, sep = "&")  # 复合型
    } else {
      dat$Type[i] <- top1_class  # 单一型
    }
  }
  
  return(dat)
}
##分类亚群
merge_type <- function(type) {
  # Define the mapping between cell types and broader categories
  type_mapping <- list(
    AT = c("AT1", "AT2"),
    B = c("B", "GCB", "pro_B", "Plasma", "pre_B"),
    Mo = c("Ace_Mo", "inf_Mo", "Mo", "Trem2"),
    DC = c("DC1", "DC2", "DC3", "pro_DC", "pDC", "MoDC"),
    AM = c("AM", "inf_AM", "pro_AM"),
    T = c("CTL1", "CTL2", "ILC2", "ILC3", "pro_T", "T4", "Treg", "NK", "RDT", "TRM8", "KIR_T8", "Naive_T4", "Naive_T8"),
    IM = c("IM1", "IM2", "Prg4_M"),
    Neu = c("Mast", "Neu"),
    epi = c("Ciliate", "Club"),
    Endo = c("Endo1", "Endo2", "Endo3", "LEC", "Mes", "Pericyte"),
    Mus = c("Fat", "SMC", "Mus"),
    Fib = "Fib",
    ISG_Endo = "ISG_Endo"
  )
  
  lookup <- setNames(rep(names(type_mapping), lengths(type_mapping)), 
                     unlist(type_mapping))
  lookup["unkown"] <- "unkown"  # Default value
  new=lookup[as.character(type)]
  # new=as.factor(new)
  return(new)
}
#分析变量type2和距离变量之间的关系
dist=function(data3){
  # > head(data3
  #        +      )
  # # A tibble: 6 × 15
  # pla    dist_inf_AM dist_inf_endo dist_BALT dist_Fib dist_Small_vessel dist_lym_vessel dist_Artery dist_vein dist_trachea group id            num type2  type2_ordered
  # <chr>        <dbl>         <dbl>     <dbl>    <dbl>             <dbl>           <dbl>       <dbl>     <dbl>        <dbl> <chr> <chr>       <int> <chr>  <ord>        
  #   1 pla_1       20000          5541.     6055.    4793.              956.          14958.      17378.    10460.        8292. d001  d001_pla_1     20 single single       
  # 2 pla_15      20000          1957.     6844.    1111.              536.           9392.      18222.     2374.        2359. d001  d001_pla_15    20 single single       
  # 3 pla_16      20000          4397.     6948.     265.             1272.           8172.      17802.     2361.         310. d001  d001_pla_16    20 single single       
  # 4 pla_17      20000          3621.    10922.    1158.              987.           1681.      20000      8411.        5659. d001  d001_pla_17    13 single single       
  # 5 pla_2        8218.         2082.     4369.    1400.              558.           2409.       8320.     4428.        2820. d001  d001_pla_2     12 single single       
  # 6 pla_27      20000          4132.     2547.    1069.             2390.            512.       8311.      379.         978. d001  d001_pla_27    12 single single   
  # # 将数据转换为长格式
  # 创建长格式数据并计算p值
  plot_data <- data3 %>%
    select(type2,starts_with("dist_")) %>%
    pivot_longer(-type2, names_to = "var") %>%
    group_by(var) %>%
    mutate(
      y_pos = max(value, na.rm = TRUE) * 1.05,  # 计算p值标签位置
      p_label = ifelse(
        wilcox.test(value ~ type2)$p.value < 0.001, "***",
        ifelse(wilcox.test(value ~ type2)$p.value < 0.01, "**",
               ifelse(wilcox.test(value ~ type2)$p.value < 0.05, "*", "ns"))
      )
    )
  
  # 绘制图表
  p1=ggplot(plot_data, aes(type2, value, fill = type2)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, aes(color = type2)) +
    geom_text(aes(x = 1.5, y = y_pos, label = p_label), size = 5, 
              data = . %>% distinct(var, .keep_all = TRUE)) +
    geom_hline(yintercept = 20000, linetype = "dashed", color = "red", alpha = 0.3) +
    facet_wrap(~var, scales = "free_y") +
    labs(x = "Type2", y = "Distance") +
    theme_minimal() +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    scale_fill_manual(values = c("#FFD700", "#87CEFA"))
  return(p1)
}
##分组变量和数值变量间关系
dist_plot <- function(data, num_var, group_var="type2", title = "", xlab = group_var, ylab = num_var) {
  # 检查分组变量是否有两个水平
  if (n_distinct(data[[group_var]]) != 2) stop("分组变量必须包含两个水平")
  
  # 计算p值和标签位置
  p_val <- wilcox.test(reformulate(group_var, num_var), data = data)$p.value
  p_label <- c("ns", "*", "**", "***")[findInterval(p_val, c(0, 0.001, 0.01, 0.05), rightmost.closed = TRUE) + 1]
  y_pos <- max(data[[num_var]], na.rm = TRUE) * 1.05
  
  # 创建标签数据框
  label_df <- data.frame(x = 1.5, y = y_pos, label = p_label)
  
  # 绘制图表
  ggplot(data, aes(x = .data[[group_var]], y = .data[[num_var]], fill = .data[[group_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2, aes(color = .data[[group_var]])) +
    annotate("text", x = 1.5, y = y_pos, label = p_label,  # 关键修复
             size = 6, color = "black") + 
    labs(title = title, x = xlab, y = ylab, 
         caption = "显著性: *** p<0.001; ** p<0.01; * p<0.05; ns 不显著") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.major.x = element_blank()) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    scale_fill_manual(values = c("#FFD700", "#87CEFA"))
}
###多组变量
dist_plot2 <- function(data, num_var, group_var = "clu", 
                      title = "", xlab = group_var, ylab = num_var) {
  
  # 检查分组变量水平数
  n_groups <- n_distinct(data[[group_var]])
  if (n_groups < 2) stop("分组变量至少需要两个水平")
  
  # 根据分组数量选择检验方法
  if (n_groups == 2) {
    test_result <- wilcox.test(reformulate(group_var, num_var), data = data)
  } else {
    test_result <- kruskal.test(reformulate(group_var, num_var), data = data)
  }
  
  # 获取p值并转换为标签
  p_val <- test_result$p.value
  p_label <- c("ns", "*", "**", "***")[findInterval(p_val, c(0, 0.001, 0.01, 0.05), 
                                                    rightmost.closed = TRUE) + 1]
  
  # 计算标签位置
  y_pos <- max(data[[num_var]], na.rm = TRUE) * 1.05
  x_pos <- ifelse(n_groups == 2, 1.5, (n_groups + 1) / 2)
  
  # 创建基础图表
  p <- ggplot(data, aes(x = .data[[group_var]], y = .data[[num_var]], 
                        fill = .data[[group_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2, aes(color = .data[[group_var]])) +
    labs(title = title, x = xlab, y = ylab,
         caption = ifelse(n_groups == 2,
                          "Wilcoxon秩和检验显著性",
                          "Kruskal-Wallis检验显著性")) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.major.x = element_blank())
  
  # 添加显著性标签
  p <- p + annotate("text", x = x_pos, y = y_pos, label = p_label,
                    size = 6, color = "black")
  
  # 动态生成颜色 (处理任意分组数量)
  if (n_groups <= 9) {
    # 使用RColorBrewer的Set1调色板 (支持最多9组)
    p <- p + 
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1")
  } else {
    # 超过9组时使用viridis调色板
    p <- p + 
      scale_fill_viridis_d() +
      scale_color_viridis_d()
  }
  
  return(p)
}

##plasma#################################################################################
plal=list()
for (i in 1:length(spa)) {
  spa4=spa[[i]]
  l=names(spa)[i]
  data=data.frame(spa4$x,spa4$y)
  colnames(data)=c("x","y")
  data$cell=spa4$nich2
  table(data$cell)
  head(data,100)
  pla=data[data$cell=="Plasma",]
  pla=find_nei(pla,k=10,eps = 50)
  pla$nichn=paste0("pla_",pla$nichn)
  ##写入spa
  fea=pla$nichn
  names(fea)=row.names(pla)
  
  spa4$pla=add_feature(spa4,fea)
  table(spa4$pla)
  spa4$pl=paste0(l,spa4$pla)
  length(table(spa4$pla))
  cor = c("grey80",colorRampPalette(brewer.pal(8, "Set1"))(length(table(spa4$pla))-1))
  DimPlot(spa4, group.by = "pla", label = TRUE, repel = TRUE,reduction = "spatial",
          cols =cor,raster=F,pt.size=0.3) + ggtitle("pla clusters") +NoLegend()
  ggsave(paste0(l,"_plaumap.png"), width = 20, height = 20,path = "pla")
  pla=subset(spa4, pla!="other")
  Idents(pla) <- pla$pla
  DotPlot(pla,features = c("Jchain","Ighm","Igha","Ighg2c","Ighg2b","Ighg1","Ighg3","Ighd"),
          cols = c("lightgrey", "blue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("pla clusters") 
  ggsave(paste0(l,"_pladot.png"), width = 10, height = 10,path = "pla")
  plal[[i]]=pla
}

##计算平均表达量
plal_avg <- lapply(plal, function(x) {
  AverageExpression(x, slot = "data")
})
## 合并所有平均表达量数据框####################################################################################
df =plal_avg[[1]][["Spatial"]]
df=as.data.frame(df)
group=rep(names(spa)[1], each = ncol(df))
for (i in 2:length(plal_avg)) {
  a=plal_avg[[i]][["Spatial"]]
  a=as.data.frame(a)
  gene=intersect(rownames(df),rownames(a))
  a=a[gene,]
  df=df[gene,]
  df=cbind(df,a)
  group=c(group,rep(names(spa)[i], each = ncol(a)))
}
df=log(df + 1)  # 对数据进行log10转换
dat=df[c("Jchain","Ighm","Igha","Ighg2c","Ighg2b","Ighg1","Ighg3","Ighd"),]
dat=as.data.frame(t(dat))
dat$group=group                  
plot(density(dat$Ighg2c))
library(Rtsne)
rownames(dat)=paste0(dat$group,"_",rownames(dat))
dat=dat[,-ncol(dat)]  # 删除group列
#dat=as.data.frame(t(dat))  # 转置数据框
head(dat, 10)  # 查看前10行数据
# 检查是否有完全重复的行
duplicated_rows <- duplicated(dat) | duplicated(dat, fromLast = TRUE)
if (any(duplicated_rows)) {
  cat("发现重复行，共", sum(duplicated_rows), "行\n")
  dat_unique <- dat[!duplicated(dat), ]  # 删除重复行
} else {
  dat_unique <- dat  # 如果没有重复行，直接使用原数据
}

tsne_out = Rtsne(
  dat_unique,
  dims = 2,
  pca = T,
  max_iter = 1000,
  theta = 0.4,
  perplexity = 20,
  verbose = F
) # 进行t-SNE降维分析
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
##
tsne_result$group=group[!duplicated(dat)]
tsne_result=cbind(tsne_result,dat_unique)

colnames(tsne_result)
ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(color=group)) +
#  ggtext(tsne_result,label = "class8")+
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )
ggsave("pla_tsne.png", width = 10, height = 10,path = "pla")
for (i in 4:ncol(tsne_result)) {
  ggplot(tsne_result, aes_string(x = "tSNE1", y = "tSNE2")) +
    geom_point(aes_string(color = colnames(tsne_result)[i]), size = 2) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = paste("Expression of", colnames(tsne_result)[i])) +
    theme_minimal() +
    theme(legend.position = "right") 
 ggsave(paste0("pla_tsne_", colnames(tsne_result)[i], ".png"), width = 10, height = 10,path = "pla")
}

##构建seurat对象###########################################################################################
colnames(df)=paste0(group, "_", colnames(df))
pbmc=CreateSeuratObject(df)
##
pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc=ScaleData(pbmc, features = rownames(pbmc))
pbmc=FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=RunPCA(pbmc, features = VariableFeatures(pbmc))
ElbowPlot(pbmc, ndims = 50) # 查看肘部图，选择合适的主成分数
pbmc=FindNeighbors(pbmc,dims = 1:10)
pbmc=FindClusters(pbmc, resolution = 2) # 调整分辨率以获得更好的聚类结果
pbmc=RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP of Plasma Cells")|
  DimPlot(pbmc, reduction = "umap", group.by = "clu", label = TRUE, repel = TRUE) + ggtitle("UMAP of Plasma Cells by Group")
pbmc$clu2=Idents(pbmc)

FeaturePlot(pbmc,features = c("Cfd"))
ma=FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
#top10 <- ma %>% group_by(cluster) %>% slice_min(p_val_adj, n = 10)
top10 <- ma %>% group_by(cluster) %>% slice_min(p_val_adj, n = 10)
top10 <- ma %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC)
#top10 %>% write.csv(file="7_Cell_Markers_smy.csv")
top10_genes <- top10 %>% select(gene)
markers=unstack(top10_genes[,c(2,1)])
FeaturePlot(pbmc,features = c("Jchain","Ighm","Igha","Ighg2c","Ighg2b","Ighg1","Ighg3","Ighd"))
FeaturePlot(pbmc,features = c("Cfd","Fabp4","Car3","Thrsp"))
pbmc$group=group
DimPlot(pbmc,group.by = "group")
colnames(subset(pbmc,idents = "2"))
new.cluster.ids <- c("other", "other", "artery_p", "secrete_p", "other", "trechea_p", "other")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = 'umap', label = TRUE, pt.size =3) + NoLegend()
###
# 提取 UMAP 坐标
umap_coords <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
umap_coords$cell_id <- rownames(umap_coords)  # 添加细胞 ID（barcode）
colnames(umap_coords) <- c("UMAP_1", "UMAP_2", "cell_id")  # 重命名列
library(ggplot2)
exp=df[c("Cfd","Fabp4","Car3","Thrsp"),]
exp=as.data.frame(t(exp))
ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, label = cell_id)) +
  geom_point(size = 1, color = "blue", alpha = 0.6) +  # 绘制点
  geom_text(size = 2, check_overlap = TRUE, nudge_y = 0.1) +  # 标注文本（避免重叠）
  theme_minimal() +
  labs(title = "UMAP 可视化（带细胞 ID 标注）")
ggsave("umap_with_labels.png", width = 10, height = 10, path = "pla")
placlu=pbmc@active.ident
#########参数获取
##细胞数量
num=list()
for (i in 1:length(spa)) {
  a=plal[[i]]
  a1=table(a$pla)
  a1=as.data.frame(a1)
  a1$group=names(spa)[i]
  a1$id=paste0(a1$group,"_",a1$Var1)
  num[[i]]=a1
}
num=bind_rows(num)

##空间定位分析############################################################################################
ldist=list()
for (i in 1:length(spa)) {
  spa4=spa[[i]]
  l=names(spa)[i]
  data=data.frame(spa4$x,spa4$y)
  colnames(data)=c("x","y")
  data$cell=spa4$nich3
  table(spa4$nich3)
  pla=plal[[i]]
  data$pla=spa4$pla
  table(data$cell)
  ######血管和气管
  pla2=subset(data, cell=="Plasma")
  for (m in c("inf_AM","inf_endo",  "BALT", "Fib", 
              "Small_vessel","lym_vessel","Artery","vein",
              "trachea")) {
  d1=knn_cell(data,5,"Plasma",m)$nearest_dist
  pla2=cbind(pla2,d1)
  colnames(pla2)[ncol(pla2)] <- paste0("dist_", m)
  }
  pla2=pla2[rownames(pla2) %in% colnames(pla), ]
  pla2$pla=pla$pla
  head(pla2,10)
  dist_columns <- grep("^dist_", names(pla2), value = TRUE)  # 提取所有 dist_* 列名
  pla_dis <- pla2 %>%
    group_by(pla) %>%                          # 按 pla 分组
    summarise(across(all_of(dist_columns), mean, na.rm = TRUE))  # 计算所有 dist_* 列的平均值  
  pla_dis$group <- l  # 添加 group 列
  pla_dis$id= paste0(pla_dis$group, "_", pla_dis$pla)
  ldist[[i]] <- pla_dis  # 将结果存储到列表中
} 
ldist <- bind_rows(ldist)  # 合并所有数据框
##汇总
head(num)
ldist$num=num$Freq
head(placlu)
head(ldist)
head(dat)
placlu_df <- placlu %>%
  enframe(name = "id", value = "clu") %>%    # 转换为两列：id 和 clu
  mutate(id = gsub("-", "_", id))            # 将 "-" 替换为 "_"（如 d001_pla-1 → d001_pla_1）
dat_df <- dat %>%
  rownames_to_column(var = "id") %>%
  mutate(
    id = gsub("\\.", "_", id)  # 将第一个点替换为下划线（如 d001_pla.1 → d001_pla_1）
  )
dat_df <- dat_df %>%
  separate(id, into = c("part1", "part2", "part3"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(id = ifelse(!is.na(part3), paste(part1, part2, sub("_.*", "", part3), sep = "_"), id)) %>%
  select(-part1, -part2, -part3)
data <- ldist %>%
  left_join(placlu_df, by = "id") %>%      # 合并 placlu 数据
  left_join(dat_df, by = "id")             # 合并 dat 数据
ldist
#data=data[data$id%in%id,]###排除不合格后补的
table(data$clu,data$group)
##定位量化分析##########################################################library(tidyverse)

# 提取数值变量列名（排除非数值列）
num_vars <- names(data)[sapply(data, is.numeric)]  # 如 dist_inf_AM, Jchain, Ighg2c, num 等

# 执行 Kruskal-Wallis 检验（所有数值变量 vs clu）
kruskal_results <- data %>%
  select(all_of(num_vars), clu) %>%
  gather(key = "variable", value = "value", -clu) %>%
  group_by(variable) %>%
  kruskal_test(value ~ clu) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

library(corrplot)

# 计算相关系数矩阵（Spearman）
cor_matrix <- data %>%
  select(all_of(num_vars)) %>%
  cor(method = "spearman", use = "complete.obs")

# 可视化热图
corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  col = colorRampPalette(c("blue", "white", "red"))(200),  # 蓝-白-红色阶
  tl.col = "black",
  tl.cex = 0.7,
  addCoef.col = "black",
  number.cex = 0.6,
  order = "hclust",          # 使用层次聚类
  hclust.method = "average", # 可选，默认是"complete"
  addrect = 3                # 绘制3个矩形框，即分为3类，这个数字可以根据需要调整
)

################
# 获取所有距离变量列名
data=sort_id(data, id_column = "id")  # 使用自定义函数对data进行排序
dist_cols <- names(data)[grepl("^dist_", names(data))]
# 对所有距离变量应用截断处理
data_plot <- data
data_plot[dist_cols] <- map(data_plot[dist_cols], ~ if_else(.x > 20000, 20000, .x))
ggplot(data_plot, aes(x = dist_vein, y = dist_lym_vessel)) +
  geom_point(aes(color = clu), size = 2) +
 # geom_text(aes(label = id), size = 2, check_overlap = TRUE, nudge_y = 0.1) +  # 添加ID标签
  labs(title = "Distance to vein vs Lymph Vessel", 
       x = "Distance to vein", 
       y = "Distance to Lymph Vessel") +
  scale_x_continuous(limits = c(0, 20000)) +
  scale_y_continuous(limits = c(0, 20000)) +
  theme_minimal() +
  theme(legend.position = "right")
##
ggplot(data_plot, aes(x = dist_trachea, y = dist_Artery)) +
  geom_point(aes(color = clu), size = 2) +
  geom_text(aes(label = id), size = 2, check_overlap = TRUE, nudge_y = 0.1) +  # 添加ID标签
  labs(title = "Distance to BALT vs Lymph Vessel", 
       x = "Distance totrachea", 
       y = "Distance to Artery") +
  scale_x_continuous(limits = c(0, 20000)) +
  scale_y_continuous(limits = c(0, 20000)) +
  theme_minimal() +
  theme(legend.position = "right")
##IG分类#######################################################################################

# 应用分类函数到您的数据
data2 <- classify_immunoglobulin(data_plot)

table(data2$Type, data2$group)  # 查看分类结果与聚类的关系
table(data2$group)
group2=substr(data2$group, 1, 3)  # 去掉前缀
group2[group2=="d01"]="d00"
data2$group2=group2
table(data2$clu,data2$group2)
##抗体种类关系分析######################################################################################
 data2=data2[!data2$Type=="未分类",]
data2$type2=ifelse(data2$Type %in% c("IgG", "IgA", "IgM/D"), "single", "complex")
table(data2$type2, data2$group2)
table(data2$Type, data2$group2)
head(data3,20)
data3=data2[,c(1:13)]
data3$type2=data2$type2
table(data2$type2, data2$clu)
table(data2$group2, data2$clu2)
library(ggplot2)
library(tidyr)
# 将type2转为有序因子
data3$type2_ordered <- factor(data3$type2, levels = c("single", "complex"), ordered = TRUE)

# 计算Spearman相关系数
cor_results <- sapply(
  grep("^dist_", names(data3), value = TRUE),
  function(var) {
    cor.test(
      x = as.numeric(data3$type2_ordered),
      y = data3[[var]],
      method = "spearman"
    )$estimate
  }
)

data4=data2[!data2$group2%in%c("d00"),]
dist(data2) 
table(data2$Type, data2$group2)
table(data4$type2, data4$clu)
table(data4$type2)
table(data2$type2, data2$group2)
plot_data <- data3 %>%
  select(type2, n,starts_with("dist_")) 

dist_plot(data2,num_var = "num")|dist_plot(data4,num_var = "num")
ma=FindMarkers(pbmc,ident.1 = "2",ident.2 = "5")
VlnPlot(pbmc,features = c("Itga1","Ace","Id2"))
ma=FindMarkers(pbmc,ident.1 = "2",ident.2 = "5")
##差异基因###############################################################################
cell=colnames(pbmc)
cell=gsub("-","_",cell)
ind=cell%in%data2$id
pbmc=pbmc[,ind]
pbmc
####重新排序
data2=sort_id(data2, id_column = "id")  # 使用自定义函数对data2进行排序
pbmc@meta.data=cbind(pbmc@meta.data,data2)
pbmc$type2=as.factor(pbmc$type2)
table(pbmc$type2)
pbmc2=subset(pbmc,subset = group2!="d00")
VlnPlot(pbmc,features = c("Ccl19","Tnfsf13b","Tnfsf13","Ccl21a","Cxcl12","S1pr1"),group.by = "type2",split.by = "group2")
VlnPlot(pbmc,features = c("Cxcl10","Cxcl9","Cxcr3","Ccl28","Aicda"),group.by = "type2",split.by = "group2")
VlnPlot(pbmc,features = c("Itga10","Itgb1","Vcam1"),group.by = "type2",split.by = "group2")
B=subset(spam, subset = type2 == "B")
DotPlot(spam,features =  c("Ccl19","Tnfsf13b","Tnfsf13","Ccl21a","Cxcl12","S1pr1","Vcam1","Mki67"),group.by = "nich3")
DotPlot(B,features =  c("Ccl19","Tnfsf13b","Tnfsf13","Ccl21a","Cxcl12","S1pr1","Vcam1","Mki67"),group.by = "nich3")
ma=FindMarkers(pbmc2,group.by = "type2", ident.1 = "single", ident.2 = "complex", logfc.threshold = 0.25, min.pct = 0.1)
ma <- FindMarkers(pbmc, group.by = "type2", ident.1 = "single", ident.2 = "complex", logfc.threshold = 0.25, min.pct = 0.1, assay = "RNA")
DimPlot(pbmc,group.by = "clu",label = T)
data2$clu2=Idents(pbmc)
table(data2$clu2, data2$type2)
table(data2$clu2, data2$Type)
table(data2$clu2, data2$group2)
###提取细胞构成#######################################################################
cellt=list()
for (i in 1:length(plal)) {
  pla=plal[[i]]
  l=names(spa)[i]
  #cell_type=merge_type(pla$type)#用合并群
  cell_type=pla$type #直接使用原始细胞类型
  names(cell_type) <- colnames(pla)
  pla$type2=cell_type
  a=table(pla$pla,pla$type2)
  a=as.data.frame(a)
  
  a1 <- a %>%
    pivot_wider(
      names_from = Var2,   # 细胞类型作为列名
      values_from = Freq,  # 计数值作为单元格值
      values_fill = 0      # 缺失值用0填充
    ) %>%
    as.data.frame()        # 转换为数据框
  
  rownames(a1) <- a1$Var1
  a1<- a1[, -1]  # 移除Var1列
  rownames(a1) <- paste0(l, "_", rownames(a1)) 
  cellt[[i]] <- a1  # 将结果存储到列表中
}
merged <- cellt %>%
  map_dfr(~ as_tibble(., rownames = "Sample"), .id = "Source") %>%
  mutate(across(-c(Sample, Source), ~ replace_na(., 0))) %>%
  select(Sample, Source, sort(setdiff(names(.), c("Sample", "Source")))) %>%
  arrange(Sample)

cellt <- merged %>% 
  select(-Source) %>% 
  column_to_rownames("Sample") %>% 
  as.data.frame()
head(cellt, 10)
cellt%>%mutate(id= row.names(.)) %>%
  sort_id(id_column = "id") %>%
  select(-id) -> cellt
#cellt=cellt[cellt$B>0,]  # 只保留B细胞类型大于0的行
##回到pbmc剔除不合格细胞######
# cell=colnames(pbmc)
# cell=gsub("-","_",cell)
# ind=cell%in%rownames(cellt)
# table(ind)
# pbmc=pbmc[,ind]
# pbmc
# data2=sort_id(data2, id_column = "id")  # 使用自定义函数对data2进行排序
# data2=data2[data2$id%in%rownames(cellt),]  # 确保data2的id在pbmc中存在
# table(data2$clu,data2$group2)
####################################################
#ind=rownames(cellt) %in% data2$id
#cellt=cellt[ind,]
cellt$sum=rowSums(cellt)  # 计算每行的总和
cellt/cellt$sum -> cellt  # 计算每个细胞类型的比例
cellt=cellt[,-ncol(cellt)]  # 删除sum列
cellt=signif(cellt, digits = 3)  # 保留三位有效数字
cellt=cellt*100
#data2=cbind(data2, cellt)  # 将细胞构成数据合并到data2中
cellt=cellt[rownames(cellt) %in% data2$id, ]  # 确保cellt的行名在data2的id中
cellt$type2=data2$type2
rownames(cellt)== data2$id
dist_plot(data2,num_var = "num", group_var = "type2")
ind=data2$group2 %in% c("d00","d01")
cellt2=cellt[!ind,]
dist_plot(cellt, num_var = "Fat", group_var = "type2", title = "Cell Type Composition by Type2")
dist_plot2(cellt, num_var = "Mus", group_var = "type2", title = "Cell Type Composition by Type2")
plot_data <-cellt2 %>%
  pivot_longer(-type2, names_to = "var") %>%
  group_by(var) %>%
  mutate(
    y_pos = max(value, na.rm = TRUE) * 1.05,  # 计算p值标签位置
    p_label = ifelse(
      wilcox.test(value ~ type2)$p.value < 0.001, "***",
      ifelse(wilcox.test(value ~ type2)$p.value < 0.01, "**",
             ifelse(wilcox.test(value ~ type2)$p.value < 0.05, "*", "ns"))
    )
  )

# 绘制图表
ggplot(plot_data, aes(type2, value, fill = type2)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, aes(color = type2)) +
  geom_text(aes(x = 1.5, y = y_pos, label = p_label), size = 5, 
            data = . %>% distinct(var, .keep_all = TRUE)) +
  geom_hline(yintercept =80, linetype = "dashed", color = "red", alpha = 0.3) +
  facet_wrap(~var, scales = "free_y") +
  labs(x = "Type2", y = "Distance") +
  theme_minimal() +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#FFD700", "#87CEFA"))
###分析细胞比例与clu关系############################
cellt$type2=NULL
cellt$clu=data2$clu
dist_plot2(cellt, num_var = "Mus", group_var = "clu",
            title = "Cell Type Composition by Clu", xlab = "Clu", ylab = "DC")
library(tidyverse)
cellt=cellt[rownames(cellt) %in% data2$id, ]  # 确保cellt的行名在data2的id中
cellt$clu=data2$clu
for (i in 1:50) {
  var= colnames(cellt)[i]
  dist_plot2(cellt, num_var = var, group_var = "clu", title = "Cell Type Composition by Type2")
  ggsave(paste0("cell_type_composition_", var, ".png"), width = 10, height = 6, path = "pla/clu")
}

##
# id=rownames(cellt[cellt$B>0,])
# pbmc$id
# id2=pbmc$id
# id2=gsub("-","_",id2)
# id2
# ind=id2 %in% id
# table(ind)
# pbmc=pbmc[,ind]
# DimPlot(pbmc)
####整理剔除###############################################################################
save(pbmc, file = "D:/atlas/pla/bulk_pbmc.rdata")


for (i in 1:length(plal)) {
  pla=plal[[i]]
  l=names(spa)[i]
  data=data2[data2$group==l, ]
  id=data$pla
  id
  ind=pla$pla%in% id
  pla=pla[,ind]
  plal[[i]] <- pla  # 更新plal列表中的元素
}
save(data2,cellt,plal, file = "D:/atlas/pla/data.rdata")

Qpla_id=c()
source=c()
cell_id=c()
for (i in 1:length(plal)) {
  pla=plal[[i]]
  l=names(spa)[i]
  cell_id=c(cell_id, paste0(l, "_", pla$object_id))
  source=c(source, rep(l, ncol(pla)))
  pla_id=c(pla_id, pla$pla)
}

#load("D:/atlas/spam.Rdata")
cell_id2=colnames(spam)
ind=pla_id %in% cell_id2
df=data.frame(pla_id,source,cell_id)
df$id=paste0(df$source, "_", df$pla_id)
df
df <- merge(df, data2, by = "id", all.x = TRUE)
data2
head(df, 10)
write.csv(df, file = "D:/atlas/pla/pla_id.csv")
d211=spa[[5]]
leiden=read.csv("D:/atlas/leiden.csv", header = TRUE, stringsAsFactors = FALSE)
##做图##############################################################################
DimPlot(pbmc,label = T)
ma=FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
table(data2$clu)
data5=filter(data2,clu=="secrete_p")
data5$id

data2$clu=data$clu
# 获取所有距离变量列名
dist_cols <- names(data)[grepl("^dist_", names(data))]
# 对所有距离变量应用截断处理
data_plot <- data2
data_plot[dist_cols] <- map(data_plot[dist_cols], ~ if_else(.x > 20000, 20000, .x))
ggplot(data_plot, aes(x = dist_Artery, y = dist_trachea)) +
  geom_point(aes(color = clu), size = 2) +
 # geom_text(aes(label = id), size = 2, check_overlap = TRUE, nudge_y = 0.1) +  # 添加ID标签
  labs(title = "Distance to Artery vs  trachea", 
       x = "Distance to Artery", 
       y = "Distance trachea") +
  scale_x_continuous(limits = c(0, 20000)) +
  scale_y_continuous(limits = c(0, 20000)) +
  theme_minimal() +
  theme(legend.position = "right")
meta=pbmc@meta.data
meta=sort_id(meta, id_column = "id")  # 使用自定义函数对meta进行排序
data2$clu2=meta$clu
data4=read.csv("D:/atlas/pla/pla_id.csv")
data4$clu==data2$clu
data3=data2
data2=data3
DotPlot(spam,features = c("Cfd","Fabp4","Car3","Thrsp"),group.by = "type")
table(data2$type2,data2$clu)
table(data2$Type,data2$clu)        
dist_plot2(data2, num_var = "num", group_var = "group2",
           title = "Cell Type Composition by Clu", xlab = "Clu", ylab = "DC")
