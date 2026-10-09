
  
library(Seurat)
library(dplyr)
library(RANN)
library(Matrix)
library(future.apply)
library(nichenetr)
library(tidyverse)

lr_network <- readRDS("nichenet/lr_network_mouse_21122021.rds")

unique_lr <- lr_network %>%
  dplyr::select(from, to) %>%
  distinct()

workers_count <- 6
plan("multisession", workers = workers_count)
options(future.globals.maxSize = 7000 * 1024^2)

radius <- 100
max_k <- 100
batch_size <- 100

for(i in seq_along(spa)){
  message("====== 正在处理第 ", i, " 个切片 ======")
  
  spa_i <- spa[[i]]
  DefaultAssay(spa_i) <- "Spatial"
  
  coords <- spa_i@meta.data[, c("x", "y")]
  expr_mat <- GetAssayData(spa_i, assay = "Spatial", slot = "data")
  genes <- rownames(expr_mat)
  
  valid_lr_i <- unique_lr %>%
    filter(from %in% genes, to %in% genes)
  
  message("可计算 L-R 对数：", nrow(valid_lr_i))
  
  nn_res <- nn2(coords, coords, k = max_k, searchtype = "radius", radius = radius)
  
  idx_row <- rep(seq_len(nrow(nn_res$nn.idx)), times = max_k)
  idx_col <- as.vector(nn_res$nn.idx)
  keep <- idx_col > 0
  
  adj_mat <- sparseMatrix(
    i = idx_row[keep],
    j = idx_col[keep],
    x = 1,
    dims = c(nrow(coords), nrow(coords))
  )
  
  rs <- rowSums(adj_mat)
  rs[rs == 0] <- 1
  adj_mat_norm <- adj_mat / rs
  
  ligand_mat <- expr_mat[unique(valid_lr_i$from), , drop = FALSE]
  receptor_mat <- expr_mat[unique(valid_lr_i$to), , drop = FALSE]
  ligand_smooth <- ligand_mat %*% t(adj_mat_norm)
  
  n_pairs <- nrow(valid_lr_i)
  num_batches <- ceiling(n_pairs / batch_size)
  
  lr_list <- future_lapply(seq_len(num_batches), function(b){
    s <- (b - 1) * batch_size + 1
    e <- min(b * batch_size, n_pairs)
    
    L <- valid_lr_i$from[s:e]
    R <- valid_lr_i$to[s:e]
    
    score <- sqrt(
      ligand_smooth[L, , drop = FALSE] *
        receptor_mat[R, , drop = FALSE]
    )
    
    drop0(as(score, "dgCMatrix"))
  }, future.seed = TRUE)
  
  lr_mat <- do.call(rbind, lr_list)
  rownames(lr_mat) <- paste(valid_lr_i$from, valid_lr_i$to, sep = "_")
  colnames(lr_mat) <- colnames(spa_i)
  
  spot_totals <- colSums(lr_mat)
  spot_totals[spot_totals == 0] <- 1
  
  lr_mat_norm <- lr_mat %*% Diagonal(x = 1 / spot_totals)
  rownames(lr_mat_norm) <- rownames(lr_mat)
  colnames(lr_mat_norm) <- colnames(lr_mat)
  
  spa_i[["LR"]] <- CreateAssayObject(data = lr_mat_norm)
  DefaultAssay(spa_i) <- "LR"
  
  spa[[i]] <- spa_i
  
  rm(spa_i, expr_mat, ligand_mat, receptor_mat, ligand_smooth, lr_list, lr_mat, lr_mat_norm)
  gc()
}

plan(sequential)
message("全部切片计算完成，结果已直接存入 LR Assay。")

######################################################################################
lr_markers <- FindAllMarkers(spa1, group.by = "nich4",
                             only.pos = TRUE, 
                             min.pct = 0.01, 
                             logfc.threshold = 0.5,
                             test.use = "wilcox")

print(paste("✅ 找出了", nrow(lr_markers), "个显著富集的 L-R 关系！"))
# 假设 spa211 是你的 Seurat 对象，LR_Radius 是你刚建好的通讯 Assay
DefaultAssay(spa1) <- "LR"

# 提取你空间数据中原本的 RNA counts 总数，作为深度的代理变量
# 回归掉 nCount_Spatial 带来的全局性抬升
spa211 <- ScaleData(spa1, 
                    features = rownames(spa1), 
                    vars.to.regress = "nCount_Spatial")
m21_corrected =FindMarkers(spa1,ident.1 = "CN1_AM",ident.2 = "CN2_ISG",group.by = "nich4")

# 3. 提取 CN2 的特异性特征
CN2_specific <- m21_corrected[m21_corrected$avg_log2FC < 0, ]
# 使用逻辑回归进行差异分析，并回归掉空间 spot 的总 count 数
spatial_LR_markers_LR <- -FindAllMarkers(spa1, group.by = "nich4",
                                         only.pos = TRUE, 
                                         min.pct = 0.01, 
                                         logfc.threshold = 0.5,
                                         test.use = "wilcox")

# 查看结果
top_LR_model <- spatial_LR_markers_LR %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
# ==========================================
# 2. 格式化重构 (完美复刻你发的图片格式)
# ==========================================
print("正在整理数据为宽表格式...")

# A. 提取每个组排名前 20 的 L-R 对 (根据 logFC 排序)
top20_lr <- lr_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 20, order_by = avg_log2FC) %>%
  mutate(
    # 生成排名 (1, 2, 3...)
  rank = row_number(),
    # 替换下划线为连字符，并拼接排名 (例如: "1 Arid3a-A1cf")
    gene_format = gsub("_", "-", gene),
   # formatted_name = paste(rank, gene_format)
  ) %>%
  ungroup()

# B. 将“长表”转换为图片中的“宽表”
wide_lr_table <- top20_lr %>%
  select(rank, cluster,gene_format) %>%
  pivot_wider(
    names_from = cluster, 
    values_from =gene_format
  ) %>%
  arrange(rank)
FeaturePlot(spa211,features = "Tgfbi-Itga4")
pair=rownames(spa211)
pair=as.data.frame(pair)
spa211=spa[[5]]
spam <- merge(
  x = spa[[1]], 
  y = spa[2:8], 
  add.cell.ids = names(spa)
)
DotPlot(spa211,features = c("Tgfbi-Itga4","Ccl5-Ccr5","Ccl19-Ccr7","Ifng-Ifngr1","H2-M2-Klrc1",
                            "H2-Q4-Klrc1","Cxcl10-Cxcr3","Cxcl9-Cxcr3"),group.by = "nich4")
DotPlot(pbmc,features = c("H2-M2","H2-Q4"))       
DotPlot(spam,features = c("Tgfbi-Itga4","Ccl5-Ccr5","Ccl19-Ccr7","Ifng-Ifngr1","H2-M2-Klrc1",
                            "H2-Q4-Klrc1","Cxcl10-Cxcr3","Cxcl9-Cxcr3"),group.by = "source")
DotPlot(spam,features = c("Cd40lg-Cd40","Cd80-Ctla4","Icosl-Icos"),group.by = "source")
DotPlot(spa211,features = c("Tgfbi-Itga4","Ccl5-Ccr5","Ccl19-Ccr7","Ifng-Ifngr1","H2-M2-Klrc1",
                            "H2-Q4-Klrc1","Cxcl10-Cxcr3","Cxcl9-Cxcr3"),group.by = "nich4")
FeaturePlot(pbmc,features = c("F10","F7","F3","Itgb2"))
ma12=FindMarkers(spa1,ident.1 = "CN1_AM",ident.2 = "CN2_ISG",group.by = "nich4")
#########################################################################################
library(ggplot2)
head(names(niche_col),5)

spam2=subset(spam,nich4%in%head(names(niche_col),5))
DotPlot(
  spam, 
  features = c("H2-M2-Klrc1","H2-M2-Cd8a","H2-M2-Klrd1","F10-F7","F7-F3","F10-F3","Cxcl10-Cxcr3","Cxcl9-Cxcr3","Ccl4-Ccr5","Tgfbi-Itga4","Tnf-Tnfrsf1a",
               "Lta-Tnfrsf1b","Cd40lg-Cd40","Icosl-Icos","Cd80-Ctla4","Ccl21a-Ccr7"), 
  group.by = "nich4"
) + 
 # coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 旋转现有的 X 轴标签防止重叠



library(dplyr)
library(tibble)

cn2_idx <- which(spa1$nich4 == "CN2_ISG")

df <- tibble(
  cell = colnames(spa1)[cn2_idx],
  type = spa1$type[cn2_idx],
  LR_Cd274_Pdcd1 = as.numeric(
    GetAssayData(spa1, assay = "LR_Radius", slot = "data")["Cd274-Pdcd1", cn2_idx]
  )
)

df %>%
  filter(type == "TRM_T8") %>%
  summarise(
    n = n(),
    LR_sum = sum(LR_Cd274_Pdcd1),
    LR_mean = mean(LR_Cd274_Pdcd1),
    LR_positive_n = sum(LR_Cd274_Pdcd1 > 0),
    LR_pct_positive = mean(LR_Cd274_Pdcd1 > 0),
    LR_max = max(LR_Cd274_Pdcd1)
  )


for (i in 1:8) {
  spa1=spa[[i]]
  spa1@assays$LR_Radius=NULL
  spa[[i]]=spa1
}
