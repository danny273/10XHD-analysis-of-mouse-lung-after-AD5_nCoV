library(dplyr)
library(sf)
library(dbscan)
library(igraph)
library(Seurat)

# --- 1. 定义期望的微环境因子顺序 (CN1 到 CN11) ---
cn_levels <- c("CN1_AM", "CN2_ISG", "CN3_BALCs", "CN4_Lymphatic", "CN5_Plasma", 
               "CN6_AdvFib", "CN7_Airway", "CN8_Artery", "CN9_SMC", 
               "CN10_Capillary", "CN11_Alveolar")

# --- 2. 空间邻域寻找函数 ---
find_nei <- function(data, k=10, eps=50){
  coords <- st_coordinates(st_as_sf(data, coords = c("x", "y")))
  nn <- frNN(coords, eps = eps, sort = FALSE)
  edges <- data.frame(from = rep(1:length(nn$id), lengths(nn$id)), to = unlist(nn$id))
  g <- graph_from_data_frame(edges[edges$from < edges$to, ], directed = FALSE, vertices = data.frame(name = 1:nrow(data)))
  data$pla <- paste0("pla_", igraph::components(g)$membership)
  return(data[data$pla %in% names(which(table(data$pla) > k)), ]) # 直接过滤小群
}

# --- 3. 核心执行流 ---
immune_cells <- c("B", "T", "DC", "Mo", "Mac", "Gran")
plal <- list() 

for (i in seq_along(spa)) {
  obj <- spa[[i]]
  
  # 【关键修复】：备份旧标签为 nich5，将 nich4 转为字符型以备修改
  obj$nich5 <- obj$nich4 
  obj$nich4 <- as.character(obj$nich4)
  target_cells <- character(0)
  
  # 提取 niche16 坐标
  df_niche16 <- data.frame(x = obj$x, y = obj$y, row.names = colnames(obj))[obj$niche == "niche16", ]
  
  if (nrow(df_niche16) > 0) {
    df_pla <- find_nei(df_niche16, k = 10, eps = 50)
    sub_obj <- subset(obj, cells = rownames(df_pla))
    sub_obj$pla <- df_pla$pla
    plal[[i]] <- sub_obj 
    
    # 计算比例并提取目标细胞
    high_plas <- sub_obj@meta.data %>% 
      group_by(pla) %>% 
      summarise(pct = sum(type2 %in% immune_cells) / n() * 100) %>% 
      filter(pct >= 50) %>% 
      pull(pla)
    
    target_cells <- rownames(sub_obj@meta.data)[sub_obj$pla %in% high_plas]
    
    # 【关键修复】：在最新的 nich4 上进行修改，把满足条件的转为 CN3
    if (length(target_cells) > 0) {
      obj$nich4[target_cells] <- "CN3_BALCs"
    }
  }
  
  # 【关键修复】：将 nich4 重新恢复为按照 CN1-CN11 排序的 Factor
  obj$nich4 <- factor(obj$nich4, levels = cn_levels)
  
  # 写回列表
  spa[[i]] <- obj
  cat(sprintf("[%s] 成功将 %d 个高纯度细胞修正为 CN3_BALCs。nich4 已恢复 Factor 排序。\n", names(spa)[i], length(target_cells)))
}
