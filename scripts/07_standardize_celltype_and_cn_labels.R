library(Seurat)
library(ggplot2)
library(RColorBrewer)
setwd("D:/atlas")
DimPlot(pbmc,label = T,raster=FALSE)
pbmc$o_clu=pbmc@active.ident
pbmc=RenameIdents(pbmc,"PNEC"="Club")
pbmc=RenameIdents(pbmc,"TRM8"="TRM_T8")
pbmc=RenameIdents(pbmc,"CTL2"="Tex_T8")
pbmc=RenameIdents(pbmc,"CTL1"="CTL","T4"="Th")
pbmc=RenameIdents(pbmc,"Trem2"="Trem2_M")
pbmc=RenameIdents(pbmc,"Mo"="cMo")
pbmc=RenameIdents(pbmc,"Mes"="Meso")
DimPlot(pbmc,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(52),raster=FALSE)
##长变量改顺序
lev=levels(pbmc@active.ident)
writeLines(lev, "original_order2.txt")
reordered_vec <- readLines("original_order2.txt")
#file.remove("original_order2.txt")


clu=factor(pbmc@active.ident,levels = reordered_vec)
pbmc$clu=clu

DimPlot(pbmc,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(51),raster=FALSE,group.by = "clu")

ggsave("fig/pbmc_umap.pdf",dpi=300,height = 12,width = 20)
table(spa1$type2,spa1$type)
FeaturePlot(pbmc,features = c("Pecam1","Epcam"))
#################################################
###重命名空转数据###################################################
for (l in 1:8) {
  spa1=spa[[l]]
  spa1[["umap"]] <- spa1[["spatial"]]
  iden=Idents(spa1)
  spa1@active.ident=spa1$type
  spa1=RenameIdents(spa1,"PNEC"="Club")
  spa1=RenameIdents(spa1,"TRM8"="TRM_T8")
  spa1=RenameIdents(spa1,"CTL2"="Tex_T8")
  spa1=RenameIdents(spa1,"CTL1"="CTL","T4"="Th")
  spa1=RenameIdents(spa1,"Trem2"="Trem2_M")
  spa1=RenameIdents(spa1,"Mo"="cMo")
  spa1=RenameIdents(spa1,"Mes"="Meso")
  
  clu=factor(spa1@active.ident,levels = reordered_vec)
  spa1$clu=clu
  spa1$type=spa1$clu
  DimPlot(spa1,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(52),raster=FALSE,group.by = "clu",reduction = "spatial")
  ##
  ##生成细胞大类
  type2=data.frame(ClusterID=unique(spa1@meta.data$type),
                   type2='unkown')  #构建数据框
  
  ## 判断亚群ID属于那类细胞
  type2[type2$ClusterID %in% c("AT1","AT2"),2]='AT'
  type2[type2$ClusterID %in% c("naive_B","GCB","pro_B","Plasma","pre_B","mem_B"),2]='B'
  type2[type2$ClusterID %in% c("Ace_Mo","inf_Mo","cMo"),2]='Mo' 
  type2[type2$ClusterID %in% c("DC1","DC2","DC3","pro_DC","pDC","MoDC"),2]='DC' 
  type2[type2$ClusterID %in% c("AM","inf_AM","pro_AM","IM1","IM2","Trem2_M"),2]='Mac'
  type2[type2$ClusterID %in% c("CTL","Tex_T8","ILC2","ILC3","pro_T","Th","Treg","NK","RDT","TRM_T8","KIR_T8","Naive_T4","Naive_T8"),2]='T'
  #type2[type2$ClusterID %in% c("IM1","IM2"),2]='IM'
  type2[type2$ClusterID %in% c("Mast","Neu"),2]='Gran'
  type2[type2$ClusterID %in% c("Ciliate","Club"),2]='Airway'
  type2[type2$ClusterID %in% c("Endo1","Endo2","Endo3","LEC","ISG_Endo"),2]='Endo'
  type2[type2$ClusterID %in% c("Fat","SMC","Mus","Meso","Pericyte"),2]='Stroma'
  type2[type2$ClusterID %in% c("Fib"),2]='Fib'
  #type2[type2$ClusterID %in% c("ISG_Endo"),2]='ISG_Endo'
  table(type2$ClusterID)
  ## 重新赋值
  spa1@meta.data$type2 = "NA"
  for(i in 1:nrow(type2)){
    spa1@meta.data[which(spa1@meta.data$type == type2$ClusterID[i]),'type2'] <- type2$type2[i]}
  table(spa1@meta.data$type2)
  #        B       DC     Mono Platelet        T 
  #      342       32      642       14     1608 
  spa1$type2=factor(spa1$type2,levels = c("B","T","DC","Mo","Mac","Gran","AT","Airway","Endo","Stroma","Fib"))
  DimPlot(spa1, reduction = 'spatial', group.by = 'nich2',raster=FALSE,
          label = TRUE, pt.size = 0.2)
  table(spa1$nich2,spa1$type2)
  
  ####更换niche名称
  spa1$nich2
  table(spa1$nich2)
  # 1. 定义字典 (如果还没运行过)
  dict <- c(
    "inf_AM"       = "CN1_AM",
    "inf_endo"     = "CN2_ISG",
    "BALT"         = "CN3_iBALT",
    "lym_vessel"   = "CN4_Lymphatic",
    "Plasma"       = "CN5_Plasma",
    "Fib"          = "CN6_AdvFib",
    "trachea"      = "CN7_Airway",
    "Artery"       = "CN8_Artery",
    "SMC"          = "CN9_SMC",
    "Small_vessel" = "CN10_Capillary",
    "stromal"      = "CN11_Alveolar"
  )
  
  # 2. 修正后的赋值代码 (使用 @meta.data 绕过 Seurat 的名称检查)
  spa1@meta.data$nich4 <- factor(dict[as.character(spa1$nich2)], levels = dict)
  #spa1@meta.data$nich4<-(dict[as.character(spa1$nich2)])
  # 3. 检查
  table(spa1$nich4)
  # DimPlot(spa1, reduction = 'spatial', group.by = 'type',raster=FALSE,
  #         label = TRUE, pt.size = 0.2)
  # DimPlot(spa1, reduction = 'spatial', group.by = 'nich4',raster=FALSE,
  #         label = TRUE, pt.size = 0.2)
  spa[[l]]=spa1
}
spa=spa[c(1:8)]
# for (i in 1:8) {
#   spa1=spa[[i]]
#   spa1@meta.data$nich4<- factor(dict[as.character(spa1$nich2)], levels = dict)
#   spa[[i]]=spa1
# }
save(spa,file = "spa.rdata")
library(dplyr)

# 直接修改元数据列
spam@meta.data$nich4 <- recode(spam@meta.data$nich4, "CN3_iBALT" = "CN3_BALCs")
spam@meta.data$nich4 <- recode(spam@meta.data$nich4, "BALCs" = "CN3_BALCs")
# 检查结果
table(spam@meta.data$nich4)
########################################################
meta=spam@meta.data
table(meta$nich4,meta$type2)
save(spa,file="spa.rdata")
