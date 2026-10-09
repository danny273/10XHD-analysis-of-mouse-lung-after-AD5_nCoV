library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library("scviR")
library(reticulate)
# BiocManager::install("SeuratWrappers")
# devtools::install_github("satijalab/seurat-wrappers", ref = "seurat5", force = TRUE, upgrade = "never")
# devtools::install_local("D:/软件/seurat-wrappers-master.zip")
# for (i in seq_along(dat3)) {
#   # 设置RNA assay的Key（通常自动生成）
#   pbmc=dat3[[i]]
#   mat=dat3[[i]][["TAG"]]@counts
#   # 设置TAG assay的Key
#   pbmc@assays[["TAG"]]=CreateSeuratObject(counts = mat, assay = "TAG")
#   dat3[[i]]=pbmc
# }
# dat4=list()
# dat5=list()
# for (i in seq_along(dat3)) {
#   a=dat3[[i]]
#   a@assays$TAG=NULL
#  dat4[[i]]=a
#  b=dat3[[i]]
#  b@assays$RNA=NULL
#  dat5[[i]]=b
# }
sce.big <- merge(dat3[[1]], 
                 y = dat3[c(2:length(dat3))],
                 add.cell.ids = file, 
                 #collapse =T,
                 project = "pbmc3")
table(sce.big$tag)
sum(sce.big$S)
# sce.big2 <- merge(dat5[[1]],assay = "TAG",
#                  y = dat5[c(2:length(dat5))],
#                  add.cell.ids = file, 
#                  #collapse =T,
#                  project = "pbmc3")
#sce.big@assays[["TAG"]] <-sce.big2@assays[["TAG"]]
rm(dat,dat2,dat3,d,pbmc)
gc()
sce.big
sce.big[["percent.Mt"]] <- PercentageFeatureSet(sce.big, pattern = "^Mt")
max(sce.big[["percent.Mt"]])
#sce.big[["percent.Ercc"]] <- PercentageFeatureSet(sce.big, pattern = "^Ercc")
sce.big[["percent.Ribosome"]] <- PercentageFeatureSet(sce.big, pattern = "^RP[SL]")
VlnPlot(sce.big, features = c("nFeature_RNA", "Contamination",
                           "percent.Mt",  "scDblFinder.score"), ncol = 4,pt.size=0,group.by = "tag")
######
# scb=FetchData(sce.big,vars = c("Contamination","scDblFinder.score"))
# scb$tag=sce.big$tag
# table(scb$tag)
# mul=scb[scb$tag=="Tag1",]
# ggplot(scb,aes(x=Contamination,y=scDblFinder.score,color=tag))+geom_point()+
#   theme_bw()+xlab("Contamination")+ylab("scDblFinder.score")+theme(legend.position = "none")
# ggplot(mul,aes(x=Contamination,y=scDblFinder.score,color=tag))+geom_point()


table(sce.big[["nFeature_RNA"]]<450)
sce.big<- subset(sce.big, subset = percent.Mt<2&nFeature_RNA>450&scDblFinder.score<0.75&Contamination<0.75)
sce.big<- NormalizeData(sce.big)
sce.big<- FindVariableFeatures(sce.big)
sce.big<- ScaleData(sce.big)
sce.big<- RunPCA(sce.big,features = VariableFeatures(object = sce.big))


# #sce.big<- IntegrateLayers(object = sce.big, method = scVIIntegration, new.reduction = "integrated.scvi", conda_env = "D:/Miniconda3/", verbose = F)
# sce.big<- IntegrateLayers(object = sce.big, method = RPCAIntegration,
#                           orig.reduction = "pca", new.reduction = "integrated.rpca",
#                           verbose = FALSE)
ElbowPlot(sce.big,100)
#################################
# sce.big <- FindNeighbors(sce.big, reduction = "integrated.rpca", dims = 1:32)
# sce.big  <- FindClusters(sce.big , resolution = 0.2)
# sce.big  <- RunUMAP(sce.big , reduction = "integrated.rpca", dims = 1:30)
sce.big <- FindNeighbors(sce.big, reduction = "pca", dims = 1:32)
sce.big  <- FindClusters(sce.big , resolution = 0.2)
sce.big  <- RunUMAP(sce.big , reduction = "pca", dims = 1:32)
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
pbmc=sce.big
pbmc$chip=substr(colnames(pbmc),0,2)
table(pbmc$chip)
DimPlot(pbmc,reduction = "umap",label = TRUE,pt.size = 0.5,split.by = "chip")
DimPlot(pbmc,reduction = "umap",label = TRUE,pt.size = 0.5,split.by = "tag",raster=FALSE,ncol = 4)
FeaturePlot(pbmc,features = c("Cd3e","Cd4","Cd8a","Nkg7"))
FeaturePlot(pbmc,features = "Top2a")
DimPlot(pbmc,split.by = "time")

#####################################################################################
cell=colnames(pbmc)
cell=str_split(cell,"_")

cell=t(cell)
cell=as.data.frame(cell)
pbmc$time=cell$V1
pbmc$id=cell$V2
group=cell$V3
group[group=="tag1"]="iliac"
group[group=="tag3"]="med"
group[group=="tag2"]="cer"
group[group=="tag4"]="lung"
group[group=="tag5"]="muscle"
group[group=="tag7"]="nose"
table(group)
pbmc$group=group
DimPlot(pbmc,split.by = "group",label = T)
sam=paste0(pbmc$time,"_",pbmc$group)
table(sam)
pbmc$sam=sam
DimPlot(pbmc,split.by = "sam",ncol = 6,label = T)
save(pbmc,file = "node1.rdata")
# ####################################################################################
# DefaultAssay(pbmc)="ADT"
# pbmc<- IntegrateLayers(object = sce.big, method = RPCAIntegration,
#                           orig.reduction = "pca", new.reduction = "integrated.rpca",
#                           verbose = FALSE)
 pbmc=JoinLayers(pbmc)
# sam=paste0(pbmc$chip,"_",pbmc$tag)
# pbmc$sam=sam
for (i in unique(sam)) {
  sub=subset(pbmc,subset = sam==i)
  feaplot=FeaturePlot(sub,features = c("Cxcl9","Trem2","Cd8a","Pdcd1"))
  ggsave(paste0(i,"_adt.png"),width = 10,height = 8,path = "fig")
}



################################################################################
# #############样品相关性
# DefaultAssay(pbmc)="RNA"
# mdata=AverageExpression(pbmc,group.by = "sam")
# mdata=mdata[["RNA"]]
# mdata=as.data.frame(mdata)
# d=cor(mdata)
# pheatmap::pheatmap(d)
# save(d,file = "sam_rela.csv")
# indcer=grep("d0_cer",sam)
# indmed=grep("d0_med",sam)
# sam[indcer]="d0_med"
# table(sam)
# sam[indmed]="d0_cer"
# pbmc$sam=sam
# pbmc$group[indcer]="med"
# pbmc$group[indmed]="cer"
# ##################################################################################
# time=pbmc$time
# time=factor(time,levels = c("d0","d1","d2","d3","d7","d14","d21","d28"))
# pbmc$time=time
# DimPlot(pbmc,split.by = "time",ncol = 4)
