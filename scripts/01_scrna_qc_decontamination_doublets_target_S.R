# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# install.packages("DoubletFinder")
# update.packages("Seurat")
# BiocManager::install("Seurat")
# BiocManager::install('glmGamPoi')
# library(DoubletFinder)
library(dplyr)
library(patchwork)
library(tidyverse)
library(Seurat)
library(patchwork)
#install.packages("cowplot")
library("cowplot")
#install.packages("stringr")
library("stringr")
library("parallel") 
library(doParallel)
library(foreach)
library(decontX)
library("scDblFinder")
rm(list=ls())
setwd("D:/atlas/xgy/")
file=list.files(path = "RNA",pattern = "RNA")
#ind=grepl(pattern = "tag[1/2/3]",x=file)
#file=file[!ind]
sample=paste0("RNA/",file)
a=detectCores(logical = F)
cl <- makeCluster(a)
registerDoParallel(cl)
dat=foreach(l = 1:length(file),.packages = c("Seurat")) %dopar% {
 Read10X(sample[l])
  
}
names(dat)=file
# for (i in 1:length(dat)) {
#   if(is.list(dat[[i]])){
#     dat[[i]]=dat[[i]]$`Gene Expression`
#   }
# }
stopImplicitCluster()
stopCluster(cl)
#################################################################
######################
dat2=list()
#######################################################
########################################################
fun <- function(i,dat=dat,file=file){
  a=dat[[i]]
  cbmc <- CreateSeuratObject(counts = a)
  Assays(cbmc)
  # #############################去ig
  # gene=rownames(cbmc)
  # gene=as.data.frame(gene)
  # gene2=grepl(pattern = "^Tr[ab][vdj]",gene$gene)
  # gene$gene2=gene2
  # geneB=grepl(pattern = "^Ig[hlk][vdj]",gene$gene)
  # gene$geneB=geneB
  # gene$ind=gene$gene2|gene$geneB
  # genei=gene$gene[gene$ind]
  # gene2=gene$gene[!gene$ind]
  # #a2=a
  # #a=a[,gene2]
  # dat3=cbmc@assays$RNA@layers$counts
  # dimnames(dat3)=dimnames(cbmc)
  # dim(dat3)
  # dati=dat3[gene$ind,]
  # dat2=dat3[!gene$ind,]
  # dim(dat2)
  # table(gene$ind)
  # a3=CreateSeuratObject(dat2)
  # head(rownames(a3))
  # head(colnames(a3))
  # a3[["imm"]] <-CreateAssayObject(counts =dati)  
  # cbmc=a3
  # #############标准流程
  # ##蛋白
  # adt_assay <- CreateAssay5Object(counts = a$`Antibody Capture`)
  # cbmc[["ADT"]] <- adt_assay
  # DefaultAssay(cbmc) <- "ADT"
  # cbmc=NormalizeData(cbmc, normalization.method = 'CLR', margin = 2) %>%
  #   ScaleData()
  #########################################
  DefaultAssay(cbmc) <- "RNA"
  cbmc <- NormalizeData(cbmc)
  cbmc <- FindVariableFeatures(cbmc)
  cbmc <- ScaleData(cbmc)
  cbmc <- RunPCA(cbmc, verbose = FALSE)
  cbmc <- FindNeighbors(cbmc, dims = 1:30)
  cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
  cbmc <- RunUMAP(cbmc, dims = 1:30)
  p3=DimPlot(cbmc, label = TRUE)
  ##去RNA污染和双细胞
  decontX_results <- decontX(GetAssayData(cbmc)) 
  cbmc$Contamination =decontX_results$contamination
  
  sce <- scDblFinder(GetAssayData(cbmc, slot="counts"),
                     clusters=Idents(cbmc),
                     dbr = 0.05)
  cbmc$scDblFinder.score <- sce$scDblFinder.score
  p1=FeaturePlot(cbmc,features = "scDblFinder.score")
  p2=FeaturePlot(cbmc,features = "Contamination")
  #ma5=FindMarkers(cbmc,ident.1 = "5",assay = "RNA")
  #FeaturePlot(cbmc,features = "Cd4")
  # Normalize ADT data,
  #FeatureScatter(cbmc,feature1 = "scDblFinder.score",feature2 ="Contamination" )
  #ma6=FindMarkers(cbmc,ident.1 = 4)
  p4=FeaturePlot(cbmc,features = "Top2a")
  
  wrap_plots(list(p3,p4,p1,p2),ncol = 2)
  ggsave(paste0(file[i],".png"))
  #dat2[[i]]=cbmc
  print(file[i])
  return(cbmc)
}
a=detectCores(logical = F)
cl <- makeCluster(20)
registerDoParallel(cl)
dat2=foreach(i = 1:length(file),.packages = c("Seurat","dplyr","decontX","scDblFinder","ggplot2","patchwork")) %dopar% {
  fun(i,dat = dat,file = file)
}
stopImplicitCluster()
stopCluster(cl)
names(dat2)=names(dat)

#########################################################################################
#读取tag信息##########################
ltag=list()
file2=list.files(path = "TAG")

sam2=paste0("TAG/",file2)
for (i in 1:length(file2)) {
  a=read.csv(sam2[i],header = T,sep = "\t")
  print(table(a$tag))
  ltag[[i]]=a
  names(ltag)[i]=file2[i]
}
lS=list()
file2=list.files(path = "Target")

sam2=paste0("Target/",file2)
for (i in 1:length(file2)) {
  a=Read10X(sam2[i])
  a=as.data.frame(a)
  a=as.data.frame(t(a))
  b=a$AD5_COV_R2
  names(b)=rownames(a)
  lS[[i]]=b
}

##合并所有信息
dat3=list()
for (i in 1:length(file)) {
  a=dat2[[i]]
  b=ltag[[i]]
  b2=b
  rownames(b)=b$X
  b$X=NULL
  sam=b$tag
  b2=b
  b2$tag=NULL
  table(sam)
  tag=sam
  tag[tag=="Multiplet"]="Mult"
  tag[tag=="Undetermined"]="und"
  tag[tag=="Sample_Tag1_Mouse"]="tag1"
  tag[tag=="Sample_Tag2_Mouse"]="tag2"
  tag[tag=="Sample_Tag3_Mouse"]="tag3"
  tag[tag=="Sample_Tag4_Mouse"]="tag4"
  tag[tag=="Sample_Tag5_Mouse"]="tag5"
  tag[tag=="Sample_Tag6_Mouse"]="tag6"
  tag[tag=="Sample_Tag7_Mouse"]="tag7"
  tag[tag=="Sample_Tag8_Mouse"]="tag8"
  tag[tag=="Sample_Tag9_Mouse"]="tag9"
  tag[tag=="Sample_Tag10_Mouse"]="tag10"
  tag[tag=="Sample_Tag11_Mouse"]="tag11"
  tag[tag=="Sample_Tag12_Mouse"]="tag12"
 table(tag)
  
  a$tag=tag
 # a@assays$TAG=CreateAssayObject(counts = b2)
  a$S=lS[[i]]
  dat3[[i]]=a
  a$chip=rep(file[i],ncol(a))
  }
names(dat3)=file
save(dat3,file = "raw_pbmc.rdata")
