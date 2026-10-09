library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(RColorBrewer)
DimPlot(pbmc,label = T)
FeaturePlot(pbmc,features = c("Siglecf","Cxcl9","Trem2","C1qa"))
FeaturePlot(pbmc,features = "Epcam")
S=pbmc$S
sum(S)
load("D:/space/HD1/pbmc_rctd.rdata")
for (i in 1:ncol(markers)) {
  gene=markers[,i]
  VlnPlot(pbmc2,features = gene,ncol = 2)
  ggsave(filename =paste0(colnames(markers)[i],".png"),width = 16,height = 12)
}
new.cluster.ids=c('Endo','T','Fib',"B","T","Endo","Neu","T",'IM_Mo','DC',"AM","NK","Endo","Endo",
                  "AT","SMC","Endo","Plasma","Mes","ciliate","AT","Mus")
new.cluster.ids=c('Endo','T','Fib',"B","T","Endo","Neu","T",'IM_Mo','DC',"AM","T","Endo","Endo",
                  "Epi","SMC","Endo","B","Epi","Epi","Epi","Mus")
new.cluster.ids
pbmc@active.ident=pbmc$clu1
pbmc$clu1=pbmc@active.ident
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc,label = T)
S=pbmc$S
s=ifelse(pbmc$S>0,"P","N")
table(s)
pbmc$S=s
DimPlot(pbmc,label = T,split.by = "S")
##ifelse()##打标签
tag=pbmc$tag
chip=pbmc$chip
table(chip)
chip[chip%in%c("A1","A2")]="chip1"
chip[chip%in%c("A5","A6")]="chip2"
pbmc$chip2=chip
sam=paste0(chip,"_",tag)
table(sam)
sam[sam%in%c("chip1_tag1","chip1_tag2","chip1_tag3","chip1_tag4")]="d0"
sam[sam%in%c("chip1_tag5","chip1_tag6","chip1_tag7","chip1_tag8")]="d2"
sam[sam%in%c("chip1_tag9","chip1_tag10","chip1_tag11","chip1_tag12")]="d45"
sam[sam%in%c("chip2_tag1","chip2_tag2","chip2_tag3","chip2_tag4")]="d10"
sam[sam%in%c("chip2_tag5","chip2_tag6","chip2_tag7","chip2_tag8")]="d21"
sam[sam%in%c("chip1_Mult","chip1_und","chip2_Mult","chip2_und")]="und"
table(sam)
pbmc$sam=sam
DimPlot(pbmc,label = T,group.by = "sam")
DimPlot(pbmc,label = T,split.by = "sam",ncol = 3)
pbmc$clum=pbmc@active.ident
pbmc$lev1=pbmc@active.ident
####T细胞再分选##############

nfkb=subset(pbmc,idents = c("T"))
nfkb=NormalizeData(nfkb)
nfkb=FindVariableFeatures(nfkb,selection.method = "vst",nfeatures = 2000)
nfkb=RunPCA(nfkb,features = VariableFeatures(nfkb))
ElbowPlot(nfkb,50)
nfkb=FindNeighbors(nfkb, dims = 1:20)
nfkb=RunUMAP(nfkb,dims = 1:20)
nfkb=FindClusters(nfkb,resolution =0.5)
DimPlot(nfkb,label = T)
DimPlot(nfkb,label = T,group.by = "clu2")
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
FeaturePlot(nfkb,features =markers$X8)
FeaturePlot(nfkb,features = c("Ly6c2"))
FeaturePlot(nfkb,features = "Cd4")
FeaturePlot(nfkb,features = "Cd8a")

FeaturePlot(nfkb,features = "Ccr7")
FeaturePlot(nfkb,features = "Itgae")
FeaturePlot(pbmc,features = "1700097N02Rik")
FeaturePlot(nfkb,features = "Rora")
FeaturePlot(nfkb,features = "Trgc2")
FeaturePlot(nfkb,features = "Ccr6")
FeaturePlot(nfkb,features = "Cd3d")
FeaturePlot(pbmc,features = "Klrb1c")
FeaturePlot(nfkb,features = "Tnfsf11")

FeaturePlot(nfkb,features = "Rora")
ma14=FindMarkers(nfkb,ident.1 = "14",ident.2 = "2",min.pct=0.25)
VlnPlot(nfkb, features = c("nFeature_RNA", "Contamination",
                              "percent.Mt",  "scDblFinder.score"), ncol = 4)
T1=nfkb
nfkb=T1

VlnPlot(nfkb, features = c("Ifng", "Irf7","Ifit3"), ncol = 4)
ma17=FindMarkers(nfkb,ident.1 = "1",ident.2 = "7",min.pct=0.25)
new=c("CTL1","T4","TRM8","Naive_T","NK","pro_T","pro_T","CTL2","KIR_T8","ILC3","ILC2","CTL1","RDT","Treg","TRM8","B","NK")
#CTL2具有耗竭倾向，ILC3为LTi
names(new) <- levels(nfkb)
nfkb <- RenameIdents(nfkb, new)
DimPlot(nfkb,label = T)|FeaturePlot(nfkb,features = "Cd8a")
##拆naiveT
mo=subset(nfkb,idents = c("Naive_T"))
mo=FindClusters(mo,resolution =0.4)
DimPlot(mo,label = T)
FeaturePlot(mo,features = c("Cd4","Cd8a"),blend = T)
new=c("Naive_T4","Naive_T8","Naive_T8","Naive_T4","Naive_T8","Naive_T4")
names(new) <- levels(mo)
mo <- RenameIdents(mo, new)
DimPlot(mo,label = T)
Idents(nfkb,cells=colnames(mo))=Idents(mo)
DimPlot(nfkb,label = T)

save(nfkb,file = "D:/atlas/xgy/T2.rdata")
Idents(pbmc,cells=colnames(nfkb))=Idents(nfkb)

DimPlot(pbmc,label = T)
########################################################################################
####=髓系细胞再分析##############
nfkb=subset(pbmc,idents = c("DC","IM_Mo","AM","Neu"))
nfkb=NormalizeData(nfkb)
nfkb=FindVariableFeatures(nfkb,selection.method = "vst",nfeatures = 2000)
nfkb=RunPCA(nfkb,features = VariableFeatures(nfkb))
ElbowPlot(nfkb,50)
nfkb=FindNeighbors(nfkb, dims = 1:20)
nfkb=RunUMAP(nfkb,dims = 1:20)
nfkb=FindClusters(nfkb,resolution =0.6)
DimPlot(nfkb,label = T)+NoLegend()|DimPlot(nfkb,label = T,group.by = "clu2")+NoLegend()
DimPlot(nfkb,label = T,group.by = "clu2")+NO_Legend()
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
FeaturePlot(nfkb,features =markers$X13)
FeaturePlot(pbmc,features =markers$X20)
FeaturePlot(nfkb,features = c("Ly6c2"))
FeaturePlot(nfkb,features = "Cxcl10")
FeaturePlot(nfkb,features = "Trem2")
FeaturePlot(nfkb,features = "Ccr7")
FeaturePlot(nfkb,features = "Itgae")
FeaturePlot(pbmc,features = "Cyp11a1")
FeaturePlot(nfkb,features = "Cd7")
FeaturePlot(nfkb,features = "Cd55")
FeaturePlot(nfkb,features = "Prg4")
FeaturePlot(nfkb,features = "Cd3d")
FeaturePlot(nfkb,features = "Mki67")
FeaturePlot(nfkb,features = "Ms4a1")
FeaturePlot(nfkb,features = "Cyp11a1")
VlnPlot(nfkb, features = c("nFeature_RNA", "Contamination",
                           "percent.Mt",  "scDblFinder.score"), ncol = 4)
VlnPlot(nfkb, features = c( "nFeature_RNA","Contamination"))
VlnPlot(nfkb, features = c("Ifng", "Irf7","Ifit3"), ncol = 4)
ma19=FindMarkers(nfkb,ident.1 = "19",min.pct=0.25)
sam2=paste0(nfkb$chip,nfkb$tag)
table(sam2)
new=c("Neu","AM","IM1","DC3","Mo","MoDC","DC2","DC1","Trem2","low","inf_AM","Ace_Mo","Neu","Neu","pro_DC","pro_AM","IM2","low","pDC","Mast","low","low")
#CTL2具有耗竭倾向，ILC3为LTi
names(new) <- levels(nfkb)
nfkb <- RenameIdents(nfkb, new)
DimPlot(nfkb,label = T)
mo=subset(nfkb,idents = c("Mo"))
mo=FindClusters(mo,resolution =0.1)
DimPlot(mo,label = T)
FeaturePlot(mo,features = c("Ccr2","Cxcl9","Ly6c2","Plac8","Irf7","Cxcl10"))
new=c("inf_Mo","Mo")
names(new) <- levels(mo)
mo <- RenameIdents(mo, new)
DimPlot(mo,label = T)
Idents(nfkb,cells=colnames(mo))=Idents(mo)
DimPlot(nfkb,label = T)
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
save(nfkb,file = "D:/atlas/xgy/M2.rdata")
Idents(pbmc,cells=colnames(nfkb))=Idents(nfkb)
DimPlot(pbmc,label = T)
########################################################################################
####B细胞再分析##############
nfkb=subset(pbmc,idents = c("B"))
#########去ig
gene=rownames(nfkb)
gene=as.data.frame(gene)
gene2=grepl(pattern = "^Tr[ab][vdj]",gene$gene)
gene$gene2=gene2
geneB=grepl(pattern = "^Ig[hlk][vdj]",gene$gene)
gene$geneB=geneB
gene$ind=gene$gene2|gene$geneB
genei=gene$gene[gene$ind]
gene2=gene$gene[!gene$ind]
#a2=a
#a=a[,gene2]
dat3=nfkb@assays$RNA@layers$counts
dimnames(dat3)=dimnames(nfkb)
dim(dat3)
dati=dat3[gene$ind,]
dat2=dat3[!gene$ind,]
dim(dat2)
table(gene$ind)
a3=CreateSeuratObject(dat2)
head(rownames(a3))
head(colnames(a3))
a3[["imm"]] <-CreateAssayObject(counts =dati)
nfkb=a3
nfkb1=subset(pbmc,idents = c("B"))
nfkb@meta.data=cbind(nfkb@meta.data,nfkb1@meta.data)
##############################################
nfkb=NormalizeData(nfkb)
nfkb=FindVariableFeatures(nfkb,selection.method = "vst",nfeatures = 2000)
nfkb=ScaleData(nfkb)
nfkb=RunPCA(nfkb,features = VariableFeatures(nfkb))
ElbowPlot(nfkb,50)
nfkb=FindNeighbors(nfkb, dims = 1:15)
nfkb=RunUMAP(nfkb,dims = 1:15)
nfkb=FindClusters(nfkb,resolution =0.1)
DimPlot(nfkb,label = T)
DimPlot(nfkb,label = T,group.by = "clu2")
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
FeaturePlot(nfkb,features =markers$X13)
FeaturePlot(pbmc,features =markers$X6)
FeaturePlot(nfkb,features = c("Tifa"))
FeaturePlot(nfkb,features = "Cd79a")
FeaturePlot(nfkb,features = "Trem2")
FeaturePlot(nfkb,features = "Ccr7")
FeaturePlot(nfkb,features = "Cxcr5")
FeaturePlot(pbmc,features = "Vpreb3")
FeaturePlot(nfkb,features = "Vpreb3")
FeaturePlot(nfkb,features = "Aicda")
FeaturePlot(nfkb,features = "Bcl6")
FeaturePlot(nfkb,features = "Fcer2a")
FeaturePlot(nfkb,features = "Mki67")
FeaturePlot(nfkb,features = "Rag1")
FeaturePlot(nfkb,features = "Acan")
VlnPlot(nfkb, features = c("nFeature_RNA", "Contamination",
                           "percent.Mt",  "scDblFinder.score"), ncol = 4)
DimPlot(nfkb,label = T)+NoLegend()|DimPlot(nfkb,label = T,group.by = "clu2")+NoLegend()
VlnPlot(nfkb, features = c("Ifng", "Irf7","Ifit3"), ncol = 4)
ma17=FindMarkers(nfkb,ident.1 = "1",ident.2 = "0",min.pct=0.25)
ma17=FindMarkers(nfkb,ident.1 = "5",min.pct=0.25)
ma17=ma17[order(ma17$avg_log2FC),]
head(rownames(ma17),10)
tail(rownames(ma17),10)
sam2=paste0(nfkb$chip,nfkb$tag)
table(sam2)
new=c("naive_B","mem_B","Plasma","GCB","double","pre_B","double")
names(new) <- levels(nfkb)
nfkb <- RenameIdents(nfkb, new)
DimPlot(nfkb,label = T)
mo=subset(nfkb,idents = c("GCB"))
mo=FindClusters(mo,resolution =0.2)
DimPlot(mo,label = T)
DimPlot(mo,label = T,group.by = "clu2")
FeaturePlot(mo,features = c("Aicda","Fas","Cxcr5","Cxcr4","Bcl6","Cxcl10"))
new=c("pro_B","pro_B","GCB","GCB")
names(new) <- levels(mo)
mo <- RenameIdents(mo, new)
DimPlot(mo,label = T)
Idents(nfkb,cells=colnames(mo))=Idents(mo)
DimPlot(nfkb,label = T)
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
save(nfkb,file = "D:/atlas/xgy/B2.rdata")
##排除之前认为的双细胞再重新命名
cell1= colnames(nfkb)
cell=colnames(pbmc)
ind=cell1%in%cell
cell2=cell1[ind]
nfkb2=subset(nfkb,cells=cell2)
DimPlot(nfkb2,label = T)
Idents(pbmc,cells=colnames(nfkb2))=Idents(nfkb2)
DimPlot(pbmc,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(53))
length(table(Idents(pbmc)))
#########################################################################################
####内皮细胞再分选##############
nfkb=subset(pbmc,idents = c("Endo"))
nfkb=NormalizeData(nfkb)
nfkb=FindVariableFeatures(nfkb,selection.method = "vst",nfeatures = 2000)
nfkb=RunPCA(nfkb,features = VariableFeatures(nfkb))
ElbowPlot(nfkb,50)
nfkb=FindNeighbors(nfkb, dims = 1:20)
nfkb=RunUMAP(nfkb,dims = 1:20)
nfkb=FindClusters(nfkb,resolution =0.2)
DimPlot(nfkb,label = T)
DimPlot(nfkb,label = T,group.by = "clu2")
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
DimPlot(nfkb,label = T)+NoLegend()|DimPlot(nfkb,label = T,group.by = "clu2")+NoLegend()
FeaturePlot(nfkb,features =markers$X9)
FeaturePlot(pbmc,features =markers$X9)
FeaturePlot(nfkb,features = c("Ly6c2"))
FeaturePlot(nfkb,features = "Cd4")
FeaturePlot(nfkb,features = "Prox1")
FeaturePlot(nfkb,features = "Col1a1")
FeaturePlot(nfkb,features = "Spns2")
FeaturePlot(pbmc,features = "Spns2")
FeaturePlot(nfkb,features = "Gpr182")
FeaturePlot(nfkb,features = "Ifitm3")
FeaturePlot(nfkb,features = "Tmem252")
FeaturePlot(nfkb,features = "Clca3a1")
FeaturePlot(pbmc,features = "Clca3a1")
FeaturePlot(nfkb,features = "Pecam1")
FeaturePlot(nfkb,features = "Cd74")
VlnPlot(nfkb, features = c("nFeature_RNA", "Contamination",
                           "percent.Mt",  "scDblFinder.score"), ncol = 4)

VlnPlot(nfkb, features = c("Ifng", "Irf7","Ifit3"), ncol = 4)
ma8=FindMarkers(nfkb,ident.1 = "8",min.pct=0.25)
##endo6
new=c("Endo1","Endo2","Endo1","ISG_Endo","Endo1","Endo3","Pericyte","LEC","Endo1","Fat","Endo1")
#endo1为经典内皮、endo2为新生成内皮，endo3为气管内皮
names(new) <- levels(nfkb)
nfkb <- RenameIdents(nfkb, new)
DimPlot(nfkb,label = T)
save(nfkb,file = "D:/atlas/xgy/Endo2.rdata")
Idents(pbmc,cells=colnames(nfkb))=Idents(nfkb)
DimPlot(pbmc,label = T)
####上皮细胞再分选#################################################################################
nfkb=subset(pbmc,idents = c("Epi"))
nfkb=NormalizeData(nfkb)
nfkb=FindVariableFeatures(nfkb,selection.method = "vst",nfeatures = 2000)
nfkb=RunPCA(nfkb,features = VariableFeatures(nfkb))
ElbowPlot(nfkb,50)
nfkb=FindNeighbors(nfkb, dims = 1:15)
nfkb=RunUMAP(nfkb,dims = 1:15)
nfkb=FindClusters(nfkb,resolution =0.2)
DimPlot(nfkb,label = T)+NoLegend()|DimPlot(nfkb,label = T,group.by = "clu2")+NoLegend()
DimPlot(nfkb,label = T)
DimPlot(nfkb,label = T,group.by = "clum")
DimPlot(nfkb,label = T,ncol = 3,split.by = "sam")
FeaturePlot(nfkb,features =markers$X9)
FeaturePlot(pbmc,features =markers$X6)
FeaturePlot(nfkb,features = c("Ly6c2"))
FeaturePlot(nfkb,features = "Cd4")
FeaturePlot(nfkb,features = "Cd8a")
FeaturePlot(nfkb,features = "Cxcl10")
FeaturePlot(nfkb,features = "Cpe")
FeaturePlot(pbmc,features = "Fabp4")
FeaturePlot(nfkb,features = "Gpr182")
FeaturePlot(nfkb,features = "Spns2")
FeaturePlot(nfkb,features = "Car8")
FeaturePlot(nfkb,features = "Lbp")
FeaturePlot(pbmc,features = "Spns2")
FeaturePlot(nfkb,features = "Pecam1")
FeaturePlot(nfkb,features = "Rora")
VlnPlot(nfkb, features = c("nFeature_RNA", "Contamination",
                           "percent.Mt",  "scDblFinder.score"), ncol = 4)

VlnPlot(nfkb, features = c("Ifng", "Irf7","Ifit3"), ncol = 4)
ma17=FindMarkers(nfkb,ident.1 = "1",ident.2 = "7",min.pct=0.25)
##endo6
new=c("AT2","Mes","Ciliate","AT1","Club","low","PNEC")
names(new) <- levels(nfkb)
nfkb <- RenameIdents(nfkb, new)
DimPlot(nfkb,label = T)
save(nfkb,file = "D:/atlas/xgy/epi2.rdata")
Idents(pbmc,cells=colnames(nfkb))=Idents(nfkb)
DimPlot(pbmc,label = T)
DimPlot(pbmc2,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(49))
club=FindMarkers(pbmc2,ident.1 = "Club",only.pos = T,logfc.threshold = 0.25)
FeaturePlot(nfkb,features = head(rownames(club)))
##########################
FLS=FindMarkers(pbmc,ident.1 = "FLS",only.pos = T,logfc.threshold = 0.25)
FLS=FLS[FLS$p_val==0,]
rownames(FLS)
VlnPlot(pbmc,features = c("Prg4"),pt.size=0)
pbmx3=pbmc
DimPlot(pbmc)
FeaturePlot(pbmc,features = c("Prg4","C1qa"))
pbmx3 <- RenameIdents(object = pbmx3,  'FLS' = 'Prg4_M')
pbmc=pbmx3
pbmc3=subset(pbmc,idents =c("low","double"), invert = TRUE)
DimPlot(pbmc3,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(51))
pbmc=pbmc3
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:32)
#pbmc  <- FindClusters(pbmc , resolution = 0.2)
pbmc  <- RunUMAP(pbmc , reduction = "pca", dims = 1:32)
DimPlot(pbmc,label = T,cols=colorRampPalette(brewer.pal(9,'Set1'))(52),raster=FALSE)
ggsave(filename = "pbmc_UMAP.pdf",width = 16,height = 12,dpi = 300,units = "in")
FeaturePlot(pbmc,features = c("Tgfb1","Tgfb2","Tgfb3"),raster=FALSE,cols = c("grey", "red"))
save(pbmc,file = "D:/atlas/xgy/pbmc.rdata")            
`table(pbmc@active.ident)
DimPlot(pbmc,split.by = "sam",ncol = 3,label = T)+NoLegend()
DimPlot(pbmc)
table(spa1$type)
######################################################################################
sample.markersm  <- FindAllMarkers(object = pbmc, only.pos = TRUE,
                                   min.pct = 0.05,
                                   thresh.use = 0.5)
save(sample.markersm,file = "xgy/marker.rdata")
top3 <- sample.markersm%>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 
write.csv(top3,file = "xgy/top5.csv")
lev=levels(pbmc$clu)
lev
top1=  mutate(top1,cluster = factor(cluster, levels = lev)) %>% # 将 cluster 转换为 factor 并指定 levels
  arrange(cluster)
markers_51 <- c(
  "Ighd", "Cd80", "Mki67","Aicda", "Jchain", "Rag1", "Lef1", "Ccr7",  "Gzmk", 
  "Lag3", "Itgae", "Klrc1", "Cd4", "Foxp3", "Rorc", "Gata3", "Trdc", "Nkg7", "Top2a", 
  "Xcr1", "Cd1c", "Fscn1", "Siglech", "Cd14", "F13a1", "Cxcl10", "Treml4", "C1qb", "Folr2", 
  "Marco", "Anln", "Cfb", "Trem2", "S100a8", "Ms4a2", "Rtkn2", "Sftpc", "Caps", "Scgb1a1", 
  "Sema3g", "Ednrb", "Slc6a2", "Ifit3", "Ccl21a", "Fabp4", "Mb", "Acta2", "Notch3", "Upk3b", "Col1a1"
)


# 确保 Ident 顺序正确
Idents(pbmc) <- factor(Idents(pbmc), levels = lev)

# 绘图
DotPlot(pbmc, features = markers_51 , cols = c("lightgrey", "#B33900"), dot.scale = 8) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
spam$type=factor(spam$type,levels =lev)
ggsave(filename = "fig/pbmc_marker.pdf")
table(spam$type)
DotPlot(spam, features = markers_51 , cols = c("lightgrey", "#B33900"), dot.scale = 8,group.by = "type") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
ggsave(filename = "fig/hd_marker.pdf")
DimPlot(pbmc,raster=T,
       group.by = "clu",raster.dpi = c(600, 600),
        cols =min_col,
        label = TRUE)
ggsave(filename = "fig/pbmc_umap.pdf") 
# 群太多时建议关掉legend，直接看label
# 先生成基础的矢量图对象（不要加 raster=TRUE）
p <- DimPlot(pbmc, raster = FALSE, group.by = "clu", cols = min_col, label = TRUE)

# 1. 第一步：导出高清 TIFF (作为底图，dpi 设为 600 保证点足够细腻)
ggsave(filename = "fig/pbmc_umap.tif", plot = p, device = "tiff", dpi = 600, width = 7, height = 6)

# 2. 第二步：黑科技！删掉图像中的“点”图层，只保留坐标轴、背景和图例
p_vector <- p
p_vector$layers[[1]] <- NULL  # 强制删除包含细胞点的数据图层

# 3. 导出空壳的纯矢量 PDF (由于没有了几万个点，文件只有几 KB，完全不卡)
ggsave(filename = "fig/pbmc_umap_frame.pdf", plot = p, width = 8, height = 6)
