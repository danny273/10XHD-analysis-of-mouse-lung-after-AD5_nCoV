##表达模式变化
spa2=list()
for (i in 1:length(spa)){
  a=spa[[i]]
  a=NormalizeData(a)
  a$source=names(spa)[i]
  spa2[[i]]=a
}
spa=spa2
rm(spa2)
##合并数据)
names(spa)=c("d001","d011","d101","d102","d211","d212","d451","d452")
for (i in 1:length(spa)){
  a=spa[[i]]
  a$source=names(spa)[i]
  spa[[i]]=a
}
spam=merge(spa[[1]],y=spa[2:length(spa)],add.cell.ids = names(spa),project = "spa")
gc()
ind=leiden$type==spam$type
cell_id=colnames(spam)
table(cell_id)
view=data.frame(leiden,spam$sam,spam$type2,ind)
vi2=view[view$ind==F,]
table(vi2$source)
table(vi2$spam.sam)
table(ind)
table(spam$nich2)

leiden2=leiden
table(leiden2$source)
unique(leiden2$source)
leiden2$source=factor(leiden2$source,levels = c("d001","d011","d101","d102","d211","d212","d451","d452"))
leiden2=leiden2[order(leiden2$source),]
ind=leiden2$type==spam$type
table(ind)
rownames(leiden2)=colnames(spam)
spam$type2=factor(spam$type2,levels = c("B","T","DC","Mo","Mac","Gran","AT","Airway","Endo","Stroma","Fib"))
save(spam,file = "D:/atlas/spam.Rdata")
####
leiden2$x=spam$x
leiden2$y=spam$y
write.csv(leiden2,file = "D:/atlas/leiden.csv",row.names = F)
###############################################################################################
Tc=subset(spam,subset=nich2%in%c("CTL","BALT","inf_AM","inf_endo"))
table(spam$nich2)
TcT=subset(Tc,subset=type2%in%c("T"))
TcT=NormalizeData(TcT)
VlnPlot(TcT,features = "Ccr7",group.by = "nich2")
ma=FindMarkers(TcT,ident.1 = "BALT",logfc.threshold = 0.25,group.by = "nich2",min.pct=0.05)
ma=FindAllMarkers(TcT,logfc.threshold = 0.25,group.by = "nich2",min.pct=0.05)
colnames(pbmc@meta.data)
table(pbmc$lev1)
table(pbmc$clu2)
pbmc@active.ident=pbmc$lev1
marker=FindAllMarkers(pbmc,logfc.threshold = 0.25,min.pct=0.05,only.pos = T)
Tmarker=marker[marker$cluster=="T",]
tgene=Tmarker$gene
ma=ma[rownames(ma)%in%tgene,]
ma=ma[ma$p_val<0.05,]
DoHeatmap(spam,features = c("Cxcl10","Cxcl9","Ccl5"), group.by = "nich3", label = TRUE, size =0, raster = FALSE,slot = "data")

ggsave()