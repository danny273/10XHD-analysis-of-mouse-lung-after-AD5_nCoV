# Load required libraries
library(spacexr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(progeny)
#library(Rfast)
setwd("~/space/HD1")
library("arrow")
# Data Preprocessing
# 解析命令行参数
# Load the Spatial data
spa1 <- readRDS("~/space/HD1/d21_1_b2c4rds_spatial.rds")
#pdac <- readRDS("/path/to/pdac_most_updated.rds")
spa1=d212_spatial 
# Load the Single cell data
#load("~/space/HD1/pbmc.rdata")
load("~/space/pbmc.rdata")
DimPlot(pbmc,label = T)
sc=subset(pbmc,idents = "PNEC", invert = T)
sc[["RNA"]] <- as(object = sc[["RNA"]], Class = "Assay")
sc$clu=Idents(sc)
# Prepare data for RCTD
counts <- sc@assays$RNA@counts
#Idents(sc) <- "SCT"
cluster <- as.factor(sc$clu)
names(cluster) <- colnames(sc)
nUMI <- sc$nCount_RNA
names(nUMI) <- colnames(sc)

# Create the reference object
reference <- Reference(counts, cluster, nUMI)
spa1=pdac
# Prepare spatial transcriptomics data
x=spa1@assays[["Spatial"]]@counts@x
min(x)
k=1/min(x)
x=x*k
x=round(x)
dim(spa1@assays[["Spatial"]]@counts)
#x=as.integer(x)
spa1@assays[["Spatial"]]@counts@x=x
DefaultAssay(spa1) <- "Spatial"
counts <- spa1@assays$Spatial@counts
dimnames(counts) <- dimnames(spa1)
# pos <- GetTissueCoordinates(spa)
# # coordinates_list <- lapply(spa, function(image_name) {
# #   pos <- GetTissueCoordinates(spa, image = image_name)
# #   colnames(pos) <- c('x','y')
# #   return(pos)
# # })
# coords <- pos
# rownames(coords) <- gsub("^.+\\.", "", rownames(coords))
# coords$cell=NULL
coords1 <-data.frame(spa1$array_col,spa1$array_row)
colnames(coords1) <- c('x','y')
#rownames(coords1) <- paste0("cell",rownames(coords1))
nUMI <- colSums(counts) 
head(nUMI)
query <- SpatialRNA(coords1, counts, colSums(counts))
RCTD <- create.RCTD(query, reference, max_cores = 24)
RCTD.full <- run.RCTD(RCTD, doublet_mode = "doublet")
save(RCTD.full,file = "RCTD_cut212_2.rdata")
