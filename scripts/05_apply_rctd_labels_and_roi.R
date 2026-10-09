library(Banksy)
library(spacexr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(SeuratWrappers)
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
    Med = c("Fat", "SMC", "Mus"),
    Fib = "Fib",
    ISG_Endo = "ISG_Endo"
  )
  
  lookup <- setNames(rep(names(type_mapping), lengths(type_mapping)), 
                     unlist(type_mapping))
  lookup["unkown"] <- "unkown"  # Default value
  new=lookup[as.character(type)]
  return(new)
}
server <- function(input, output, session) {
  output$"plot-text" <- renderText("Plot: ")
  output$plot <- renderPlotly(p)
  output$selected <- renderPrint({
    d <- event_data("plotly_brushed")
    if (!is.null(d)) {
      saveRDS(d, "lassso_boundary.rds")
    }
  })
}
###################################################################################################
setwd("D://atlas/rctd//")
num="102"
sam=paste0("d",num)
load(paste0("D:/atlas/rctd/RCTD_cut",num,"_2.rdata"))
spa=readRDS(paste0("D:/atlas/b2c/",sam,"_spatial.rds"))

spa$sam=sam
dir.create(paste0("D:/atlas/b2c/",sam,"/"))
##上下翻转坐标
#DimPlot(spa,reduction = "spatial",raster = F)
spat=spa@reductions[["spatial"]]
emb=spat@cell.embeddings
emb[,2]=max(emb[,2])-emb[,2]+1
spat@cell.embeddings=emb
spa@reductions[["spatial"]]=spat
spa$x=emb[,1]
spa$y=emb[,2]
type=RCTD.full@results[["results_df"]]
spa$type=type$first_type
spa$sing=type$spot_class
table(spa$type)
table(spa$sing)
spa3=subset(spa,subset = sing=="reject",invert=T)
###
type2=merge_type(spa3$type)
names(type2)=names(spa3$type)
spa3$type2=type2
DimPlot(spa3,reduction = "spatial",group.by = "type2",raster = F)



##不规则区域#########################################################################
# 手动绘制多边形顶点（左键点击，右键结束）
library(shiny)
library(plotly)
map=data.frame(spa3$x,spa3$y)
colnames(map)=c("x","y")
ind=sample(rownames(map),10000)
tem=map[ind,]
###
ui <- fluidPage(
  plotlyOutput("plot"),
  verbatimTextOutput("selected")
)

p <-
  plot_ly(showlegend = T) %>%
  add_markers(data =tem, x=~x, y =~y, tooltip = c('x', 'y'))

shinyApp(ui, server)
lasso_boudnary <- readRDS("lassso_boundary.rds")
head(lasso_boudnary)
point_in_poly_idx <-
  sp::point.in.polygon(
    point.x =map$x,
    point.y = map$y,
    pol.x = lasso_boudnary$x,
    pol.y = lasso_boudnary$y
  )
points_selected <- map[point_in_poly_idx != 0, ]
spa3=subset(spa3,cells=rownames(points_selected))

DefaultAssay(spa3) <- "SCT"
FeaturePlot(spa3,features = "Cxcl10",reduction = "spatial",raster = F)
ggsave(filename = paste0("D:/atlas/b2c/",sam,"/Cxcl10.png"),width = 20,height = 20)
#FeaturePlot(spa3,features = "Cxcl9",reduction = "spatial",raster = F)
#ggsave(filename = paste0("D:/atlas/b2c/",sam,"/Cxcl9.png"),width = 20,height = 20)
DimPlot(spa3,reduction = "spatial",group.by = "type2",raster = F)
ggsave(filename = paste0("D:/atlas/b2c/",sam,"/type2.png"),width = 20,height = 20)
save(spa3,file = paste0("D:/atlas/",sam,".rdata"))
####################################################################################
library(sceasy)
library(reticulate)
#py_available()#[1] FALSE  #检查您的系统是否安装过Python
use_condaenv("D:\\miniconda3\\envs\\py39/")

py_config()#安装的python版本环境查看，显示anaconda和numpy的详细信息。
# version:        3.9.21 (main, Dec 11 2024, 16:35:24) [MSC v.1929 64 bit (AMD64)]
# Architecture:   64bit
# numpy:          D:/miniconda3/envs/py39/Lib/site-packages/numpy
# numpy_version:  1.26.4
py_available()#[1] TRUE   #检查您的系统是否安装过Python
py_module_available("pandas")#检查“pandas”是否安装

numpy <- import("numpy")
sceasy::convertFormat(spa3, from="seurat", to="anndata", assay = "Spatial", main_layer = "count",
                      outFile=paste0("D:/atlas/",sam,".h5ad"))     
table(spa3$type2)

