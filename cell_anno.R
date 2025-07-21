library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(ggsci)
library(patchwork) 

rm(list = ls()) 
gc()

data <- read.csv("/Users/mac/Desktop/08.1.scm6A_data/counts_Ctrl_MtoH_Scm6A_pred_results.csv", header = T,row.names=4)
data[1:3,1:8]
data <- data[,-c(1:6)]
data[1:3,1:3]
data <- as(data, "sparseMatrix")

merged <- CreateSeuratObject(counts = data, project = "sample")

merged <- NormalizeData(merged, 
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)

merged <- FindVariableFeatures(merged, 
                               selection.method = "vst",
                               nfeatures = 2500)

merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 40
                 , features = VariableFeatures(object = merged)
                 , verbose = F)

merged <- JackStraw(merged, num.replicate = 100)

merged <- ScoreJackStraw(merged, dims = 1:20) 

merged <- FindNeighbors(merged, dims = 1:20)

merged <- FindClusters(merged, resolution = 0.5)

head(Idents(merged), 5)


merged <- RunUMAP(merged, dims = 1:20) 
merged <- RunTSNE(merged, dims = 1:20)

cell_info_all <- read.csv("/Users/mac/Desktop/08.1.scm6A_data/cellanno_Ctrl.csv", header = TRUE,row.names = 1)

row.names(cell_info_all) <- gsub("-", ".", row.names(cell_info_all))

merged@meta.data <- cbind (merged@meta.data,cell_info_all[match(row.names(merged@meta.data),row.names(cell_info_all)),'cell_type_anno_new'])
B <- as.data.frame(merged@meta.data)    
B
colnames(merged@meta.data)[6]<- c("Ctrl m6A")
P <- DimPlot(merged, group.by= "Ctrl m6A", 
             raster=FALSE, pt.size=0.5,alpha = 0.8,reduction = "umap");P
ggsave("Ctrl_UMAP_celltype_anno_SingleR_m6A.pdf",P,width=8,height=6)


