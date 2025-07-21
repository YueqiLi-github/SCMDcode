##加载包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

rm(list=ls())

UCEC <- Load10X_Spatial(data.dir = "/Users/mac/Documents/An/单细胞m6A数据库/plot/Fig3.空转/UCEC",
                        assay = "Spatial",
                        filename = "GSM6177623_NYU_UCEC3_Vis_processed_filtered_feature_bc_matrix.h5",
                        slice = "slice1",filter.matrix = TRUE,to.upper = FALSE
                        )

counts_matrix <- GetAssayData(UCEC, layer = "counts", assay = "Spatial")
counts_matrix_dense <- as.matrix(counts_matrix)
write.csv(counts_matrix_dense,'UCEC_count.csv')


library(Matrix)
csv_file_path = "/Users/mac/Documents/An/单细胞m6A数据库/plot/Fig3.空转/UCEC_matrix_Scm6A_pred_results.csv"
A <- read.csv(csv_file_path, header = TRUE, row.names = 4, check.names = FALSE)
m6A_df <- A[ , -c(1:6)]
rm(list=c('A','csv_file_path'))

sparse_matrix <- as(as.matrix(m6A_df), "dgCMatrix")
new_assay <- CreateAssayObject(counts = sparse_matrix)

UCEC[["m6A"]] <- new_assay
rm(list=c('sparse_matrix','new_assay','m6A_df'))

UCEC <- NormalizeData(UCEC, assay = "m6A", verbose = FALSE)
UCEC <- FindVariableFeatures(UCEC, assay = "m6A", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
UCEC <- ScaleData(UCEC, assay = "m6A", verbose = FALSE)

UCEC <- RunPCA(UCEC, assay = "m6A", verbose = FALSE)

ElbowPlot(UCEC,ndims = 50)

UCEC <- FindNeighbors(UCEC, assay = "m6A", dims = 1:50, verbose = FALSE, graph.name = "m6A_snn")
UCEC <- FindClusters(UCEC, graph.name = "m6A_snn", resolution = 0.1, verbose = FALSE)
UCEC <- RunUMAP(UCEC, assay = "m6A", dims = 1:50, verbose = FALSE)


p1 <- DimPlot(UCEC, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(UCEC, label = TRUE, label.size = 3)
p1+p2

SpatialFeaturePlot(UCEC, features = c(), 
              #     alpha = c(0.3, 1),
                   pt.size.factor = 2)

SpatialFeaturePlot(UCEC, features = c("METTL3", "METTL14","WTAP","FTO"), alpha = c(0.8, 1), pt.size.factor = 2)


csv_file_path = "/Users/mac/Documents/An/单细胞m6A数据库/plot/Fig3.空转/UCEC_matrix_RBP.csv"
RBP_df <- read.csv(csv_file_path, header = TRUE, row.names = 1, check.names = FALSE)

RBP_df[is.na(RBP_df)] <- 0
sparse_matrix <- as(as.matrix(RBP_df), "dgCMatrix")
new_assay <- CreateAssayObject(counts = sparse_matrix)
UCEC[["RBP"]] <- new_assay
UCEC <- NormalizeData(UCEC, assay = "RBP", verbose = FALSE)
UCEC <- FindVariableFeatures(UCEC, assay = "RBP", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
UCEC <- ScaleData(UCEC, assay = "RBP", verbose = FALSE)
UCEC <- RunPCA(UCEC, assay = "RBP", verbose = FALSE)
ElbowPlot(UCEC,ndims = 50)

UCEC <- FindNeighbors(UCEC, assay = "RBP", dims = 1:30, verbose = FALSE, graph.name = "RBP_snn")
UCEC <- FindClusters(UCEC, graph.name = "RBP_snn", resolution = 0.5, verbose = FALSE)
UCEC <- RunUMAP(UCEC, assay = "RBP", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(UCEC, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(UCEC, label = TRUE, label.size = 3)
p1+p2





UCEC
head(UCEC@meta.data)
str(UCEC)
head(UCEC@images[["slice1"]]@image)

plot1 <- VlnPlot(UCEC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(UCEC, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(UCEC, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(UCEC, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

UCEC[["percent.mt"]] <- PercentageFeatureSet(UCEC, pattern = "^mt[-]")
plot1 <- VlnPlot(UCEC, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(UCEC, features = "percent.mt") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

UCEC1 <- subset(UCEC, 
                 subset = nFeature_Spatial > 200 & 
                   nFeature_Spatial <7500 & 
                   nCount_Spatial > 1000 & 
                   nCount_Spatial < 60000 & 
                   percent.mt < 25)
UCEC1

plot1 <- VlnPlot(UCEC1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(UCEC1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

UCEC <- SCTransform(UCEC, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(UCEC, features = c("METTL3", "METTL14","WTAP","FTO"))

expression_matrix <- UCEC@assays[["SCT"]]@counts
expression_matrix_dense <- as.matrix(expression_matrix)
write.csv(expression_matrix_dense,'UCEC_matrix.csv')

p1 <- SpatialFeaturePlot(UCEC, features = c("CCL5","CD7","PTPRC","CXCR4","CD3D","GNLY","CST7","SYTL3","CD8A","NKG7","TRBC2","CYTIP","CD2","TSC22D3","KLRD1","SAMSN1","ARHGDIB","CLEC2D","SMAP2","CD96"), 
                         pt.size.factor = 2)
p2 <- SpatialFeaturePlot(UCEC, features = "METTL3", alpha = c(0.1, 1))
p1+p2

pdf("MRGGraph.pdf",width=4,height=4)
p1
dev.off()


UCEC <- RunPCA(UCEC, assay = "SCT", verbose = FALSE)
ElbowPlot(UCEC,ndims = 50)

UCEC <- FindNeighbors(UCEC, reduction = "pca", dims = 1:50)
UCEC <- FindClusters(UCEC, verbose = FALSE)
UCEC <- RunUMAP(UCEC, reduction = "pca", dims = 1:50)
p1 <- DimPlot(UCEC, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(UCEC, label = TRUE, label.size = 3)
p1+p2

SpatialDimPlot(UCEC, cells.highlight = CellsByIdentities(object = UCEC, idents = c(0, 1, 2, 3, 4, 5, 6)), facet.highlight = TRUE, ncol = 3)

ISpatialDimPlot(UCEC)
ISpatialFeaturePlot(UCEC,feature = "FTO")
LinkedDimPlot(UCEC)

de_markers <- FindMarkers(UCEC, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = UCEC, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
UCEC <- FindSpatiallyVariableFeatures(UCEC, assay = "SCT", features = VariableFeatures(UCEC)[1:1000], selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(UCEC, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(UCEC, features = top.features, ncol = 3, alpha = c(0.1, 1))

cortex <- subset(UCEC, idents = c(0, 1, 2, 3, 4, 5, 6))

cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1+p2
