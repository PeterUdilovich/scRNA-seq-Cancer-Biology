install.packages('Seurat')
install.packages("dplyr")
install.packages("patchwork")
install.packages("SeuratObject")

#Additional recommended packages (speed and performance)
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

#Load the libraries
library(dplyr)
library(Seurat)
library(patchwork)


# Load the datasets ####
benign_breast.data <- Read10X_h5("/Users/test/Desktop/Compbio/main_rna_project/project_data/benign_breast/filtered_feature_bc_matrix.h5")
tumor_breast.data <- Read10X(data.dir = "/Users/test/Desktop/Compbio/main_rna_project/project_data/malignant_breast_sample")
benign_icc.data <- read.csv(
  "//Users/test/Desktop/Compbio/main_rna_project/project_data/adjacent_tissue_sample.csv",
  row.names = 1,        # first column = gene names
  check.names = FALSE   # keep barcodes exactly as they are
)
# 2. turn into matrix
benign_icc.data <- as.matrix(benign_icc.data)


tumor_icc.data <- read.csv(
  "//Users/test/Desktop/Compbio/main_rna_project/project_data/icc_tumor_tissue_sample.csv",
  row.names = 1,        # first column = gene names
  check.names = FALSE   # keep barcodes exactly as they are
)
tumor_icc.data <- as.matrix(tumor_icc.data)


#Seurat Objects ####
benign_breast <- CreateSeuratObject(counts = benign_breast.data, project = "bb", min.cells = 3, min.features = 200)
tumor_breast <- CreateSeuratObject(counts = tumor_breast.data, project = "tb", min.cells = 3, min.features = 200)
benign_icc <- CreateSeuratObject(counts = benign_icc.data, project = "bicc", min.cells = 3, min.features = 200)
tumor_icc <- CreateSeuratObject(counts = tumor_icc.data, project = "ticc", min.cells = 3, min.features = 200)

#Initial QC ####
benign_breast[["percent.mt"]] <- PercentageFeatureSet(benign_breast, pattern = "^MT-")
tumor_breast[["percent.mt"]] <- PercentageFeatureSet(tumor_breast, pattern = "^MT-")
benign_icc[["percent.mt"]] <- PercentageFeatureSet(benign_icc, pattern = "^MT-")
tumor_icc[["percent.mt"]] <- PercentageFeatureSet(tumor_icc, pattern = "^MT-")

#Violin Plot
VlnPlot(benign_breast, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(tumor_breast, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(benign_icc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(tumor_icc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Scatter Plots
plot1 <- FeatureScatter(benign_breast, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(benign_breast, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(tumor_breast, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tumor_breast, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(benign_icc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(benign_icc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(tumor_icc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tumor_icc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filter based on QC
benign_breast <- subset(benign_breast, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
tumor_breast <- subset(tumor_breast, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
benign_icc <- subset(benign_icc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
tumor_icc <- subset(tumor_icc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize the Data ####

benign_breast <- NormalizeData(benign_breast, normalization.method = "LogNormalize", scale.factor = 10000)
tumor_breast <- NormalizeData(tumor_breast, normalization.method = "LogNormalize", scale.factor = 10000)
benign_icc <- NormalizeData(benign_icc, normalization.method = "LogNormalize", scale.factor = 10000)
tumor_icc <- NormalizeData(tumor_icc, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable features ####
benign_breast <- FindVariableFeatures(benign_breast, selection.method = "vst", nfeatures = 2000)
tumor_breast <- FindVariableFeatures(tumor_breast, selection.method = "vst", nfeatures = 2000)
benign_icc <- FindVariableFeatures(benign_icc, selection.method = "vst", nfeatures = 2000)
tumor_icc <- FindVariableFeatures(tumor_icc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
bbtop10 <- head(VariableFeatures(benign_breast), 10)
tbtop10 <- head(VariableFeatures(tumor_breast), 10)
bicctop10 <- head(VariableFeatures(benign_icc), 10)
ticctop10 <- head(VariableFeatures(tumor_icc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(benign_breast)
plot2 <- LabelPoints(plot = plot1, points = bbtop10, repel = TRUE)
plot1 + plot2
plot1 <- VariableFeaturePlot(tumor_breast)
plot2 <- LabelPoints(plot = plot1, points = tbtop10, repel = TRUE)
plot1 + plot2
plot1 <- VariableFeaturePlot(benign_icc)
plot2 <- LabelPoints(plot = plot1, points = bicctop10, repel = TRUE)
plot1 + plot2
plot1 <- VariableFeaturePlot(tumor_icc)
plot2 <- LabelPoints(plot = plot1, points = ticctop10, repel = TRUE)
plot1 + plot2

# Scaling the data ####
bb_all.genes <- rownames(benign_breast)
benign_breast <- ScaleData(benign_breast, features = bb_all.genes)
tb_all.genes <- rownames(tumor_breast)
tumor_breast <- ScaleData(tumor_breast, features = tb_all.genes)
bicc_all.genes <- rownames(benign_icc)
benign_icc <- ScaleData(benign_icc, features = bicc_all.genes)
ticc_all.genes <- rownames(tumor_icc)
tumor_icc <- ScaleData(tumor_icc, features = ticc_all.genes)

head(tumor_icc[["RNA"]]$scale.data)

#Perform linear dimensional reduction ####

benign_breast <- RunPCA(benign_breast, features = VariableFeatures(object = benign_breast))
tumor_breast <- RunPCA(tumor_breast, features = VariableFeatures(object = tumor_breast))
benign_icc <- RunPCA(benign_icc, features = VariableFeatures(object = benign_icc))
tumor_icc <- RunPCA(tumor_icc, features = VariableFeatures(object = tumor_icc))

DimPlot(benign_breast, reduction = "pca") + NoLegend()
DimPlot(tumor_breast, reduction = "pca") + NoLegend()
DimPlot(benign_icc, reduction = "pca") + NoLegend()
DimPlot(tumor_icc, reduction = "pca") + NoLegend()

#Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(benign_breast)
ElbowPlot(tumor_breast)
ElbowPlot(benign_icc)
ElbowPlot(tumor_icc)

#Cluster the cells ####

benign_breast <- FindNeighbors(benign_breast, dims = 1:12)
benign_breast <- FindClusters(benign_breast, resolution = 0.5)

tumor_breast <- FindNeighbors(tumor_breast, dims = 1:10)
tumor_breast <- FindClusters(tumor_breast, resolution = 0.5)

benign_icc <- FindNeighbors(benign_icc, dims = 1:12)
benign_icc <- FindClusters(benign_icc, resolution = 0.5)

tumor_icc <- FindNeighbors(tumor_icc, dims = 1:9)
tumor_icc <- FindClusters(tumor_icc, resolution = 0.5)


# Run non-linear dimensional reduction (UMAP/tSNE) ####
benign_breast <- RunUMAP(benign_breast, dims = 1:12)
DimPlot(benign_breast, reduction = "umap")

tumor_breast <- RunUMAP(tumor_breast, dims = 1:10)
DimPlot(tumor_breast, reduction = "umap")

benign_icc <- RunUMAP(benign_icc, dims = 1:12)
DimPlot(benign_icc, reduction = "umap")

tumor_icc <- RunUMAP(tumor_icc, dims = 1:9)
DimPlot(tumor_icc, reduction = "umap")

#Save/Export Seurat Objects ####

saveRDS(benign_breast, file = "../main_rna_project/seurat_objects/clustered/bbclustered.rds")
saveRDS(tumor_breast, file = "../main_rna_project/seurat_objects/clustered/tbclustered.rds")
saveRDS(benign_icc, file = "../main_rna_project/seurat_objects/clustered/bicc_clustered.rds")
saveRDS(tumor_icc, file = "../main_rna_project/seurat_objects/clustered/ticc_clustered.rds")















