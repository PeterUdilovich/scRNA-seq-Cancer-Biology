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


# Load the PBMC dataset
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


#Seurat Objects
benign_breast <- CreateSeuratObject(counts = benign_breast.data, project = "bb", min.cells = 3, min.features = 200)
tumor_breast <- CreateSeuratObject(counts = tumor_breast.data, project = "tb", min.cells = 3, min.features = 200)
benign_icc <- CreateSeuratObject(counts = benign_icc.data, project = "bicc", min.cells = 3, min.features = 200)
tumor_icc <- CreateSeuratObject(counts = tumor_icc.data, project = "ticc", min.cells = 3, min.features = 200)

benign_breast[["percent.mt"]] <- PercentageFeatureSet(benign_breast, pattern = "^MT-")
tumor_breast[["percent.mt"]] <- PercentageFeatureSet(tumor_breast, pattern = "^MT-")
benign_icc[["percent.mt"]] <- PercentageFeatureSet(benign_icc, pattern = "^MT-")
tumor_icc[["percent.mt"]] <- PercentageFeatureSet(tumor_icc, pattern = "^MT-")

VlnPlot(benign_breast, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(tumor_breast, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(benign_icc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(tumor_icc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

benign_breast <- subset(benign_breast, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
benign_breast <- NormalizeData(benign_breast)

benign_breast <- FindVariableFeatures(benign_breast, selection.method = "vst", nfeatures = 2000)

"""benign_breast <- FindVariableFeatures(benign_breast, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(benign_breast), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(benign_breast)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2"""



all.genes <- rownames(benign_breast)
benign_breast <- ScaleData(benign_breast, features = all.genes)



benign_breast <- RunPCA(benign_breast, features = VariableFeatures(object = benign_breast))
print(benign_breast[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(benign_breast, dims = 1:2, reduction = "pca")

DimPlot(benign_breast, reduction = "pca") + NoLegend()

benign_breast <- RunUMAP(benign_breast, dims = 1:10)
DimPlot(benign_breast, reduction = "umap")














