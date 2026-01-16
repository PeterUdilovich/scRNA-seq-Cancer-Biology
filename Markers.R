library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

benign_breast <- readRDS("../main_rna_project/seurat_objects/clustered/bbclustered.rds")
tumor_breast <- readRDS("../main_rna_project/seurat_objects/clustered/tbclustered.rds")
benign_icc <- readRDS("../main_rna_project/seurat_objects/clustered/bicc_clustered.rds")
tumor_icc <- readRDS("../main_rna_project/seurat_objects/clustered/ticc_clustered.rds")

#Finding differentially expressed features (cluster biomarkers) ####

benign_breast.markers <- FindAllMarkers(benign_breast, only.pos = TRUE)
tumor_breast.markers <- FindAllMarkers(tumor_breast, only.pos = TRUE)
benign_icc.markers <- FindAllMarkers(benign_icc, only.pos = TRUE)
tumor_icc.markers <- FindAllMarkers(tumor_icc, only.pos = TRUE)

#Find top 12 genes for each cluster ####

#bb
top12_listbb <- benign_breast.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 12, with_ties = FALSE) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

top12_listbb

#tb
top12_listtb <- tumor_breast.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 12, with_ties = FALSE) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

top12_listtb

#bicc
top12_listbicc <- benign_icc.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 12, with_ties = FALSE) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

top12_listbicc

#ticc
top12_listticc <- tumor_icc.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 12, with_ties = FALSE) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

top12_listticc




#Assigning cell type identity to clusters ####

#bb
bb.new.cluster.ids <- c("Neuronal", "Epithelial (luminal/secretory)", "Blood vascular endothelial", "Epithelial (basal)", "T cells", "Myeloid APC",
                     "Unassigned", "Probably doublets", "Fibroblast", "Unassigned (cycling)", "Lymphatic endothelial")
names(bb.new.cluster.ids) <- levels(benign_breast)
benign_breast <- RenameIdents(benign_breast, bb.new.cluster.ids)
DimPlot(benign_breast, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#tb
tb.new.cluster.ids <- c("Cytotoxic T", "Malignant epithelial", "Malignant epithelial", "Cancer-associated fibroblasts", "Blood vascular endothelial", "Tumor-associated macrophages",
                        "Unassigned")
names(tb.new.cluster.ids) <- levels(tumor_breast)
tumor_breast <- RenameIdents(tumor_breast, tb.new.cluster.ids)
DimPlot(tumor_breast, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#bicc

bicc.new.cluster.ids <- c("CD4 T", "Cytotoxic lymphocytes", "NK cells", "Unassigned", "Neutrophils", "Unassigned (cycling)",
                        "Dendritic cells", "Epithelial", "Unassigned (cycling)", "cDC1 dendritic cells", "B cells", "Blood vascular endothelial", "Plasmacytoid DC", "Hepatocyte-like contamination")
names(bicc.new.cluster.ids) <- levels(benign_icc)
benign_icc <- RenameIdents(benign_icc, bicc.new.cluster.ids)
DimPlot(benign_icc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#ticc

ticc.new.cluster.ids <- c("Fibroblasts", "Malignant epithelial (mucinous/secretory)", "Tumor-associated macrophages", "Malignant epithelial", "Dendritic cells")
names(ticc.new.cluster.ids) <- levels(tumor_icc)
tumor_icc <- RenameIdents(tumor_icc, ticc.new.cluster.ids)
DimPlot(tumor_icc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#Double Plots ####

p1 <- DimPlot(benign_breast, reduction="umap", label=TRUE, repel=TRUE, pt.size=0.5) + NoLegend() +
  coord_cartesian(clip="off") + theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

p2 <- DimPlot(tumor_breast, reduction="umap", label=TRUE, repel=TRUE, pt.size=0.5) + NoLegend() +
  coord_cartesian(clip="off") + theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

p1 | p2



p3 <- DimPlot(benign_icc, reduction="umap", label=TRUE, repel=TRUE, pt.size=0.5) + NoLegend() +
  coord_cartesian(clip="off") + theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

p4 <- DimPlot(tumor_icc, reduction="umap", label=TRUE, repel=TRUE, pt.size=0.5) + NoLegend() +
  coord_cartesian(clip="off") + theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

p3 | p4




