library(SingleCellExperiment) 
library(simpleSingleCell) 
library(umap) 
library(Seurat) 
library(Matrix) 
library(scran) 
library(scater) 
library(harmony)
library(patchwork)
library(reticulate)

# "/home/krushna/Documents/RStudio/adipose_result/HFD_8W_1R_HFD_8W_2R_HFD_8W_3R_HFD_16W_1R_HFD_16W_2R_HFD_16W_3R_HFD_24W_1R_HFD_24W_2R_HFD_24W_3R/2021_08_26_22_48_08/eWAT_annotated_preprocessed.Rds"
# "/home/krushna/Documents/RStudio/adipose_result/RC_8W_1R_RC_8W_2R_RC_8W_3R_RC_14W_1R_RC_14W_2R_RC_14W_3R_RC_22W_1R_RC_22W_2R_RC_22W_3R_RC_30W_1R_RC_30W_2R_RC_30W_3R/2021_08_26_23_13_06/eWAT_annotated_preprocessed.Rds"
RDS_file_name  = "/home/krushna/Documents/RStudio/adipose_result/HFD_8W_1R_HFD_8W_2R_HFD_8W_3R_HFD_16W_1R_HFD_16W_2R_HFD_16W_3R_HFD_24W_1R_HFD_24W_2R_HFD_24W_3R/2021_08_26_22_48_08/eWAT_annotated_preprocessed.Rds"

folder = paste0(dirname(RDS_file_name),'/')
sink(paste0(folder,'Integration.txt'))
pdf(file = paste0(folder,'Integration.pdf'))

eWAT <- readRDS(RDS_file_name)
# eWAT <- subset(x = eWAT, subset = cluster_name == "Immune")

eWAT <- RenameAssays(object = eWAT, originalexp = 'RNA')
eWAT <- FindVariableFeatures(eWAT, nfeature=2000)
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony", reduction.name = 'harmony_umap', reduction.key = 'HARMONYUMAP_')
#eWAT <- RunTSNE(eWAT, dims=1:20, reduction="harmony")

eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 0.05, algorithm = 1) 
DimPlot(eWAT, group.by="Replicate", reduction = "harmony_umap") + plot_annotation(title = 'With_harmony')
DimPlot(eWAT, group.by="Week", reduction = "harmony_umap") + plot_annotation(title = 'With_harmony')
DimPlot(eWAT, group.by='seurat_clusters',reduction = "harmony_umap", label=TRUE)  + plot_annotation(title = 'With_harmony')

for (week_val in unique(eWAT$Week)){
  print(week_val)
  print(DimPlot(subset(eWAT, Week == week_val), reduction = "harmony_umap") + plot_annotation(title = paste('With_harmony',week_val)))
}
#0.005 -> 15
# no of points in each clusters...
print("number of points in cluster with harmony")
print("cluster : points")
for (x in 0:(length(unique(eWAT$seurat_clusters))-1)) {
  print(sprintf("%d : %d ",x,sum(eWAT$seurat_clusters == x)))
}

# saveRDS(eWAT, paste0(folder,'eWAT_with_embeddings.Rds'))

sink()
dev.off()

pdf(file = paste0(folder,'Multiresolution.pdf'))
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
for (res in c(0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5)){
  eWAT <- FindClusters(eWAT, resolution = res) ;
  print(DimPlot(eWAT, group.by='seurat_clusters', reduction = 'harmony_umap') + plot_annotation(title = res));
}
dev.off()
#subset(x = eWAT, subset = seurat_clusters == 16 & eWAT@reductions[["harmony_umap"]]@cell.embeddings[,1] <5)