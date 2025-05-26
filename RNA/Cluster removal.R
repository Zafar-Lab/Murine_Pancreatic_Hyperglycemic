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
library(Seurat) 
library(scran) 
library(scater) 
library(dplyr)
library(patchwork)


#reading the dataset
RDS_file_name  = "C:/Users/admin/Desktop/Islet Project/RDS File/Endothelial/RC_Endothelial.Rds"
# folder = paste0(dirname(RDS_file_name),'/')
# sink(paste0(folder,'Integration.txt'))
# pdf(file = paste0(folder,'Integration.pdf'))

eWAT <- readRDS(RDS_file_name)
# eWAT <- subset(x = eWAT, subset = cluster_name == "Immune")

#to get the embeddings
eWAT <- RenameAssays(object = eWAT, originalexp = 'RNA')
eWAT <- FindVariableFeatures(eWAT, nfeature=2000)
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony", reduction.name = 'harmony_umap', reduction.key = 'HARMONYUMAP_')
#eWAT <- RunTSNE(eWAT, dims=1:20, reduction="harmony")
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")

#define the resolution that you want
eWAT <- FindClusters(eWAT, resolution = 0.2, algorithm = 1) 
DimPlot(eWAT, group.by='seurat_clusters',reduction = "harmony_umap", label=TRUE)  + plot_annotation(title = 'With_harmony')

# to remove some of the cells based on cluster number and umap_axis
cells_to_remove <- subset(x = eWAT, subset = seurat_clusters == 4)
#& eWAT@reductions[["harmony_umap"]]@cell.embeddings[,1] <5)
saveRDS(cells_to_remove,"C:/Users/admin/Desktop/Islet Project/RDS File/Endothelial/RC_endocrine4")
#axis x eWAT@reductions[["harmony_umap"]]@cell.embeddings[,1]
#axis Y eWAT@reductions[["harmony_umap"]]@cell.embeddings[,2]
eWAT <- eWAT[,!colnames(eWAT) %in% colnames(cells_to_remove)]

#DE analysis
eWAT <- FindClusters(eWAT, resolution = 0.2, algorithm = 1) 
DimPlot(eWAT, group.by='seurat_clusters',reduction = "harmony_umap", label=TRUE)  + plot_annotation(title = 'With_harmony')
feature_plots_genes_dot <- function(eWAT,markers_genes_to_plot,REDUCTION){
  for (ab in markers_genes_to_plot){
    i<-1
    genes <- ab[[2]]
    while (i+9 <= length(ab[[2]])){
      values<-unlist(genes[i:(i+9)])
      try({plt <-DotPlot(eWAT, features = values) + RotatedAxis() +plot_annotation(ab[[1]]);
      print(plt)})
      i <- i+10
    }
    if (i<=length(ab[[2]])){
      values<-unlist(genes[i:length(ab[[2]])])
      try({plt <-DotPlot(eWAT, features = values) + RotatedAxis() +plot_annotation(ab[[1]]);
      print(plt)})
    }
  }
}
pbmc.markers <- FindAllMarkers(eWAT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(pbmc.markers, paste0(folder,'markers.Rds'))
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top50 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(top50,paste0(folder,'top50.csv'), row.names = FALSE)
DoHeatmap(eWAT, features = top20$gene) + plot_annotation('Heap map of top 10 markers over all clusters')
for (x in 0:(length(unique(eWAT$seurat_clusters))-1)) {
  gene_to_print <- as.list(top20[top20$cluster == x,]$gene)
  feature_plots_genes_dot(eWAT,list(list(paste('top10 for cluster',x),gene_to_print)),REDUCTION)
}