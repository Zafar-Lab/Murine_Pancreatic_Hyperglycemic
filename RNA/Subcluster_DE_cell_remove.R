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
RDS_file_name  = "/home/krushna/Documents/Rstudio/adipose_result/Immune/after_removal/Immune_HFD_after_removal_merged_reremoval.RDS"

folder = paste0(dirname(RDS_file_name),'/')
sink(paste0(folder,'Integration.txt'))
pdf(file = paste0(folder,'Integration.pdf'))

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
eWAT <- FindClusters(eWAT, resolution = 0.4, algorithm = 1) 
DimPlot(eWAT, group.by='seurat_clusters',reduction = "harmony_umap", label=TRUE)  + plot_annotation(title = 'With_harmony')

# to remove some of the cells based on cluster number and umap_axis
cells_to_remove <- subset(x = eWAT, subset = seurat_clusters == 16 & eWAT@reductions[["harmony_umap"]]@cell.embeddings[,1] <5)
#axis x eWAT@reductions[["harmony_umap"]]@cell.embeddings[,1]
#axis Y eWAT@reductions[["harmony_umap"]]@cell.embeddings[,2] 

eWAT <- eWAT[,!colnames(eWAT) %in% colnames(cells_to_remove)]


#DE analysis
# eWAT <- FindClusters(eWAT, resolution = 0.2, algorithm = 1) 
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
# saveRDS(pbmc.markers, paste0(folder,'markers.Rds'))
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top50 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
# write.csv(top50,paste0(folder,'top50.csv'), row.names = FALSE)
DoHeatmap(eWAT, features = top20$gene) + plot_annotation('Heap map of top 10 markers over all clusters')
for (x in 0:(length(unique(eWAT$seurat_clusters))-1)) {
  gene_to_print <- as.list(top20[top20$cluster == x,]$gene)
  feature_plots_genes_dot(eWAT,list(list(paste('top10 for cluster',x),gene_to_print)),REDUCTION)
}



#week wise plots
for (week_val in unique(eWAT$Week)){
  print(week_val)
  print(DimPlot(subset(eWAT, Week == week_val),group.by="seurat_clusters", reduction = 'harmony_umap') + plot_annotation(title = week_val))
}

print("number of points in cluster with harmony")
print("cluster : points")
for (x in 0:(length(unique(eWAT$seurat_clusters))-1)) {
  print(sprintf("%d : %d ",x,sum(eWAT$seurat_clusters == x)))
}


#rename the cluster
new.name <- c('Macrophage1','B cells','T cells1','Macrophage2','Macrophage3','Macrophage4','Macrophage5','T cells2','Macrophage6','NK cells','Macrophage7')
old.name <-  eWAT@meta.data[["seurat_clusters"]]
eWAT@meta.data[["cluster_name"]] <- plyr::mapvalues(x = old.name, from = levels(old.name), to = new.name)

DimPlot(eWAT, group.by='cluster_name',reduction = "harmony_umap", label=TRUE) 
for (week_val in unique(eWAT$Week)){
  print(week_val)
  print(DimPlot(subset(eWAT, Week == week_val),group.by="cluster_name", reduction = 'harmony_umap',label=TRUE) + plot_annotation(title = week_val))
}
print("number of points in cluster with harmony")
print("cluster : points")
cluster_names <- unique(eWAT$cluster_name)
for (x in cluster_names) {
  print(sprintf("%s : %d ",x,sum(eWAT$cluster_name == x)))
}


cluster_names <- unique(eWAT$cluster_name)
for (week_val in unique(eWAT$Week)){
  print(week_val)
  week_subset <- subset(eWAT, Week == week_val)
  for (x in cluster_names) {
    print(sprintf("%s : %d ",x,sum(week_subset$cluster_name == x)))
  }
}

# saveRDS(eWAT, paste0(folder,'eWAT_with_embeddings.Rds'))
# sink()
# FeaturePlot(eWAT, features = c('Ifit1'), reduction = 'harmony_umap')
dev.off()

# temp  = eWAT@assays[["RNA"]]@counts@Dimnames[[2]]
# temp[eWAT@assays[["RNA"]]@counts@Dimnames[[2]] %in% m1@assays[["RNA"]]@counts@Dimnames[[2]]]