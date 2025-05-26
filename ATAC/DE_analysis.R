library(dplyr)
library(patchwork)
library(Signac)
library(Seurat)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

#see the RNA translated data
DefaultAssay(integrated_data) <- 'RNA'
FeaturePlot(
  object = integrated_data,
  features = c('Ppy', 'Sst', 'Gcg', 'Ins1','Ins2')
)

# create a new UMAP using the integrated embeddings
integrated_data <- FindNeighbors(object = integrated_data, reduction = 'integrated_lsi', dims = 2:30)
integrated_data <- FindClusters(object = integrated_data, resolution=1, verbose = FALSE, algorithm = 3)
integrated_data <- RunUMAP(integrated_data, reduction = "integrated_lsi", dims = 2:30)

DimPlot(integrated_data, group.by = "seurat_clusters")

#DE analysis
# change back to working with peaks instead of gene activities
DefaultAssay(integrated_data) <- 'peaks'
da_peaks <- FindAllMarkers(
  object = integrated_data,
  # ident.1 = 0,
  # ident.2 = 1,
  only.pos = TRUE,
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nFeature_peaks'
)
head(da_peaks)


plot1 <- VlnPlot(
  object = integrated_data,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = unique(Idents(integrated_data))
)
plot2 <- FeaturePlot(
  object = integrated_data,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)
plot1 | plot2

# convert the regions into genes
open_all <- rownames(da_peaks)
closest_genes_all <- ClosestFeature(integrated_data, regions = open_all)
markers = cbind(da_peaks, closest_genes_all)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# write.csv(top20,"path")


#check the DE genes for annotaion
DefaultAssay(integrated_data) = 'RNA'
integrated_data <- FindVariableFeatures(integrated_data, nfeature=dim(integrated_data)[1])
DotPlot(
  object = integrated_data,
  features = c(top20$gene_name[1:10])
)
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

for (x in 0:(length(unique(integrated_data$seurat_clusters))-1)) {
  gene_to_print <- as.list(top20[top20$cluster == x,]$gene_name)
  feature_plots_genes_dot(integrated_data,list(list(paste('top10 for cluster',x),gene_to_print)),'UMAP')
}




new.name <- c("Beta","Alpha","Unknown","Unknown","Unknown","Unknown","Immune")
old.name <-  integrated_data@meta.data[["seurat_clusters"]]
integrated_data@meta.data[["cluster_name"]] <- plyr::mapvalues(x = old.name, from = levels(old.name), to = new.name)

DimPlot(
  object = integrated_data,
  group.by = 'cluster_name',
  label = TRUE,
  repel = TRUE)+ ggtitle('scATAC-seq')

#saveRDS(integrated_data,path)

CoveragePlot(
  object = integrated_data,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)