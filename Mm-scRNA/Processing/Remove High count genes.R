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

feature_plots_genes <- function(eWAT,markers_genes_to_plot,umap_embedings){
  for (ab in markers_genes_to_plot){
    i<-1
    genes <- ab[[2]]
    while (i+4 <= length(ab[[2]])){
      values<-combine(genes[i:(i+4)])
      try({plt <- FeaturePlot(eWAT, features = values,reduction = umap_embedings) + DimPlot(eWAT, group.by="seurat_clusters",reduction = umap_embedings) + plot_annotation(ab[[1]]);
      print(plt)})
      i <- i+5
    }
    if (i<=length(ab[[2]])){
      values<-combine(genes[i:length(ab[[2]])])
      try({plt <- FeaturePlot(eWAT, features = values,reduction = umap_embedings) + DimPlot(eWAT, group.by="seurat_clusters",reduction = umap_embedings) + plot_annotation(ab[[1]]);
      print(plt)})
    }
  }
}

RunSCVI <- function(eWAT){
  
  top2000 <- head(VariableFeatures(eWAT), 2000)
  eWAT2000 <- eWAT[top2000]
  print(eWAT2000)
  
  
  sc <- import('scanpy', convert = FALSE)
  scvi <- import('scvi', convert = FALSE)
  scvi$settings$progress_bar_style = 'tqdm'
  adata <- sc$AnnData(
    X   = t(as.matrix(GetAssayData(eWAT2000,slot='counts'))), #scVI requires raw counts
    obs = eWAT2000[[]],
    var = GetAssay(eWAT2000)[[]]
  )
  print(adata)
  
  scvi$data$setup_anndata(adata,batch_key = 'Dataset')
  
  # create the model
  model = scvi$model$SCVI(adata, use_cuda = TRUE)
  
  # train the model
  model$train()
  
  
  latent = model$get_latent_representation()
  
  # put it back in our original Seurat object
  latent <- as.matrix(latent)
  rownames(latent) = colnames(eWAT2000)
  eWAT[['scvi']] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(eWAT))
  return(eWAT)
}

markers_genes_to_plot <- list(
  list('Alpha markers',list('Gcg','Irx2','Arx','Neurod1','Ttr','Pcsk2','Mafb')),
  list('Beta markers', list('Ins','Ins1','Ins2','Igf1r','Pdx1','Nkx6-1','Nkx2-2','Iapp', 'Gcgr','Hspa5')),
  list('gamma markers', list('Ppy','Sertm1','Spock1','Abcc9','Slitrk6','Meis2','Etv1','Id4')),
  list('delta markers', list('Sst','Rbp4')),
  list('endothelial markers', list('Flt1','Plvap','Cd36'))
)



RDS_file_name  = "/home/krushna/Documents/RStudio/adipose_result/HFD_8W_1R_HFD_8W_2R_HFD_8W_3R_HFD_16W_1R_HFD_16W_2R_HFD_16W_3R_HFD_24W_1R_HFD_24W_2R_HFD_24W_3R/2021_08_26_22_48_08/eWAT_annotated_preprocessed.Rds"
folder = paste0(dirname(RDS_file_name),'/')
sink(paste0(folder,'remove high count.txt'),append = T)
pdf(file = paste0(folder,'remove high count.pdf'))

eWAT <- readRDS(RDS_file_name)


eWAT_sce <- as.SingleCellExperiment(eWAT)
eWAT_sce <- addPerFeatureQC(eWAT_sce)
# rownames(eWAT_sce[eWAT_sce@rowRanges@elementMetadata$mean>50])
# eWAT_sce@rowRanges@elementMetadata$mean[eWAT_sce@rowRanges@elementMetadata$mean>50]


fil_genes <- rownames(eWAT)[which(!(rownames(eWAT) %in% rownames(eWAT_sce[eWAT_sce@rowRanges@elementMetadata$mean>50])))]
eWAT_fil <- eWAT[fil_genes,]


eWAT_fil <- FindVariableFeatures(eWAT_fil, nfeature=2000)
eWAT_fil <- ScaleData(eWAT_fil)
eWAT_fil <- RunPCA(eWAT_fil)
eWAT_fil <- RunUMAP(eWAT_fil, dims=1:20, reduction="pca", reduction.name = 'pca_umap', reduction.key = 'PCAUMAP_')
#eWAT <- RunTSNE(eWAT, dims = 1:20, reduction = 'pca')6
eWAT_fil <- FindNeighbors(object = eWAT_fil, dims = 1:20, reduction = "pca")
eWAT_fil <- FindClusters(eWAT_fil, resolution = 0.05, algorithm = 1) 

DimPlot(eWAT_fil, group.by="Replicate", reduction="pca_umap")  + plot_annotation(title = 'Without_harmony')
DimPlot(eWAT_fil, group.by="Week", reduction="pca_umap")  + plot_annotation(title = 'Without_harmony')
DimPlot(eWAT_fil, group.by='seurat_clusters', reduction="pca_umap") + plot_annotation(title = 'Without_harmony')
# no of points in each clusters...
print("number of points in cluster without harmony")
print("cluster : points")
for (x in 0:(length(unique(eWAT_fil$seurat_clusters))-1)) {
  print(sprintf("%d : %d ",x,sum(eWAT_fil$seurat_clusters == x)))
}
eWAT@reductions[["pca_umap_removed"]] <- eWAT_fil@reductions[["pca_umap"]]
eWAT@meta.data[["seurat_clusters"]] <- eWAT_fil@meta.data[["seurat_clusters"]]
feature_plots_genes(eWAT,markers_genes_to_plot,"pca_umap_removed")


#with harmony
### Final embedding and clustering using ambient gene removal and Harmony 
# eWAT <- FindVariableFeatures(eWAT, nfeature=2000)
# eWAT <- ScaleData(eWAT)
# eWAT <- RunPCA(eWAT)
eWAT_fil <- RunHarmony(eWAT_fil, group.by.vars="Dataset", dims.use=1:20)
eWAT_fil <- RunUMAP(eWAT_fil, dims=1:20, reduction="harmony", reduction.name = 'harmony_umap', reduction.key = 'HARMONYUMAP_')
#eWAT <- RunTSNE(eWAT, dims=1:20, reduction="harmony")

eWAT_fil <- FindNeighbors(object = eWAT_fil, dims = 1:20, reduction = "harmony")
eWAT_fil <- FindClusters(eWAT_fil, resolution = 0.05, algorithm = 1) 
DimPlot(eWAT_fil, group.by="Replicate", reduction = "harmony_umap") + plot_annotation(title = 'With_harmony')
DimPlot(eWAT_fil, group.by="Week", reduction = "harmony_umap") + plot_annotation(title = 'With_harmony')
DimPlot(eWAT_fil, group.by='seurat_clusters',reduction = "harmony_umap")  + plot_annotation(title = 'With_harmony')

#0.005 -> 15
# no of points in each clusters...
print("number of points in cluster with harmony")
print("cluster : points")
for (x in 0:(length(unique(eWAT_fil$seurat_clusters))-1)) {
  print(sprintf("%d : %d ",x,sum(eWAT_fil$seurat_clusters == x)))
}

eWAT@reductions[["harmony_umap_removed"]] <- eWAT_fil@reductions[["harmony_umap"]]
eWAT@meta.data[["seurat_clusters"]] <- eWAT_fil@meta.data[["seurat_clusters"]]
feature_plots_genes(eWAT,markers_genes_to_plot,"harmony_umap_removed")

#scvi
eWAT_fil <- RunSCVI(eWAT_fil)
eWAT_fil <- FindNeighbors(eWAT_fil, dims = 1:10, reduction = 'scvi')
eWAT_fil <- FindClusters(eWAT_fil, resolution = 0.05)

eWAT_fil <- RunUMAP(eWAT_fil, dims = 1:10, reduction = 'scvi', reduction.name = 'scvi_umap', reduction.key = 'SCVIUMAP_')

DimPlot(eWAT_fil, group.by="Replicate", reduction = "scvi_umap") + plot_annotation(title = 'With_scvi')
DimPlot(eWAT_fil, group.by="Week", reduction = "scvi_umap") + plot_annotation(title = 'With_scvi')
DimPlot(eWAT_fil, group.by='seurat_clusters',reduction = "scvi_umap")  + plot_annotation(title = 'With_scvi')
eWAT@reductions[["scvi_umap_removed"]] <- eWAT_fil@reductions[["scvi_umap"]]
eWAT@meta.data[["seurat_clusters"]] <- eWAT_fil@meta.data[["seurat_clusters"]]
feature_plots_genes(eWAT,markers_genes_to_plot,"scvi_umap_removed")


#seurat

pancreas.list <- SplitObject(eWAT_fil, split.by = "Dataset")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

reference.list <- pancreas.list[unique(eWAT_fil$Dataset)] 
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:20)

DefaultAssay(pancreas.integrated) <- "integrated"

pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:20)

pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:20, reduction = "pca")
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.05, algorithm = 1) 

DimPlot(pancreas.integrated , group.by="Dataset") + plot_annotation(title = 'With_Seurat')
DimPlot(pancreas.integrated , group.by="Replicate") + plot_annotation(title = 'With_Seurat')
DimPlot(pancreas.integrated , group.by='seurat_clusters')  + plot_annotation(title = 'With_Seurat')
eWAT@reductions[["seurat_umap_removed"]] <- pancreas.integrated@reductions[["umap"]]
eWAT@meta.data[["seurat_clusters"]] <- pancreas.integrated@meta.data[["seurat_clusters"]]
feature_plots_genes(eWAT,markers_genes_to_plot,"seurat_umap_removed")
# FeaturePlot(pancreas.integrated , features = c('Apoe'), reduction = 'umap')# + DimPlot(eWAT, group.by="seurat_clusters") 
# FeaturePlot(pancreas.integrated , features = c('Sox9'), reduction = 'umap')# + DimPlot(eWAT, group.by="seurat_clusters") 
# FeaturePlot(pancreas.integrated , features = c('Acta2'), reduction = 'umap')# + DimPlot(eWAT, group.by="seurat_clusters") 
# FeaturePlot(pancreas.integrated , features = c('Cd34'), reduction = 'umap')# + DimPlot(eWAT, group.by="seurat_clusters") 
# FeaturePlot(pancreas.integrated , features = c('Nkx6-1'), reduction = 'umap')# + DimPlot(eWAT, group.by="seurat_clusters") 
# FeaturePlot(pancreas.integrated , features = c('Flt1'), reduction = 'umap')# + DimPlot(eWAT, group.by="seurat_clusters") 


sink()
dev.off()
