# https://satijalab.org/seurat/archive/v3.1/de_vignette.html
library(Seurat) 
library(scran) 
library(scater) 
library(dplyr)
library(patchwork)


feature_plots_genes <- function(eWAT,markers_genes_to_plot,REDUCTION){
  for (ab in markers_genes_to_plot){
    i<-1
    genes <- ab[[2]]
    while (i+4 <= length(ab[[2]])){
      values<-combine(genes[i:(i+4)])
      try({plt <- FeaturePlot(eWAT, features = values, reduction = REDUCTION) + DimPlot(eWAT, group.by="seurat_clusters",reduction = REDUCTION) + plot_annotation(ab[[1]]);
      print(plt)})
      i <- i+5
    }
    if (i<=length(ab[[2]])){
      values<-combine(genes[i:length(ab[[2]])])
      try({plt <- FeaturePlot(eWAT, features = values, reduction = REDUCTION) + DimPlot(eWAT, group.by="seurat_clusters",reduction = REDUCTION) + plot_annotation(ab[[1]]);
      print(plt)})
    }
  }
}


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

RDS_file_name  = "/home/krushna/Documents/RStudio/adipose_result/HFD_8W_1R_HFD_8W_2R_HFD_8W_3R_HFD_16W_1R_HFD_16W_2R_HFD_16W_3R_HFD_24W_1R_HFD_24W_2R_HFD_24W_3R/2021_08_26_22_48_08/eWAT_annotated_embedded.Rds"

folder = paste0(dirname(RDS_file_name),'/')
sink(paste0(folder,'DE_analysis.txt'))
pdf(file = paste0(folder,'DE_analysis.pdf'))

eWAT <- readRDS(RDS_file_name)
# eWAT <- subset(x = eWAT, subset = cluster_name == "Immune")
# eWAT <- subset(x = eWAT, subset = seurat_clusters %in% c("0","1","2","3","4","5","6","8","9","11","12","13"))

REDUCTION <- 'harmony_umap'

# eWAT <- RenameAssays(object = eWAT, originalexp = 'RNA')
# eWAT <- FindVariableFeatures(eWAT, nfeature=2000)
# eWAT <- ScaleData(eWAT)
# eWAT <- RunPCA(eWAT)
# eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
# eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony", reduction.name = 'harmony_umap', reduction.key = 'HARMONYUMAP_')

eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 0.05, algorithm = 1) 

pbmc.markers <- FindAllMarkers(eWAT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(pbmc.markers, paste0(folder,'markers.Rds'))
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top50 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(top50,paste0(folder,'top50.csv'), row.names = FALSE)
DoHeatmap(eWAT, features = top20$gene) + plot_annotation('Heap map of top 10 markers over all clusters')

#a 2D list first is the name and second is teh gene name
markers_genes_to_plot <- list(
  list('Alpha markers',list('Gcg','Irx2','Neurod1','Pcsk2','Irx2','Loxl4','Dpp4','Gc')),
  list('Beta markers', list('Ins','Ins1','Ins2','Igf1r','Pdx1','Nkx6-1','Nkx2-2','Iapp','Foxo1','Slc2a2','Gcgr','Ucn3','Pcks1','Pcks2','Herpud1','Hspa5')),
  list('gamma markers', list('Ppy','Sertm1','Spock1','Slitrk6','Meis2','Etv1')),
  list('immune markers', list('Flt3','Cd244a','Cd209a','Ccr7','Ltb','Ccl17', 
                              'Cd19','Igkc','Ighm','Cd79a', 'Cd3e','Trbc2','Cd3e','Cd4')),
  list('delta markers', list('Sst','Cckbr','Ghsr','Rbp4','Gabrb3','Casr','Ffar4','Gpr120')),
  list('epsilon markers', list('Ghrl','Oprk1','Asgr1')),
  list('Ductal markers',list('Hnf1b')),
  list('endothelial markers', list('Plvap','Esm1','Flt1','Plvap','Pecam1','Cd34')),
  list('Acinar markers', list('Cpa','Col1a1', 'Prss2','Prss3','Try4','Try5','Try10', 'Try10', 'Cel', 'Cpd1', 'Cela2a', 'Ctrb1')),
  list('Macrophage markers',list('C1qa','Emr1','Adgre1','Lyz2','Apoe','Ccl3','Ly6a','Atf3','Siglech')),
  list('Dendritic cell markers', list('Flt3','Cd244a','Cd209a','Ccr7','Ltb','Ccl17')),
  list('B cell', list('Cd19','Igkc','Ighm','Cd79a')),
  list('T cell', list('Cd3e','Trbc2','Cd3e','Cd4')),
  list('NK', list('Eomes','Ncr1','Klra8','Ly-49h','Klra1'))
  
)

feature_plots_genes(eWAT,markers_genes_to_plot,REDUCTION)
# feature_plots_genes_dot(eWAT,markers_genes_to_plot,REDUCTION)

#top 10 for each cluster
for (x in 0:(length(unique(eWAT$seurat_clusters))-1)) {
  gene_to_print <- as.list(top20[top20$cluster == x,]$gene)
  feature_plots_genes_dot(eWAT,list(list(paste('top10 for cluster',x),gene_to_print)),REDUCTION)
}



sink()
dev.off()


#rename_clusters
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 0.06, algorithm = 1) 
new.name <- c("Endocrine","Endocrine","Endothelial","Endocrine","Immune","Stromal","Progenitor","Acinar","Immune","Fibroblast")
#new.name <- c("Endocrine","Endocrine","Endocrine","Endothelial","Progenitor","Immune","Stromal","Acinar","Immune","Erythrocyte")
old.name <-  eWAT@meta.data[["seurat_clusters"]]
eWAT@meta.data[["cluster_name"]] <- plyr::mapvalues(x = old.name, from = levels(old.name), to = new.name)
DimPlot(eWAT, group.by="cluster_name",reduction = REDUCTION)
for (week_val in unique(eWAT$Week)){
  print(week_val)
  print(DimPlot(subset(eWAT, Week == week_val),group.by="cluster_name", reduction = REDUCTION) + plot_annotation(title = week_val))
}

# saveRDS(eWAT,"/home/krushna/Documents/RStudio/adipose_result/HFD/1810_Ins2_also_removed/eWAT_with_embeddings.Rds")