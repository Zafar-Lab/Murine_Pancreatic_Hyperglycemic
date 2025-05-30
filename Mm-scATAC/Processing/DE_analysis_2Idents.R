library(stringr)
library(dplyr)
library(patchwork)
library(Signac)
library(Seurat)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
set.seed(1234)
# DE analysis of alpha atac
#subset & setting Idents
integrated_data = readRDS("/home/krushna/Documents/RStudio/Atac/Endocrine/Alpha/Alpha_atac.Rds")
DimPlot(integrated_data, group.by = "seurat_clusters")
integrated_data$Diet_Week =  paste0(sapply(strsplit(integrated_data$dataset,"-",), `[`, 1),'-',sapply(strsplit(integrated_data$dataset,"-",), `[`, 2))


pair_list <- list(list("RC-8W","HFD-8W"),list("RC-8W","HFD-24W"), list("RC-8W","RC-14W"),list("RC-8W","RC-30W"),
                  list("RC-14W","RC-30W"),list("RC-14W","HFD-8W"),list("RC-14W","HFD-24W"),
                  list("RC-30W","HFD-8W"),list("RC-30W","HFD-24W"),
                  list("HFD-8W","HFD-24W"))

#,list("RC-8W","HFD-16W") list("HFD-8W","HFD-16W"),,list("RC-22W","HFD-16W"),list("RC-8W","RC-22W")
integrated_data_all = integrated_data
integrated_data <- subset(x = integrated_data_all, subset = seurat_clusters == "Alpha-1")
Idents(integrated_data) = integrated_data$Diet_Week

for (pair in pair_list){
  print(paste(pair[[1]],pair[[2]]))
  DefaultAssay(integrated_data) <- 'peaks'
  da_peaks <- FindMarkers(
    object = integrated_data,
    ident.1 = pair[[1]],
    ident.2 = pair[[2]],
    # only.pos = TRUE,
    min.pct = 0.05,
    test.use = 'LR',
    logfc.threshold = 0.1,
    latent.vars = 'nFeature_peaks'
  )
  # convert the regions into genes
  open_all <- rownames(da_peaks)
  closest_genes_all <- ClosestFeature(integrated_data, regions = open_all)
  markers = cbind(da_peaks, closest_genes_all)
  # top20 <- rbind(markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC),markers %>% group_by(cluster) %>% top_n(n = -50, wt = avg_log2FC))
  write.csv(markers,paste0("Alpha_Atac_",str_replace_all(toString(pair[[1]]),"-","_"),"_VS_",str_replace_all(toString(pair[[2]]),"-","_"),".csv"))
}

da_peaks = read.table("/home/krushna/Documents/RStudio/Atac/Endocrine/Alpha/Alplha-1_coverage/Sc_ATAC_ALPHA_1_24w.csv", header=F, sep=",",row.names = 1)
for (region in rownames(da_peaks)){
  png(file= paste0("/home/krushna/Documents/RStudio/Atac/Endocrine/Alpha/Alplha-1_coverage/",region,".png"),
      width=1200, height=600)
  print(CoveragePlot(
    object = integrated_data,
    region = region,
    extend.upstream = 40000,
    extend.downstream = 20000
  ))
  dev.off()
}
