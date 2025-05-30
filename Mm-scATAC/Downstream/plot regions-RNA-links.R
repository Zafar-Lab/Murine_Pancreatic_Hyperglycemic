library(Seurat)
library(dplyr)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
plot_regions <- function(Mm.ATAC, genes_to_plot,  save_path){
  Mm.ATAC <- LinkPeaks(
    object = Mm.ATAC,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = genes_to_plot
  )
  
  pdf(save_path)
  for (gene in genes_to_plot){
    tryCatch({ coverage_plot <- CoveragePlot(
      object = Mm.ATAC,
      region = gene,
      expression.assay = "RNA",
      extend.upstream = 5000,
      extend.downstream = 5000,
    )
    # features = gene,
    
    coverage_plot[[1]][[1]] <- coverage_plot[[1]][[1]] + 
      theme(panel.grid.major.y = element_line(color = "grey80"))
    coverage_plot <- coverage_plot + plot_annotation(title = paste(gene))
    print(
      coverage_plot
    )}, error = function(e){print(paste(gene, "gene not found"))})
    
  }
  dev.off()
  
}

Mm.ATAC <- readRDS("./Mm.ATAC-Alpha.RDS")
Mm.ATAC <- RegionStats(Mm.ATAC, genome = BSgenome.Mmusculus.UCSC.mm10)

genes_to_plot = c("Prrc2c", "Neurod1", "Itpr1", "Foxa2", "Atf3", "Golgb1", "Zdbf2", "Pax6", "Foxa2", "Arx", "Mafb")
plot_regions(Mm.ATAC,  genes_to_plot,  "ROI-RNA-Alpha1.pdf")

genes_to_plot = c("Rpl37a", "Rpl38", "Rps27", "Smim22", "Bola2", "Arf5", "Rplp1")
plot_regions(Mm.ATAC,  genes_to_plot,  "ROI-RNA-Alpha2.pdf")

Mm.ATAC_a1 <- subset(Mm.ATAC, subset = predicted.id == "Alpha-1")
Idents(Mm.ATAC_a1) <- Mm.ATAC_a1$diet.week

ICT <- c("Tomm7", "Sec61g", "Mrln", "Tmed9", "Rbp4") #intercellular transport
plot_regions(Mm.ATAC_a1,  ICT,  "ROI-RNA-Alpha1-ICT.pdf")

OP <- c("Cox7a2l", "Cox4i1", "Chchd2") #Oxidative phosphorylation
plot_regions(Mm.ATAC_a1,  OP,  "ROI-RNA-Alpha1-Oxidative_phosphorylation.pdf")

response_os <- c("Rack1", "Hsph1", "Chchd2", "Rps3") #Response to oxidative stress
plot_regions(Mm.ATAC_a1,  response_os,  "ROI-RNA-Alpha1-Response_to_os.pdf")

genes_to_plot = c("Ins1", "Ins2", "Ppy","Sst", "Gcg")
plot_regions(Mm.ATAC,  genes_to_plot,  "ROI-RNA.pdf")

