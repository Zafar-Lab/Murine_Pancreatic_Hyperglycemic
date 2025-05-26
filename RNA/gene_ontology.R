library(clusterProfiler)
require(DOSE)
library(enrichplot)
library(org.Mm.eg.db)

library(readxl)
library(GOSemSim)
EG <- pbmc.markers

EG <- read.csv("/home/krushna/Documents/Rstudio/Immune_B cells_RC_8W_VS_HFD_8W_HFD_16W_HFD_24W.csv")
genes <- EG$X
EG1 <- enrichGO(gene          = genes,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05
                )
ds<-data.frame(EG1@result)
dotplot(EG1)


CV<-EG1@result
write.csv(CV,'D://TCGA BRCA new//MSC TNBC NON TNBC//Xcell//Final Results//Postitively up regeulated genes//Pathway Non_TNBC_Cluster_9_significant_positively_upregulated_genes.csv')
typeof(EG1)
dotplot(EG1, showCategory=10, split=".sign") + facet_grid(.~.sign)






original_gene_list <- EG$logFC_Diet
names(original_gene_list) <- EG$X
original_gene_list = sort(original_gene_list, decreasing = TRUE)
gse <- gseGO(geneList=original_gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
d <- godata('org.Mm.eg.db', ont="BP")
x2 <- pairwise_termsim(gse,method="Wang", semData = d) 
emapplot(x2)
ridgeplot(gse) + labs(x = "enrichment distribution")

library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=EG$X, pathway.id="dme04130", species = "dme")

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=EG$X, pathway.id="dme04130", species = "dme", kegg.native = F)


