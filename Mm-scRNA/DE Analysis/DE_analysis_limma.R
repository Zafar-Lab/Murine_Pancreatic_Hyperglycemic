# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/dechapter.html
# https://github.com/JesperGrud/snRNAseq_eWAT/blob/main/Scripts/Main_01_Overall.md
library(edgeR)
library(Seurat) 
library(scran) 
library(scater) 
library(patchwork)
library(limma)
library(MAST)
library(harmony)
library(ComplexHeatmap)
library(dplyr)
library(stringr)

HFD <- readRDS("/home/krushna/Documents/RStudio/adipose_result/HFD/1810_Ins2_also_removed/eWAT_with_embeddings.Rds")
RC <- readRDS("/home/krushna/Documents/RStudio/adipose_result/RC/1810_Ins2_also_removed/eWAT_with_embeddings.Rds")

HFD <- subset(x = HFD, subset = cluster_name == "Immune")
RC <- subset(x = RC, subset = cluster_name == "Immune")
eWAT <- merge(HFD, y = RC)

# restrictive genes (high avg values)
high_exp_genes <- c("Chga", "Gcg", "Gm42418", "Iapp", "Ins1", "Ins2", "Malati", "Meg3", "mt-Atp6", "mt-col", "mt-C02", "mt-C03", "mt-cytb", "mt-Ndi", "mt-Nd2", "mt-Nd4",  "Pcsk2", "Ppy", "Pyy", "sc92")
eWAT <- eWAT[!(rownames(eWAT) %in% high_exp_genes),]
eWAT <- eWAT[!str_detect(rownames(eWAT),"mt-"),]  #remove mito genes
rm(HFD)
rm(RC)
gc() #free unused R memory


Exprs <- eWAT@assays$RNA@data
Exprs <- as.matrix(Exprs)


# Setting up design matrices for differential testing of diet effects
MD <- eWAT@meta.data
MD$Cell <- rownames(MD)
Mat <- data.frame(Cell = colnames(Exprs))
Mat <- merge(Mat, MD, by = "Cell", sort = F)
Diet <- factor(as.character(as.numeric(Mat$Diet == "HFD")), levels = c("0", "1"))
Design <- model.matrix(~ Diet)

# Fitting a linear model for each replicate
Fit <- lmFit(Exprs, Design)

# Computing Emperical Bayes statistics for each replicate
Fit <- eBayes(Fit, trend = T)

# Differential testing in each replicate
Result <- topTable(Fit, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")

# Setting column names
colnames(Result)[c(1, 5)] <- c("logFC_Diet", "Padj_Diet")

# Merging results from each replicate
# ResultFrame <- Result
# 
# # Annotating diet effects
# ResultFrame$Diet <- "NS"
# ResultFrame[ResultFrame$logFC_Diet > 0 & ResultFrame$Padj_Diet <= 0.05, "Diet"] <- "Enriched"
# ResultFrame[ResultFrame$logFC_Diet < 0 & ResultFrame$Padj_Diet <= 0.05 , "Diet"] <- "Exclusive"
# colnames(ResultFrame)[1] <- "Symbol"

# Re-ordering columns
# ResultFrame <- ResultFrame[, c(1, 6, 2, 3, 4, 5)]
#top5 and bottom5
top5 <- rbind(Result %>% top_n(n = 10, wt = logFC_Diet),Result %>% top_n(n = -10, wt = logFC_Diet))

tb50 <- rbind(Result %>% top_n(n = 25, wt = logFC_Diet),Result %>% top_n(n = -25, wt = logFC_Diet))
write.csv(tb50,'/home/krushna/Documents/RStudio/adipose_result/Top_HFD_vs_RC/updown_Immune.csv', row.names = TRUE)


Averages <- data.frame(gene = rownames(top5))
for (i in 1:nrow(Averages)) {
  Averages[i,2] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1],
                                                          colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Diet == "RC",])])))
  Averages[i,3] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1],
                                                          colnames(eWAT@assays$RNA@data) %in%  rownames(eWAT@meta.data[ eWAT@meta.data$Diet == "HFD",])])))
}
rownames(Averages) <- Averages[,1]
colnames(Averages)[2:3] <- c("RC","HFD")
# temp <- t(scale(t(Averages[,c(2,3)])))
temp <- Averages[,c(2,3)]
Heatmap(temp, cluster_columns=F)

# Clust1 <- ResultFrame
# 
# 
# Clust1 <- Clust1[Clust1$Marker != "NS", ]
# 


# 
# eWAT <- FindVariableFeatures(eWAT, nfeature=2000)
# eWAT <- ScaleData(eWAT)
# eWAT <- RunPCA(eWAT)
# eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
# eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony", reduction.name = 'harmony_umap', reduction.key = 'HARMONYUMAP_')
# #eWAT <- RunTSNE(eWAT, dims=1:20, reduction="harmony")
# 
# eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
# eWAT <- FindClusters(eWAT, resolution = 0.05, algorithm = 1) #eWAT <- RunTSNE(eWAT, dims = 1:20, reduction = 'pca')6
# 
# DimPlot(eWAT, reduction = "harmony_umap", label=T)
# DimPlot(eWAT , reduction = "harmony_umap",  label=T, split.by="Diet")
# DimPlot(eWAT, label=T,reduction = "harmony_umap", split.by="Replicate")
# DimPlot(eWAT, group.by="Replicate", reduction = "harmony_umap")
# DimPlot(eWAT, group.by="Diet", reduction = "harmony_umap")
# 
# 
# 
# 
# 
# 
Endocrine <- data.frame(gene = c("Chga","Chgb","Scg3","Syp"), type = "Endocrine")
Endothelial <- data.frame(gene = c("Plvap","Flt1","Esm1"), type = "Endothelial")
Immune <- data.frame(gene = c("H2-Eb1","Lyz2","C1qc"), type = "Immune")
Progenitor <- data.frame(gene = c("Mki67","Smc4","Ube2c"), type = "Progenitor")
Stromal <- data.frame(gene = c("Acta2","Myl9"), type = "Stromal")
Acinar <- data.frame(gene = c("Try4","Cela1","Try5"), type = "Acinar")
Erythrocyte <- data.frame(gene = c("Hba-a1","Hbb-bt","Alas2"), type = "Erythrocyte")

types <- c("Endocrine","Endothelial","Immune","Progenitor", "Stromal","Acinar","Erythrocyte")
# Combine results
Averages <- rbind(Endocrine,Endothelial,Immune,Progenitor, Stromal, Acinar,Erythrocyte)

# Calculate average expression for selected marker genes
datasets <- unique(eWAT@meta.data[["Week"]])
for (i in 1:nrow(Averages)) {
  for (j in 1:length(types)){
    for (k in 1:length(datasets)){
      Averages[i,3+length(datasets)*(j-1)+k-1] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$cluster_name == types[j] & eWAT@meta.data$Week == datasets[k],])])))
      print(paste(i,3+length(datasets)*(j-1)+k-1,types[j],"_",datasets[k]))
      colnames(Averages)[3+length(datasets)*(j-1)+k-1] = paste0(types[j],"_RC_",datasets[k])
      }
  
    }
}




rownames(Averages) <- Averages[,1]

# Exp1 <- Averages[,c(3,5,7,9,11,13,15,17,19,21)]
# 
# Exp2 <- Averages[,c(4,6,8,10,12,14,16,18,20,22)]
# Exp2 <- t(scale(t(Exp2)))
# Averages[,c(3,5,7,9,11,13,15,17,19,21)] <- Exp1
# Averages[,c(4,6,8,10,12,14,16,18,20,22)] <- Exp2

temp <- t(scale(t(Averages[,c(3:ncol(Averages))])))
# Plot the heatmap
Heatmap(temp, na_col = "black", cluster_columns=F,
        row_order = rownames(temp) , 
        column_order = colnames(temp))
