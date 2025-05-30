# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/dechapter.html
# https://github.com/JesperGrud/snRNAseq_eWAT/blob/main/Scripts/Main_01_Overall.md
# upregulated / downregulated adj.p.val<0.05 and logfc ><0 --- https://rdrr.io/bioc/limma/src/R/decidetests.R
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

#create col_data column
# Diet this var is true i.e 1 when its rep1 so up/down wrt rep1
# control--rep1, experiment-- rep2 
Limma_de_genes <- function(eWAT,col_data,rep1,rep2){
  eWAT <- eWAT[,eWAT@meta.data[[col_data]] == rep1 | eWAT@meta.data[[col_data]] == rep2]
  Exprs <- eWAT@assays$RNA@data
  Exprs <- as.matrix(Exprs)
  MD <- eWAT@meta.data
  MD$Cell <- rownames(MD)
  Mat <- data.frame(Cell = colnames(Exprs))
  Mat <- merge(Mat, MD, by = "Cell", sort = F)
  Diet <- factor(as.character(as.numeric(Mat[[col_data]] == rep1)), levels = c("0", "1"))
  Design <- model.matrix(~ Diet)
  Fit <- lmFit(Exprs, Design)
  
  # Computing Emperical Bayes statistics for each replicate
  Fit <- eBayes(Fit, trend = T)
  
  # Differential testing in each replicate
  Result <- topTable(Fit, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "P")
  colnames(Result)[c(1, 5)] <- c("logFC_Diet", "Padj_Diet")
  # summary(decideTests(Fit[,"Diet1"]))
  return(Result)
}

Limma_de_genes <- function(eWAT,col_data,rep1,rep2){
  filter_rep1 <- rep(FALSE,dim(eWAT@meta.data)[1])
  for (rep in rep2){
    filter_rep1 <- filter_rep1 | eWAT@meta.data[[col_data]] == rep1
  }
  filter_rep2 <- rep(FALSE,dim(eWAT@meta.data)[1])
  for (rep in rep2){
    filter_rep2 <- filter_rep2 | eWAT@meta.data[[col_data]] == rep2
  }
    eWAT <- eWAT[, filter_rep1| filter_rep2]
  Exprs <- eWAT@assays$RNA@data
  Exprs <- as.matrix(Exprs)
  MD <- eWAT@meta.data
  MD$Cell <- rownames(MD)
  Mat <- data.frame(Cell = colnames(Exprs))
  Mat <- merge(Mat, MD, by = "Cell", sort = F)
  Diet <- factor(as.character(as.numeric(filter_rep1)), levels = c("0", "1"))
  Design <- model.matrix(~ Diet)
  Fit <- lmFit(Exprs, Design)
  
  # Computing Emperical Bayes statistics for each replicate
  Fit <- eBayes(Fit, trend = T)
  
  # Differential testing in each replicate
  Result <- topTable(Fit, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "P")
  colnames(Result)[c(1, 5)] <- c("logFC_Diet", "Padj_Diet")
  # summary(decideTests(Fit[,"Diet1"]))
  return(Result)
}

Remove_high_and_mito_genes <- function(eWAT){
  high_exp_genes <- c("Chga", "Gcg", "Gm42418", "Iapp", "Ins1", "Ins2", "Malati", "Meg3", "mt-Atp6", "mt-col", "mt-C02", "mt-C03", "mt-cytb", "mt-Ndi", "mt-Nd2", "mt-Nd4",  "Pcsk2", "Ppy", "Pyy", "sc92")
  eWAT <- eWAT[!(rownames(eWAT) %in% high_exp_genes),]
  eWAT <- eWAT[!str_detect(rownames(eWAT),"mt-"),]  
  return(eWAT)
}

Merge_RC_HFD <- function(RC,HFD){
  eWAT <- merge(HFD, y = RC)
  eWAT <- eWAT[rownames(eWAT) %in% intersect(rownames(HFD),rownames(RC)),]
  return(eWAT)
}

RC <- readRDS("/home/krushna/Documents/Rstudio/adipose_result/Immune/RC_immune.Rds")
RC$Diet_Week <- paste0(RC$Diet,"-",RC$Week) #create that column
HFD <- readRDS("/home/krushna/Documents/Rstudio/adipose_result/Immune/HFD_immune.Rds")
HFD$Diet_Week <- paste0(HFD$Diet,"-",HFD$Week) #create that column

eWAT <- Merge_RC_HFD(RC,HFD)
eWAT <- Remove_high_and_mito_genes(eWAT)
top_de_genes <- list()
pair_list <- list(list("RC-8W","HFD-8W"),list("RC-8W","HFD-16W"),list("RC-8W","HFD-24W"),
                  list("RC-14W","HFD-8W"),list("RC-22W","HFD-16W"),list("RC-30W","HFD-24W"),
                  list("RC-8W","RC-14W"),list("RC-8W","RC-22W"),list("RC-8W","RC-30W"),
                  list("HFD-8W","HFD-16W"),list("HFD-8W","HFD-24W"))
cluster_name_list <- list("Macrophage", "B cells", "T cells") #as.list(as.list(levels(eWAT$seurat_clusters)))
for (clus in cluster_name_list){
  eWAT_sub <- subset(x = eWAT, subset = cluster_name == clus)
  print(clus)
  for (pair in pair_list){
    print(paste(pair[[1]],pair[[2]]))
    Result <- Limma_de_genes(eWAT_sub,'Diet_Week',pair[[1]],pair[[2]])
    top5 <- rbind(Result %>% top_n(n = 5, wt = logFC_Diet),Result %>% top_n(n = -5, wt = logFC_Diet))
    # top200 <- rbind(Result %>% top_n(n = 100, wt = logFC_Diet),Result %>% top_n(n = -100, wt = logFC_Diet))
    write.csv(Result,paste0("Immune_",clus,"_",str_replace_all(toString(pair[[1]]),", ","_"),"_VS_",str_replace_all(toString(pair[[2]]),", ","_"),".csv"))
    top_de_genes <- c(top_de_genes,rownames(top5))
  }
}



Averages <- data.frame(gene = unlist(top_de_genes))
pair_unlist <- unlist(pair_list)
for (i in 1:nrow(Averages)) {
  for (j in 2:(length(pair_unlist)+1)){
    print(paste(j,j-1))
    Averages[i,j] <-  log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1],
                                                             colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Week == pair_unlist[[j-1]],])])))
    colnames(Averages)[j] <- pair_unlist[[j-1]]
    
  }
}


rownames(Averages) <- make.names(Averages[,1], unique = TRUE)

temp <- t(scale(t(Averages[,c(2:ncol(Averages))])))
# Plot the heatmap
Heatmap(temp, na_col = "black", cluster_columns=F,
        row_order = rownames(temp) , 
        column_order = colnames(temp))


