# including libraries.....
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(BiocSingular) 
library(umap) 
library(Seurat) 
library(Matrix) 
library(scran) 
library(scater) 
library(DropletUtils) 
library(batchelor)
library(harmony)
library(ComplexHeatmap)
library(circlize)
library(MAST)
library(limma)
library(RANN)
library(biomaRt)
library(kBET)
library(lisi)
library(org.Mm.eg.db)
library(dplyr)
library(clusterProfiler)
library(rWikiPathways)
library(GEOquery)
library(edgeR)
library(glmnet)
library(velocyto.R)
library(phateR)
library(ElPiGraph.R)
library(slingshot)
library(TSCAN)
library(monocle3)
library(tradeSeq)
library(seriation)
library(patchwork)
# library(scCATCH)
# library(SoupX)
# library(scDblFinder)
library(celda)



Read_Data_Sets <- function(path="./Desktop/Projects/Project_DS/islets/",dataset_names){
  datasets = c()
  for (name in dataset_names){
    print(paste("reading....\n",path,name,sep=''))
    temp <- Read10X_h5(paste(path,name,sep=''), use.names = TRUE, unique.features = TRUE)
    datasets <- c(datasets,temp)
  }
  print('dimentions of the dataset after reading data')
  for (temp in datasets){
    print(dim(temp))
  }
  return(datasets)
}  

# Reduce_Data <- function(datasets){
#   reduced = c()
#   for (data in datasets){
#     reduced <- c(reduced,data[1:3000,1:3000])
#   }
#   return(reduced)
# }

Bring_Same_Dimentions <- function(datasets){
  # Find all non-zero genes in all conditions
  print('bring same rows')
  Genes <- c()
  for (data in datasets){
    Genes <- c(Genes,rownames(data))
  }
  Genes <- unique(Genes)
  samedim = c()
  # Insert empty lines with missing genes to get same dimensions
  for (data in datasets){
    Tmp <- as.sparse(matrix(ncol=ncol(data), nrow=length(Genes[!(Genes %in% rownames(data))]), data=0))
    rownames(Tmp) <- Genes[!(Genes %in% rownames(data))]
    data <- Matrix::rbind2(data, Tmp)
    data <- data[ order(rownames(data)),]
    samedim <- c(samedim,data)
    print(dim(data))
  }
  return(samedim)
  
}

Make_Cell_Names_Unique <- function(datasets,names){
  unique_cells <- c()
  for (i in 1:length(datasets)){
    colnames(datasets[[i]]) <- paste(names[[i]],"_",colnames(datasets[[i]]), sep="")
    unique_cells <-c(unique_cells,datasets[[i]])
  }
  return(unique_cells)
}

Create_Sce <- function(datasets){
  sce <- list()
  for (data in datasets){
    temp <- SingleCellExperiment(list(counts=data))
    sce <- c(sce,temp)
  }
  return(sce)
}

Add_QC_metric <- function(datasets_sce){
  qc_added <- list()
  for (data in datasets_sce){
    temp <- addPerCellQC(data,subsets=list(Mito=grep(pattern = "^mt-", x = rownames(data), value = FALSE)))
    qc_added <- c(qc_added,temp)
  }
  return(qc_added)
}


#check -> https://rdrr.io/bioc/scuttle/f/vignettes/overview.Rmd (todo)
# doing something wrong much cells are going away so using new
Threshold_Droplets <- function(datasets_sce){
  thesh_drop <- list()
  print('Threshold Droplets')
  for (data in datasets_sce){
    data <- data[,!(data$subsets_Mito_percent > 10)] 
    # data <- data[,!(data$sum/data$detected > 2.5)]
    data <- data[,(data$sum >= 1000 & data$detected >= 500)]
    thesh_drop <- c(thesh_drop,data)
    print(dim(data))
  }
  return(thesh_drop)
}

# http://www.bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/overview.html -- autodetect outlier
Threshold_Droplets_new <- function(datasets_sce){
  thesh_drop <- list()
  print('Threshold Droplets new')
  for (data in datasets_sce){
    is.mito <- grep(pattern = "^mt-", x = rownames(data), value = FALSE)
    per.cell <- perCellQCMetrics(data, subsets=list(Mito=is.mito))
    qc.stats <- quickPerCellQC(per.cell, percent_subsets="subsets_Mito_percent")
    colSums(as.matrix(qc.stats))
    filtered <- data[,!qc.stats$discard]
    print(dim(filtered))
    thesh_drop <- c(thesh_drop,filtered)
  }
  return(thesh_drop)
}

Threshold_Genes <- function(datasets_sce){
  low = list()
  for (data in datasets_sce){
    #low <- c(low,names(which(nexprs(data, byrow=T) <= 10)))
    low <- c(low,list(names(which(nexprs(data, byrow = T) <= 10))))
  }
  Low <- Reduce(intersect, low)
  print("low genes are: ")
  print(length(Low))
  filtered = list()
  print('Threshold Genes')
  for (data in datasets_sce){
    temp <- data[ which(!(rownames(data) %in% Low)),] 
    filtered <- c(filtered,temp)
    print(dim(temp))
  }
  return(filtered)
} 

# http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html
Normalize_The_Count <- function(datasets_sce){
  normalized <- list()
  for (data in datasets_sce){
    clus <- quickCluster(data, use.ranks=FALSE, BSPARAM=IrlbaParam())
    temp <- computeSumFactors(data, min.mean=0.1, cluster=clus)
    normalized <- c(normalized,logNormCounts(temp))
  }
  return(normalized)
}

Doublet_Scores <- function(datasets_sce){
  doublet <- list()
  for (data in datasets_sce){
    data$DoubletScore <- scDblFinder::computeDoubletDensity(data, BSPARAM=IrlbaParam())  
    doublet <- c(doublet,data)
  }
  return(doublet)
}

Doublet_Threshold <- function(datasets_sce){
  doublet <- list()
  for (data in datasets_sce){ 
    outlier <- isOutlier(data$DoubletScore,type = 'higher')
    doublets_removed <- data[,!outlier]
    print(dim(doublets_removed))
    doublet <- c(doublet,doublets_removed)
  }
  return(doublet)
}

Create_seurat <- function(datasets_sce){
  surat <- list()
  for(data in datasets_sce){
    temp <- as.Seurat(data, counts = "counts", data = "logcounts")
    surat <- c(surat,temp)
  }
  return(surat)
}


ExprsFun <- function(x) { sum(x > 0) }

Find_Ambient <- function(datasets,datasets_sce){
  empty <- list()
  for (data in datasets){
    Empty <- data[, colnames(data) %in% names(which(Matrix::colSums(data) >= 1 & Matrix::colSums(data) <= 10)),] # Find droplets that are empty, but has recovered genes (less than 10 UMI counts, but more than 1)
    Empty <- Empty[ rownames(Empty) %in% names(which(Matrix::rowSums(Empty) >= 50)),] # Find genes with at least 50 UMI counts in 'empty' droplets
    dim(Empty)
    Empty <- as.matrix(t(Empty))
    Empty <- apply(Empty,2,ExprsFun)/nrow(Empty)
    empty <- c(empty,sort(names(which(Empty >= 0.001)))) # Find genes expressed in at least 0.1% of the empty droplets. 
  }
  Amb <- unique(empty) 
  Amb <- Amb[Amb %in% rownames(datasets_sce[[1]])] 
  
  
  return(Amb)
}

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets
Remove_Empty_Droplets <- function(datasets_sce){
  removed_ed <- list()
  i <- 0 
  print('Remove empty droplets')
  for (data in datasets_sce){
    # br.out <- barcodeRanks(counts(data))
    # # Making a plot.
    # plot(br.out$rank, br.out$total, log="xy", xlab=paste("Rank",as.character(i),sep = ' '), ylab="Total")
    # o <- order(br.out$rank)
    # lines(br.out$rank[o], br.out$fitted[o], col="red")
    # 
    # abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    # abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    # legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    #        legend=c("knee", "inflection"))
    # print(metadata(br.out)$knee)
    
    #e.out <- emptyDrops(counts(data), lower=metadata(br.out)$knee) #need to set lower value
    e.out <- emptyDrops(counts(data), lower = 2000)
    temp <- data[,which(e.out$FDR <= 0.05)] # earlier value 0.05, 0.09, 0.15 (not changing)
    # 0.01
    print(dim(temp))
    removed_ed <- c(removed_ed,temp)
    i<-i+1
  }
  return(removed_ed)
}


file_name_extension_create_dir <- function(replicate_names){
  date_time <- toString(Sys.time())
  date_time <- chartr(old = "- :", new = "___", date_time)
  date_time <- paste0(date_time)
  # temp <- chartr(old= '-',new='_',paste(replicate_names,collapse =  '_'))
  temp <- strsplit(replicate_names[1],"-")[[1]][1]
  file_name <- paste0('adipose_result/',temp,'/',date_time,'/')
  dir.create(file_name,recursive = TRUE)
  return(file_name)
}


# required variables to update
path = "/home/krushna/Documents/adipose-tissue-analysis/our dataset/"
dataset_names = c(
  "HFD-8W-1R_filtered_feature_bc_matrix.h5",
  "HFD-8W-2R_filtered_feature_bc_matrix.h5",
  "HFD-8W-3R_filtered_feature_bc_matrix.h5",
  "HFD-16W-1R_filtered_feature_bc_matrix.h5",
  "HFD-16W-2R_filtered_feature_bc_matrix.h5",
  "HFD-16W-3R_filtered_feature_bc_matrix.h5",
  "HFD-24W-1R_filtered_feature_bc_matrix.h5",
  "HFD-24W-2R_filtered_feature_bc_matrix.h5",
  "HFD-24W-3R_filtered_feature_bc_matrix.h5"
)
# "RC-8W-1R_filtered_feature_bc_matrix.h5",
# "RC-8W-2R_filtered_feature_bc_matrix.h5",
# "RC-8W-3R_filtered_feature_bc_matrix.h5",
# "RC-14W-1R_filtered_feature_bc_matrix.h5",
# "RC-14W-2R_filtered_feature_bc_matrix.h5",
# "RC-14W-3R_filtered_feature_bc_matrix.h5",
# "RC-22W-1R_filtered_feature_bc_matrix.h5",
# "RC-22W-2R_filtered_feature_bc_matrix.h5",
# "RC-22W-3R_filtered_feature_bc_matrix.h5",
# "RC-30W-1R_filtered_feature_bc_matrix.h5",
# "RC-30W-2R_filtered_feature_bc_matrix.h5",
# "RC-30W-3R_filtered_feature_bc_matrix.h5"



replicate_names <- sapply(strsplit(dataset_names,"_"), `[`, 1)
file_name <- file_name_extension_create_dir(replicate_names)


#saving the output
sink(paste0(file_name,'decontxv2_preprocessed.txt'))
pdf(file = paste0(file_name,'decontxv2_preprocessed.pdf'))



replicates <- Read_Data_Sets(path,dataset_names)

# Full data is needed.
# replicates <- Reduce_Data(replicates)

## Fill in the matrices to give them all the same dimensions
## RATIONALE: Genes with 0 counts across all barcodes in a particular experiment are left out from zUMIs.
replicates <- Bring_Same_Dimentions(replicates)

## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
replicates <- Make_Cell_Names_Unique(replicates,replicate_names)
# replicates_cnu_cp <- replicates

## Create SingleCellExperiment objects to use DropUtils, scran, scater, etc.
replicates_sce <- Create_Sce(replicates)
replicates_sce_t1 <- replicates_sce
# replicates_sce_cre_cp <- replicates_sce
# AJ: 16 may
replicates_sce <- Remove_Empty_Droplets(replicates_sce)

#todo check for correctness of perCellQCMetrics as its updated or downgrade scater
## Calculate QC parameters (throws a warning, that can be ignored)
# https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/overview.html
# https://bioconductor.org/books/release/OSCA/quality-control.html
replicates_sce <- Add_QC_metric(replicates_sce)

## REMOVE: Droplets with more than 15% mitochondrial reads, less than 1000 UMIs, 
#less than 500 genes or extremely high ratio between counts and genes (low complexity))
replicates_sce <- Threshold_Droplets_new(replicates_sce)

## Threshold filtering of genes in each dataset
## REMOVE: Genes expressed in less than 10 nuclei in all datasets
replicates_sce <- Threshold_Genes(replicates_sce)

req_just_pre = TRUE
if (req_just_pre){
  replicates_sce_ <- Normalize_The_Count(replicates_sce)
  replicates_seurat_ <- Create_seurat(batchelor::multiBatchNorm(replicates_sce_))
  eWAT_ <- merge(replicates_seurat_[[1]], y=replicates_seurat_[2:length(replicates_seurat_)])
  diet_week_replicate <- sapply(strsplit(rownames(eWAT_@meta.data),"_"),`[`,1)
  eWAT_$Dataset <- diet_week_replicate
  eWAT_$Replicate <- sapply(strsplit(diet_week_replicate,"-"),`[`,3)
  eWAT_$Week <- sapply(strsplit(diet_week_replicate,"-"),`[`,2)
  eWAT_$Diet <- sapply(strsplit(diet_week_replicate,"-"),`[`,1)
  saveRDS(eWAT_, paste0(file_name,'eWAT_preprocessed.Rds'))
}

finalList <- list()
i <- 0
print ('decontX step:')
for (data in replicates_sce){
  print (i)
  #delta = c(10, 20), estimateDelta = FALSE
  df = decontX(data,delta = c(10, 20), estimateDelta = FALSE)
  print ('Contaimination is')
  print (mean(df$decontX_contamination))
  assayNames(df)<-c('count.org','counts')
  finalList <- c(finalList, df)
  i <- i + 1
}

replicates_sce <- finalList
## Detect ambient genes ## RATIONALE: Ambient genes are detected in droplets that do not contain nuclei.
Amb <- Find_Ambient(replicates,replicates_sce)
saveRDS(Amb, paste0(file_name,"Ambient.Rds"))
# Normalize the counts
replicates_sce <- Normalize_The_Count(replicates_sce)

# AJ: Step is akipped: Memory error is coming
## Calculate doublet scores. NOTE: THIS STEP IS NON-DETERMINISTI - RESULTS VARY FROM RUN TO RUN
# replicates_sce <- Doublet_Scores(replicates_sce)
# replicates_sce_DP_cp <- replicates_sce
# replicates_sce <- Doublet_Threshold(replicates_sce)
### QC by clustering - Individual datasets. NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
### RATIONaLE: Low quality nuclei may be included due to threshold effects. Deep clustering can help to reveal if there are groups of low quality nuclei that does not mix with the remaining nuclei, and thus can be removed.

## Create Seurat objects
replicates_seurat <- Create_seurat(batchelor::multiBatchNorm(replicates_sce))

# Merge the datasets
eWAT <- merge(replicates_seurat[[1]], y=replicates_seurat[2:length(replicates_seurat)])
diet_week_replicate <- sapply(strsplit(rownames(eWAT@meta.data),"_"),`[`,1)
eWAT$Dataset <- diet_week_replicate
eWAT$Replicate <- sapply(strsplit(diet_week_replicate,"-"),`[`,3)
eWAT$Week <- sapply(strsplit(diet_week_replicate,"-"),`[`,2)
eWAT$Diet <- sapply(strsplit(diet_week_replicate,"-"),`[`,1)

saveRDS(eWAT, paste0(file_name,'decontX_preprocessed.Rds'))
# eWAT <- readRDS("/home/krushna/Documents/RStudio/adipose_result/HFD/decontxv2/decontX_preprocessed.Rds")
data_seu <- eWAT
data <- as.SingleCellExperiment(data_seu)
gene_index_list <- list(
  Gcg=grep(pattern = "^Gcg$", x = rownames(data), value = FALSE),
  Ins1 = grep(pattern = "^Ins1$", x = rownames(data), value = FALSE),
  Ins2 = grep(pattern = "^Ins2$", x = rownames(data), value = FALSE),
  Ppy = grep(pattern = "^Ppy$", x = rownames(data), value = FALSE),
  Sst = grep(pattern = "^Sst$", x = rownames(data), value = FALSE)
)

data@colData@listData[["sum"]] <- NULL
data@colData@listData[["detected"]] <- NULL
data@colData@listData[["total"]] <- NULL

data <- addPerCellQC(data,subsets=gene_index_list)
gene_list = list(list('Gcg'),list("Ins1","Ins2"),list("Sst"),list("Ppy") )
noofcells = length(colnames(data))
thrange <- c(5,10,15,20)
for (threshold in thrange){
  print(threshold);
  for(i in 1:length(gene_list)){
    for (g in gene_list[[i]]){
      values <- data@colData@listData[[paste0("subsets_",g,"_percent")]];
      print(paste(g,length(values[values>threshold])));
    } 
  }
   
  
  for (i in 1:length(gene_list)){
    for (j in 1:length(gene_list[[i]])){
      for (i1 in i:length(gene_list)){
        if(i1==i){next}
        for (j1 in 1:length(gene_list[[i1]])){
          g1 = gene_list[[i]][j]
          g2 = gene_list[[i1]][j1]
          values1 <- data@colData@listData[[paste0("subsets_",g1,"_percent")]]
          values2 <- data@colData@listData[[paste0("subsets_",g2,"_percent")]]
          print(paste(g1,g2,length(values1[values1>threshold & values2>threshold]) ))
          data@colData@listData[[paste0(g1,g2,threshold)]] <- as.integer(values1>threshold & values2>threshold)
        }
      }
    }
  }
}
filter_cell = rep(c(FALSE),each=noofcells)
threshold_to_filter = 5
for (i in 1:length(gene_list)){
  for (j in 1:length(gene_list[[i]])){
    for (i1 in i:length(gene_list)){
      if(i1==i){next}
      for (j1 in 1:length(gene_list[[i1]])){
        g1 = gene_list[[i]][j]
        g2 = gene_list[[i1]][j1]

        values1 <- data@colData@listData[[paste0("subsets_",g1,"_percent")]]
        values2 <- data@colData@listData[[paste0("subsets_",g2,"_percent")]]
        filter_cell <- filter_cell | (values1>threshold_to_filter & values2>threshold_to_filter)
      }
    }
  }
}
eWAT <- as.Seurat(data, counts = "counts", data = "logcounts")
eWAT <- RenameAssays(object = eWAT, originalexp = 'RNA')
eWAT <- FindVariableFeatures(eWAT, nfeature=2000)
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony", reduction.name = 'harmony_umap', reduction.key = 'HARMONYUMAP_')
REDUCTION = 'harmony_umap'

for (g in gene_list){
  print(FeaturePlot(eWAT, features =  c(paste0("subsets_",g,"_percent")), reduction = REDUCTION) + plot_annotation(g));
  print(FeaturePlot(eWAT,features = c(g), reduction = REDUCTION));
} 

eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony") 
eWAT <- FindClusters(eWAT, resolution = 0.05, algorithm = 1) 
DimPlot(eWAT, group.by='seurat_clusters',reduction = "harmony_umap")   + plot_annotation(title = 'With_harmony')

for (threshold in thrange){
  for (i in 1:size){
    g1 = gene_list[i]
    for(j in i:size){
      if(i==j){next}
      g2 = gene_list[j]
      print(FeaturePlot(eWAT,features = paste0(g1,g2,threshold), reduction = REDUCTION))
    }
  }
}


saveRDS(data_seu[,!filter_cell], paste0(file_name,'decontX_multi_gene_removed_preprocess.Rds'))


sink()
dev.off()

# ay = 'logcounts'
# vioplot(temp@assays@data@listData[[ay]]['Ins1',]/colSums(temp@assays@data@listData[[ay]]))
