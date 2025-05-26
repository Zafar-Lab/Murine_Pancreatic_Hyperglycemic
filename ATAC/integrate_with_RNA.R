library(dplyr)
library(patchwork)
library(Signac)
library(Seurat)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
set.seed(1234)


integrated_data <- readRDS("/home/krushna/Documents/RStudio/Atac/All/integrated_data.Rds")
integrated_data = subset(integrated_data, predicted.id_1 == "Endocrine" )
data_rna <- readRDS("/home/krushna/Documents/RStudio/All/HFD_RC_combined.Rds")

DefaultAssay(integrated_data) = 'RNA'
transfer.anchors <- FindTransferAnchors(
  reference = data_rna,
  query = integrated_data,
  reduction = 'cca',
  features = rownames(integrated_data)
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = data_rna$cluster_name,
  weight.reduction = integrated_data[['integrated_lsi']],
  dims = 2:30
)
predicted.labels2 <- TransferData(
  anchorset = transfer.anchors,
  refdata = data_rna$cluster_name_rename,
  weight.reduction = integrated_data[['integrated_lsi']],
  dims = 2:30
)
colnames(predicted.labels) = paste0(colnames(predicted.labels),"_1")
colnames(predicted.labels2) = paste0(colnames(predicted.labels2),"_2")

integrated_data <- AddMetaData(object = integrated_data, metadata = predicted.labels)
integrated_data <- AddMetaData(object = integrated_data, metadata = predicted.labels2)

plot1 <- DimPlot(
  object = data_rna,
  group.by = 'cluster_name',
  label = TRUE,
  repel = TRUE,
  reduction = 'harmony_umap') + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = integrated_data,
  group.by = 'predicted.id_1',
  label = TRUE,
  repel = TRUE)+ ggtitle('scATAC-seq')

plot1 + plot2
#saveRDS(integrated_data,path)




# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(data_rna)
refdata <- GetAssayData(data_rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = integrated_data[["integrated_lsi"]],
                           dims = 2:30)
integrated_data[["RNA_imputed"]] <- imputation

DefaultAssay(integrated_data) <- "RNA_imputed"
coembed <- merge(x = data_rna, y = integrated_data)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "cluster_rename"))
