#ref https://yoseflab.github.io/VISION/articles/VISION-vignette.html
library(Seurat)
library(SeuratDisk)
library(VISION)
library(ggplot2)
library(ComplexHeatmap)
data = readRDS("/home/krushna/Documents/adipose-tissue-analysis/vision/data/RC_HFD_Endocrine_Alpha_after_merging_reremoval.Rds")
data = subset(data,subset = cluster_name == 'Mm-Alpha 1')
data = subset(data,subset = Diet_Week != "HFD_24W")


DefaultAssay(data)
data@meta.data = data@meta.data[, c("cluster_name", "Replicate", "Week", "Diet", "Diet_Week")]

options(mc.cores = 50)

signatures = list.files(path="/home/krushna/Documents/adipose-tissue-analysis/vision/data/signatures", pattern=NULL, all.files=FALSE,full.names=TRUE)
vision.obj <- Vision(data, signatures = signatures, latentSpace=data@reductions[["harmony"]]@cell.embeddings, dimRed='harmony_umap')
vision.obj <- analyze(vision.obj)
viewResults(vision.obj)

sigScores <- getSignatureScores(vision.obj)[rownames(data@meta.data), ]
data@meta.data = cbind(data@meta.data, sigScores)
result <- apply(sigScores, 2, function(x) (x-min(x))/(max(x)-min(x))) # scaling to 0-1
result <- aggregate(result, by = list(group = data@meta.data$Diet_Week), FUN = mean)
row_names = result[,1]
result <- result[,2:dim(result)[2]]
rownames(result) <- row_names
result = t(as.matrix(result))
Heatmap(result, row_names_gp = gpar(fontsize = 6), name = "Normalized Scores",
        cluster_columns = FALSE)
col_order = c("RC_8W", "RC_14W", "RC_22W", "HFD_8W", "HFD_16W")
Heatmap(result, row_names_gp = gpar(fontsize = 6), name = "Normalized Scores",
        cluster_columns = FALSE, column_order = col_order)


# DotPlot(data, features = colnames(sigScores)[2:5], scale = FALSE) & coord_flip()

# DimPlot(data, reduction = 'harmony_umap', group.by = 'cluster_name')


# # Display autocorrelation coefficients, p-values for signatures
# head(getSignatureAutocorrelation(vision.obj))


# Plot signature scores for a signature of interest
# umap <- getProjections(vision.obj)[["Seurat_harmony_umap"]]

# ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores) + geom_point()

# VlnPlot(data,features = "scores", group.by = 'cluster_name')