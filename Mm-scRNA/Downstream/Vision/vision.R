#ref https://yoseflab.github.io/VISION/articles/VISION-vignette.html
library(Seurat)
library(SeuratDisk)
library(VISION)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(stringr)
data = readRDS("./data/RC_HFD_MERGED_BETA.Rds")
# data = subset(data,subset = cluster_name == 'Mm-Alpha 1')
# data = subset(data,subset = Diet_Week != "HFD_24W")
# data = subset(data,subset = cluster_name == "Mm-Beta 1")
data = subset(data,subset = cluster_name != "Mm-Beta 5")
data = subset(data,subset = cluster_name != "Mm-Beta 4")
# data = subset(data,subset = Diet_Week != "HFD_24W")



DefaultAssay(data)
data@meta.data = data@meta.data[, c("cluster_name", "Replicate", "Week", "Diet", "Diet_Week")]

options(mc.cores = 50)

signatures = list.files(path="/home/krushna/Documents/adipose-tissue-analysis/misc/vision/data/signatures", pattern=NULL, all.files=FALSE,full.names=TRUE)
vision.obj <- Vision(data, signatures = signatures, latentSpace=data@reductions[["harmony"]]@cell.embeddings, dimRed='harmony_umap')
vision.obj <- analyze(vision.obj)
viewResults(vision.obj)
# saveRDS(vision.obj, "./results/Beta-clusters-vision.obj.RDS")
sigScores <- getSignatureScores(vision.obj)[rownames(data@meta.data), ]
data@meta.data = cbind(data@meta.data, sigScores)
result <- apply(sigScores, 2, function(x) (x-min(x))/(max(x)-min(x))) # scaling to 0-1
result <- aggregate(result, by = list(group = data@meta.data$cluster_name), FUN = mean)
# result <- aggregate(result, by = list(group = data@meta.data$Diet_Week), FUN = mean)

row_names = result[,1]
result <- result[,2:dim(result)[2]]
rownames(result) <- row_names
result = t(as.matrix(result))
# col_order = c("RC_8W", "RC_14W", "RC_22W", "RC_30W", "HFD_8W", "HFD_16W")
col_order = c("Mm-Beta 1","Mm-Beta 2","Mm-Beta 3")
rownames(result) <-sapply(rownames(result), function(x) gsub("(.{30})", "\\1\n", x) )

tiff("./results/Beta clsuterwise.tiff", width = 8, height = 15, units = "in", res = 300, compression = "lzw")
Heatmap(result, row_names_gp = gpar(fontsize = 8), name = "Normalized Scores",
        cluster_columns = FALSE, column_order = col_order,show_heatmap_legend = T,
        row_names_max_widt= unit(10, "cm"))

dev.off()
# heatmap_legend_param = list(at = c(0,75,150))
# DotPlot(data, features = colnames(sigScores)[2:5], scale = FALSE) & coord_flip()

# DimPlot(data, reduction = 'harmony_umap', group.by = 'cluster_name')


# # Display autocorrelation coefficients, p-values for signatures
# head(getSignatureAutocorrelation(vision.obj))

               
# Plot signature scores for a signature of interest
# umap <- getProjections(vision.obj)[["Seurat_harmony_umap"]]

# ggplot() + aes(x=umap[, 1], y=umap[, 2], color=sigScores) + geom_point()

# VlnPlot(data,features = "scores", group.by = 'cluster_name')

