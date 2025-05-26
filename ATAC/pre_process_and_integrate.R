#running on all data
library(dplyr)

library(patchwork)

library(Signac)
library(Seurat)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

Read_Dataset_With_Counts <-  function(counts_path,fragments_path){
  # low = list()
  counts <- Read10X_h5(filename = counts_path)
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = fragments_path
  )
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
  )
  print(pbmc)
  return(pbmc)
}

Read_Dataset_with_granges <- function(fragments_paths,genome_ranges){
  datasets = list()
  for(frag_path in fragments_paths){
    fragcounts <- CountFragments(fragments = frag_path)
    atac.cells <- fragcounts[fragcounts$frequency_count > 2000, "CB"]
    atac.frags <- CreateFragmentObject(path = frag_path, cells = atac.cells)
    counts <- FeatureMatrix(
      fragments = atac.frags,
      features = genome_ranges,
      cells = atac.cells
    )
    atac.assay <- CreateChromatinAssay(
      counts = counts,
      sep = c(":", "-"),
      genome = 'mm10',  #updated here
      min.features = 200,
      fragments = atac.frags
    )
    pbmc2 <- CreateSeuratObject(counts = atac.assay, assay = "peaks")
    print('before processing')
    print(pbmc2)
    
    pbmc2 <- FindTopFeatures(pbmc2 , min.cutoff = 10)
    pbmc2 <- RunTFIDF(pbmc2)
    pbmc2 <- RunSVD(pbmc2)
    
    temp = strsplit(frag_path,"/")[[1]]
    temp = temp[length(temp)]
    dataset_name = strsplit(temp,"_")[[1]][1]
    pbmc2$dataset = dataset_name
    print(dataset_name)
    print(VlnPlot(
      object = pbmc2,
      features = c('nCount_peaks', 'nFeature_peaks'),
      pt.size = 0.1,
      ncol = 5
    ) +  plot_annotation(title = dataset_name))
    pbmc2 <- subset(pbmc2, nCount_peaks > 2000 & nCount_peaks < 30000)
    print(pbmc2)
    datasets = c(datasets,pbmc2)
  }
  return(datasets)
}

Annotate_Datasets <- function(datasets){
  datasets_annotated = list()
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  # # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  
  for (dataset in datasets){
    print(unique(dataset$dataset))
    Annotation(dataset) <- annotations
    gene.activities <- GeneActivity(dataset)
    
    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    dataset[['RNA']] <- CreateAssayObject(counts = gene.activities)
    dataset <- NormalizeData(
      object = dataset,
      assay = 'RNA',
      normalization.method = 'LogNormalize',
      scale.factor = median(dataset$nCount_RNA)
    )
    datasets_annotated = c(datasets_annotated,dataset)  
  }
  return(datasets_annotated)
}

counts_f1 <- "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run0038/10x_outputs/HFD-8W-1_filtered_peak_bc_matrix.h5"
frag_f1 <- "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run0038/bam/HFD-8W-1_fragments.tsv.gz"
dataset1 = Read_Dataset_With_Counts(counts_f1,frag_f1)

fragments_paths = c(
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run037/bam/RC-8W-1_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run037/bam/RC-8W-2_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run037/bam/RC-8W-3_fragements.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run037/bam/RC-14W-1_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run0038/bam/RC-14W-2_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run0038/bam/RC-30W-1_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run39/bam/RC-30W-2_fragements.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run0038/bam/HFD-8W-1_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run0038/bam/HFD-8W-2_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run39/bam/HFD-24W-1_fragments.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run39/bam/HFD-24W-2_fragements.tsv.gz",
  "/home/krushna/Documents/adipose-tissue-analysis-atac-seq/dataset/Sc-ATAC seq data/run39/bam/HFD-24W-3_fragements.tsv.gz"
)


datasets = Read_Dataset_with_granges(fragments_paths, granges(dataset1))



datasets_annotated = Annotate_Datasets(datasets)

merged_data = merge(datasets_annotated[[1]], datasets_annotated[2:length(datasets_annotated)])

merged_data <- FindTopFeatures(merged_data, min.cutoff = 10)
merged_data <- RunTFIDF(merged_data)
merged_data <- RunSVD(merged_data)
merged_data <- RunUMAP(merged_data, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(merged_data, group.by = "dataset")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = datasets_annotated,
  anchor.features = rownames(datasets_annotated[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merged_data[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  k.weight = 40
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, resolution=0.1, verbose = FALSE, algorithm = 3)
DimPlot(object = integrated, label = TRUE) + NoLegend()
p2 <- DimPlot(integrated, group.by = "dataset")

print(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
integrated_data = integrated
