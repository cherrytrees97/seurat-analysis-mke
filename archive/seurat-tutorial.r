# Seurat tutorial: vignette for guided clustering
# Can be found at https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(dplyr)
library(Seurat)
library(patchwork)

# Set working directory for this script
setwd("C:/Users/Micha/Documents/data_analysis/seurat_environment")

# SECTION 1: Load the PBMC dataset
pbmc_data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
# pbmc_data <- Read10X(data.dir = "../data/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(
    counts = pbmc_data, # raw count data
    project = "pbmc3k", # project name of Seurat object
    min.cells = 3, # filter: only features with at least 3 cells
    min.features = 200, # filter: only cells with at least 200 features
)
pbmc

# SECTION 2: Standard pre-processing workflow
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Get the percent mitochondrial gene for each cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.copy(png, "results/vlnplot.png")
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(
    pbmc, 
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
)
plot2 <- FeatureScatter(
    pbmc, 
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
)
plot1 + plot2
dev.copy(png, "results/qc_scatter.png")
dev.off()

#Get cells where # genes is between 200 and 2500, and % of mitochondrial genes < 5
pbmc <- subset(
    pbmc, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)

# SECTION 3: Normalizing Data
# Features are normalized based on total expression then multiplied by scale factor.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# SECTION 4: Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.copy(png, "results/variable_genes.png")
dev.off()

# SECTION 5: Scaling the data 
# shift expression of each gene so that mean expression across cells is 0
# scale expression of each gene so that the variance across cells is 1

all.genes <- rownames(pbmc) # needed for heatmaps
pbmc <- ScaleData(pbmc, features = all.genes)

# SECTION 6: Perform linear dimensional reduction w/ PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
# Method 1
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# Method 2
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# Method 3
DimPlot(pbmc, reduction = "pca")
# Method 4
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# Method 5
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# SECTION 7: Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

# SECTION 8: Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# SECTION 9: Run non-linear dimensional reduction (UMAP/tSNE)
# Note that UMAP is the method designed for scRNA-seq datasets.
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "results/pbmc_tutorial.rds")

# SECTION 10: Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# NOTE: %>% is a forward pipe operator.
pbmc_markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

#Plot by gene expression
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

# SECTION 11: Assigning cell type identity to clusters
#c() = combine values into a vector or list
new_cluster_ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new_cluster_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_cluster_ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "results/pbmc3k_final.rds")