# Running Seurat data on our own 10X data. 
# Pipeline is adapted from the ctrl_cb3K vignette. 
# Can be found at https://satijalab.org/seurat/articles/ctrl_cb3k_tutorial.html

library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

# Set working directory for this script
setwd("C:/Users/Micha/Documents/data_analysis/seurat_environment")

# SECTION 1: Load the ctrl_cb dataset
cb_data <- Read10X(data.dir = "data/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
ctrl_cb <- CreateSeuratObject(
    counts = cb_data, # raw count data
    project = "ctrl_cb", # project name of Seurat object
    min.cells = 3, # filter: only features with at least 3 cells
    min.features = 200, # filter: only cells with at least 200 features
)
ctrl_cb

# SECTION 2: Standard pre-processing workflow
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Get the percent mitochondrial gene for each cell
ctrl_cb[["percent.mt"]] <- PercentageFeatureSet(ctrl_cb, pattern = "^MT-")

# Visualize QC metrics as a violin plot
cb_vlnplot <- VlnPlot(
    ctrl_cb, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    ncol = 3,
)
#Save plot - used for rest of figures in this script.
ggsave(
    "qc_summary_plot.png",
    plot = cb_vlnplot,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
mt_scatter <- FeatureScatter(
    ctrl_cb, 
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
)
ggsave(
    "mt_scatter.png",
    plot = mt_scatter,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

cov_depth_scatter <- FeatureScatter(
    ctrl_cb, 
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
)
ggsave(
    "feature_depth_scatter.png",
    plot = cov_depth_scatter,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

#Get cells where # genes is between 1000 and 6000, and % of mitochondrial genes < 5
ctrl_cb <- subset(
    ctrl_cb,
    subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5
)

# SECTION 3: Normalizing Data
# Features are normalized based on total expression then multiplied by scale factor.
ctrl_cb <- NormalizeData(ctrl_cb, normalization.method = "LogNormalize", scale.factor = 10000)

# SECTION 4: Identification of highly variable features (feature selection)
ctrl_cb <- FindVariableFeatures(ctrl_cb, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ctrl_cb), 10)

# plot variable features with and without labels
variable_feature <- VariableFeaturePlot(ctrl_cb)
ggsave(
    "variable_genes.png",
    plot = variable_feature,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

variable_feature_labeled <- LabelPoints(plot = variable_feature, points = top10, repel = TRUE)
ggsave(
    "variable_genes_labeled.png",
    plot = variable_feature_labeled,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)


# SECTION 5: Scaling the data 
# shift expression of each gene so that mean expression across cells is 0
# scale expression of each gene so that the variance across cells is 1

all.genes <- rownames(ctrl_cb) # needed for heatmaps
ctrl_cb <- ScaleData(ctrl_cb, features = all.genes)

# SECTION 6: Perform linear dimensional reduction w/ PCA
ctrl_cb <- RunPCA(ctrl_cb, features = VariableFeatures(object = ctrl_cb))

# Examine and visualize PCA results a few different ways
# Method 1
print(ctrl_cb[["pca"]], dims = 1:5, nfeatures = 5)
# Method 2
viz_dim_load <- VizDimLoadings(ctrl_cb, dims = 1:2, reduction = "pca")
ggsave(
    "PCA_top-genes.png",
    plot = viz_dim_load,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

# Method 3
pca_plot <- DimPlot(ctrl_cb, reduction = "pca")
ggsave(
    "PCA_plot.png",
    plot = pca_plot,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

# SECTION 7: Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#ctrl_cb <- JackStraw(ctrl_cb, num.replicate = 100)
#ctrl_cb <- ScoreJackStraw(ctrl_cb, dims = 1:20)
#JackStrawPlot(ctrl_cb, dims = 1:15)

cb_elbow <- ElbowPlot(ctrl_cb)
ggsave(
    "elbow_plot.png",
    plot = cb_elbow,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

# SECTION 8: Cluster the cells (assuming 20 dimensions for now)
ctrl_cb <- FindNeighbors(ctrl_cb, dims = 1:20)
ctrl_cb <- FindClusters(ctrl_cb, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(ctrl_cb), 5)

# SECTION 9: Run non-linear dimensional reduction (UMAP/tSNE)
# Note that UMAP is the method designed for scRNA-seq datasets.
ctrl_cb <- RunUMAP(ctrl_cb, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
umap_plot <- DimPlot(ctrl_cb, reduction = "umap")
ggsave(
    "umap_plot.png",
    plot = umap_plot,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)

saveRDS(ctrl_cb, file = "results/ctrl_cb.rds")

# SECTION 10: Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
ctrl_cb_markers <- FindAllMarkers(ctrl_cb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# NOTE: %>% is a forward pipe operator.
ctrl_cb_markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

#Plot by gene expression
umap_gene_plot <- FeaturePlot(ctrl_cb, features = c("Sox2", "Msx1", "Msx2", "Msx3"))
ggsave(
    "umap_gene_plot.png",
    plot = umap_gene_plot,
    device = "png",
    path = "results",
    units = "in",
    width = 8,
    height = 8,
    bg = "white",
)
