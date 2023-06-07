# Running Seurat data on our own 10X data.
# Pipeline is adapted from the PBMC3K vignette.
# Can be found at https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# LIBRARY imports
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

# SECTION 0: SET THESE PARAMETERS PRIOR TO STARTING
# Set working directory for this script
wd <- "~/data_analysis/2023-05-17_E10A-Seurat-analysis"
setwd(wd)
# Set results directory and create it.
results <- paste(wd, "results", sep="/")
dir.create(results)

# SECTION 1: Load the cb_data dataset
cb_data <- Read10X(data.dir = "data/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
cb_data <- CreateSeuratObject(
  counts = cb_data, # raw count data
  project = "cb_data", # project name of Seurat object
  min.cells = 3, # filter: only features with at least 3 cells
  min.features = 200, # filter: only cells with at least 200 features
)
cb_data

# SECTION 2: Standard pre-processing workflow

# Visualize QC metrics as a violin plot
cb_data[["percent.mt"]] <- PercentageFeatureSet(cb_data, pattern = "?^mt-")

cb_vlnplot <- VlnPlot(
  cb_data,
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
cov_depth_scatter <- FeatureScatter(
  cb_data, 
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
cb_vlnplot
min_nFeature_RNA <- 1000
max_nFeature_RNA <- 5000

cb_data <- subset(
  cb_data,
  subset = (nFeature_RNA > min_nFeature_RNA) & (nFeature_RNA < max_nFeature_RNA) & (percent.mt < 5)
)

# SECTION 3: Normalizing Data
# Features are normalized based on total expression then multiplied by scale factor.
cb_data <- NormalizeData(cb_data, normalization.method = "LogNormalize", scale.factor = 10000)

# SECTION 4: Identification of highly variable features (feature selection)
cb_data <- FindVariableFeatures(cb_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(cb_data), 20)
write.csv(
  top20,
  file = 'results/top_20_genes.csv'
)


# plot variable features with and without labels
variable_feature <- VariableFeaturePlot(cb_data)
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

variable_feature_labeled <- LabelPoints(plot = variable_feature, points = top20, repel = TRUE)
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

all.genes <- rownames(cb_data) # needed for heatmaps
cb_data <- ScaleData(cb_data, features = all.genes)

# SECTION 6: Perform linear dimensional reduction w/ PCA
cb_data <- RunPCA(cb_data, features = VariableFeatures(object = cb_data))

# Examine and visualize PCA results a few different ways
# Method 1
print(cb_data[["pca"]], dims = 1:5, nfeatures = 5)

# Method 2
viz_dim_load <- VizDimLoadings(cb_data, dims = 1:2, reduction = "pca")
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
pca_plot <- DimPlot(cb_data, reduction = "pca")
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
#cb_data <- JackStraw(cb_data, num.replicate = 100)
#cb_data <- ScoreJackStraw(cb_data, dims = 1:20)
#JackStrawPlot(cb_data, dims = 1:15)

cb_elbow <- ElbowPlot(cb_data)
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
cb_data <- FindNeighbors(cb_data, dims = 1:20)
cb_data <- FindClusters(cb_data, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(cb_data), 5)

# SECTION 9: Run non-linear dimensional reduction (UMAP/tSNE)
# Note that UMAP is the method designed for scRNA-seq datasets.
cb_data <- RunUMAP(cb_data, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
umap_plot <- DimPlot(cb_data, reduction = "umap")
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

saveRDS(cb_data, file = "results/cb_data.rds")

# SECTION 10: Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
cb_data_markers <- FindAllMarkers(cb_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# NOTE: %>% is a forward pipe operator.
cb_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  write.csv(
    file = "results/cluster_genes.csv"
  )