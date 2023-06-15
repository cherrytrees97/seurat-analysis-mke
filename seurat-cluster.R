# Single-cell RNA-seq - clustering - no integration because same time point.

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# SECTION 0: SET THESE PARAMETERS PRIOR TO STARTING
# SET: working directory
wd <- "~/data_analysis/cerebellum-scRNA-seq-analysis"
setwd(wd)
# SET: names of datasets to process
dataset_names <- c("E10A","E10B")
# SET: job name
job_name <- "E10_merged"
# Create results directory
results <- paste(wd, "results", job_name, sep="/")
dir.create(results)

rds_path <- paste0(results, "/", job_name, "_qc.rds")
merged_seurat_obj <- readRDS(rds_path)

#Load the cell cycle data, add it to Seurat obj metadata
load(paste0(wd, "/data/mus_musculus_g2m_s_phase_genes.RData"))
merged_seurat_obj <- CellCycleScoring(merged_seurat_obj, g2m.features = g2m_genes, s.features = s_genes)

#Normalize
merged_seurat_obj <- SCTransform(merged_seurat_obj, vst.flavor = "v2")
DefaultAssay(object = merged_seurat_obj) <- "SCT"

# Run PCA
merged_seurat_obj <- RunPCA(object = merged_seurat_obj)

# Plot PCA
PCA_by_sample <- PCAPlot(merged_seurat_obj,
        split.by = "sample")
ggsave(
  paste0(job_name, "_PCA_plot_by_sample.png"),
  plot = PCA_by_sample,
  device = "png",
  path = results,
  units = "in",
  width = 16,
  height = 16,
  bg = "white",
)


# Run UMAP
merged_seurat_obj <- RunUMAP(merged_seurat_obj, 
                             dims = 1:40,
                             reduction = "pca")

# Determine the K-nearest neighbor graph
merged_seurat_obj <- FindNeighbors(object = merged_seurat_obj, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
merged_seurat_obj <- FindClusters(object = merged_seurat_obj,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Plot the UMAP
umap_plot <- DimPlot(merged_seurat_obj,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
ggsave(
  paste0(job_name, "_umap_plot.png"),
  plot = umap_plot,
  device = "png",
  path = results,
  units = "in",
  width = 16,
  height = 16,
  bg = "white",
)

#For Joanna
markers_of_interest = c(
  "Atoh1",
  "Ptf1a",
  "Msx1",
  "Msx2",
  "Msx3",
  "Wls",
  "Lmx1a"
)

# Markers
markers = c(
  "Hes5",
  "Lhx1",
  "Lhx5",
  "Pax6",
  "Id1",
  "Msx1",
  "Msx2",
  "Msx3",
  "Ptf1a",
  "Atoh1",
  "Lmx1a",
  "Meis2",
  "Mki67"
)
msx_genes = c(
  "Msx1",
  "Msx2",
  "Msx3"
)
early_progenitor_markers = c(
  "Hes5",
  "Id1",
  "Msx3",
  "Lhx5",
  "Pax6",
  "Lmx1a",
  "Msx1",
  "Gdf7",
  "Lhx1"
)

umap_gene_plot <- FeaturePlot(merged_seurat_obj, features = markers)

ggsave(
  paste0(job_name, "_umap_gene_plot.png"),
  plot = umap_gene_plot,
  device = "png",
  path = results,
  units = "in",
  width = 16,
  height = 16,
  bg = "white",
)

#Refer to https://github.com/satijalab/seurat/issues/1623 for why we need this slot.
gene_vln_plot <- VlnPlot(merged_seurat_obj, features = markers)
ggsave(
  paste0(job_name, "_vln_plot.png"),
  plot = gene_vln_plot,
  device = "png",
  path = results,
  units = "in",
  width = 16,
  height = 16,
  bg = "white",
)

markers <- FindAllMarkers(object = merged_seurat_obj, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  write.csv(
    file = paste0(results, "/", job_name, "_cluster_genes.csv")
  )

background_genes <- merged_seurat_obj@assays$SCT@data@Dimnames[1]
background_genes %>% write.csv(file = paste0(results, "/", job_name, "_background_genes.csv"))

saveRDS(merged_seurat_obj, file = paste0(results, "/", job_name, "_clustered.rds"))