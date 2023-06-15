# LIBRARY imports
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)

save_plot <- function(name, output_path, plot_obj){
  ggsave(
    name,
    plot = plot_obj,
    device = "png",
    path = output_path,
    units = "in",
    width = 4,
    height = 4,
    bg = "white",
  )
}

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

# SECTION 1: DATA IMPORT
for (name in dataset_names) {
  dataset_path <- paste(wd, "data", name, "raw_feature_bc_matrix", sep="/")
  data <- Read10X(data.dir = dataset_path)
  seurat_obj <- CreateSeuratObject(
    counts = data, # raw count data
    project = name, # project name of Seurat object
    min.cells = 3, # filter: only features with at least 3 cells
    min.features = 200, # filter: only cells with at least 200 features
  )
  assign(name, seurat_obj)
}
#SET: all Seurat objects to merge. Object name = dataset name
merged_seurat <- merge(
  x = E10A,
  y = E10B,
  add.cell.id = dataset_names
)

head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# SECTION 2: QUALITY CONTROL
# Log(10) genes per UMI
merged_seurat$log10_genes_per_UMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# Mitochondrial ratio
merged_seurat$mito_ratio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
merged_seurat$mito_ratio <- merged_seurat@meta.data$mito_ratio / 100
# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Create sample column 
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^E10A_"))] <- "E10A"
metadata$sample[which(str_detect(metadata$cells, "^E10B_"))] <- "E10B"
# Update Seurat object with metadata
merged_seurat@meta.data <- metadata

# SECTION 3: QC Plots
# 1) Cell count plot
cell_count_plot <- ggplot(metadata, aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")
cell_count_plot
save_plot(paste(job_name, "cell_counts.png", sep = "_"), results, cell_count_plot)

# 2) UMI counts (transcripts) per cell
umi_per_cell_plot <- ggplot(metadata, aes(color=sample, x=nCount_RNA, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() + 
  ylab("Cell density") +
  geom_vline(xintercept = 500)
save_plot(paste(job_name, "UMI_per_cell.png", sep = "_"), results, umi_per_cell_plot)

# 3) Genes detected per cell
genes_per_cell_plot <- ggplot(metadata, aes(color=sample, x=nFeature_RNA, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
save_plot(paste(job_name, "genes_per_cell.png", sep = "_"), results, genes_per_cell_plot)

# 4) Transcriptome complexity
log10_genes_per_UMI_plot <- ggplot(metadata, aes(x=log10_genes_per_UMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
save_plot(paste(job_name, "transcriptome_complexity.png", sep = "_"), results, log10_genes_per_UMI_plot)

# 5) Mitochondrial counts ratio
mito_ratio_plot <- ggplot(metadata, aes(x=mito_ratio, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2)
save_plot(paste(job_name, "mito_ratio.png", sep = "_"), results, mito_ratio_plot)

#6 Scatterplot - mito-ratio
scatter_summary <- ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA, color=mito_ratio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample)
save_plot(paste(job_name, "qc_summary.png", sep = "_"), results, scatter_summary)

# SECTION 4: FILTERING
# NOTE: view the plots in previous code before continuing here. 
# Apply filtering based on parameters
filtered_seurat <- subset(
  x = merged_seurat,
  subset = (nCount_RNA >= 1200) &
    (nCount_RNA <= 15000) &
    (nFeature_RNA >= 500) &
    (log10_genes_per_UMI > 0.80) &
    (mito_ratio < 0.05)
)

#Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# SECTION 5: QC GRAPHS FOR FILTERED DATASET
#Scatterplot - mito-ratio
filtered_scatter_summary <- ggplot(filtered_seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=mito_ratio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample)
save_plot(paste(job_name, "filtered_qc_summary.png", sep = "_"), results, filtered_scatter_summary)

# Cell count plot
filtered_cell_count_plot <- ggplot(filtered_seurat@meta.data, aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells")
cell_count_plot
save_plot(paste(job_name, "filtered_cell_counts.png", sep = "_"), results, filtered_cell_count_plot)

# SECTION 6: OUTPUT
file_name = paste0(job_name, "_qc.rds")
saveRDS(filtered_seurat, file = paste(results, file_name, sep = "/"))