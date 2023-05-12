library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

# Import the data
cb_data <- readRDS("results/cb_data.rds")

#Plot by gene expression
umap_gene_plot <- FeaturePlot(cb_data, features = c(
    "Hes5",
    "Pax6",
    "Lhx1",
    "Lhx5",
    "Id1",
    "Msx3",
    "Lmx1a",
    "Msx1",
    "Msx2"
))
ggsave(
    "umap_marker_gene_plot.png",
    plot = umap_gene_plot,
    device = "png",
    path = "results",
    units = "in",
    width = 16,
    height = 16,
    bg = "white"
)
