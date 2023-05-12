library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

cb_data <- readRDS("results/cb_data.rds")

#Plot by gene expression
umap_gene_plot <- FeaturePlot(cb_data, features = c(
    "Wls",
    "Ascl1",
    "Neurog1",
    "Neurog2",
    "Atoh1",
    "Ptf1a",
    "Sox2",
    "Msx1",
    "Msx2",
    "Msx3"
))
ggsave(
    "umap_gene_plot.png",
    plot = umap_gene_plot,
    device = "png",
    path = "results",
    units = "in",
    width = 16,
    height = 16,
    bg = "white"
)
