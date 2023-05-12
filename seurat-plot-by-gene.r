library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

ctrl_cb <- readRDS('results/ctrl_cb.rds')

#Plot by gene expression
umap_gene_plot <- FeaturePlot(ctrl_cb, features = c("Sox2", "Msx1", "Msx2", "Msx3", "Atoh1", "Ptf1a"))
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
