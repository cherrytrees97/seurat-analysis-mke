# seurat-analysis-mke
 Scripts to learn how to use seurat

## Preparing environment
1. Install renv into the base R environment
2. Run `renv::init()` to generate the R virtual environment 
3. Install Seurat using `renv::install("Seurat")`
4. Run `renv::use_python()` and select the appropriate Python install you want to setup in the virual environment. 
5. Install UMAP-learn using `reticulate::py_install("umap-learn")`
6. Install limma using `renv::install("limma")` this isnt correct... need to fix. 