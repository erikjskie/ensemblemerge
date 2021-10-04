<img src = "ensemblemerge_graphics_v2.png" width = "250" height = "200" align = "right" />

<!-- badges: start -->
[![Build Status](https://travis-ci.com/erikjskie/ensemblemerge.svg?token=TzArZ5EDidamcqdAtCie&branch=main)](https://travis-ci.com/github/erikjskie/ensemblemerge)
[![Project Status: Work In Progress – The project is in a usable state and is being actively developed.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test coverage](https://codecov.io/gh/erikjskie/ensemblmerge/branch/main/graph/badge.svg)](https://codecov.io/gh/erikjskie/ensemblmerge?branch=main)
<!-- badges: end -->

# ensemblemerge
ensemblemerge is a package that implements a common work flow for several single cell RNA sequence integration methods including Seurat V3, Scanorama, Harmony, Liger, fastMNN, and bbknn. The merging process is designed to be as streamlined as possible for the user, only requiring a **SummarizedExperiment**, **SingleCellExperiment**, or **Seurat** object and a specification of which method should be used to integrate in `setParams()`.

**Instructions, documentation, and examples can be found online at: [EnsembleMerge](https://erikjskie.github.io/packages/ensemblemerge/)**

**Docker images of pre-built package environments are found online at: [DockerHub](https://hub.docker.com/repository/docker/skiex003/ensemblemerge)**

## Prerequisites
The following R packages are required for installation of ensemblemerge:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "LoomExperiment"))
devtools::install_github('satijalab/seurat-wrappers')
devtools::install_github("cellgeni/sceasy")
```

The following packages are necessary only if using the following Merge methods:
  * For Harmony:
  ```r
  devtools::install_github("immunogenomics/harmony")
  ```
  * For Scanorama:
  ```r
  system("pip3 install scanpy")
  system("pip3 install scanorama")
  ```
  * For bbknn:
  For R
  ```r
  remotes::install_github("rstudio/reticulate") #increases available memory
  ```
  Python Packages
  ```r
  system("pip3 install pynndescent") #optimizes dimension reduction
  system("pip3 install leidenalg") #optional clustering algorithm that improves bbknn performance
  ```
## Quick start guide

Once ensemblemerge is currently under reveiw for cran submission, currently ensemblemerge can be installed by:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("devtools")
devtools::install_github("erikjskie/ensemblemerge", force = FALSE)
```

Once ensemblemerge is installed, merging batches can be performed by the following:

```r
library(ensemblemerge)

data #some SingleCellExperiment or SummarizedExperiment object with batches to be integrated

params = setParams() #different integration methods can be selected by setting method = c("Seurat", "Harmony", "Liger", "Scanorama", "BBKNN")

merged_data = Merge(params, data)
```
