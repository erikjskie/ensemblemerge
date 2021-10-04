#' Run Seurat V3 functions
#'
#' @import reticulate
#' @import Seurat
#'
#' @param data SingleCellExperiment object containing single cell counts matrix
#' @param params ScanoramaParams object generated from setParams(method = "BBKNN") function
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_Scanorama <- function(params, data){
  ### load up python environment ###
  reticulate::py_config()
  sc <- reticulate::import("scanpy", delay_load = TRUE)
  sr <- reticulate::import("scanorama", delay_load = TRUE)
  ad <- reticulate::import("anndata", delay_load = TRUE, convert = FALSE)
  ### convert params to python ###
  py$min_genes = params@min_genes
  py$min_cells = params@min_cells
  py$batch_size = params@batch_size
  py$return_dense = params@return_dense
  py$knn = params@knn
  py$svd_solver = params@svd_solver
  py$npcs = params@npcs
  py$batch = params@batch
  py$nhvg = params@nhvg

  ### run BBKNN integration ###
  if(class(data) == "Seurat"){
    data <- Seurat::as.SingleCellExperiment(data, counts = "counts", data = NULL)
  }
  data = suppressWarnings(sceasy::convertFormat(data, from = "sce", to = "anndata"))
  py$adata = data
  py_run_string("import numpy as np
import scanorama
import pandas as pd
import scanpy as sc

#sc.pp.filter_cells(adata, min_genes=min_genes)
#sc.pp.filter_genes(adata, min_cells=min_cells)


groups = adata.obs.groupby(batch).indices

data = []
for group in groups:
  data.append(adata[groups[group]])

for adata in data:
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=int(nhvg), inplace=True)

for adata in data:
  adata.X = adata.X.tocsr()

adatas_cor = scanorama.correct_scanpy(data, return_dimred=True)

Data = adatas_cor[0].concatenate(adatas_cor[1:len(data)], index_unique=None)

sc.tl.pca(Data, svd_solver=svd_solver, n_comps=int(npcs))

Data.obsm['X_pca'] *= -1

Data.obsm['X_scanorama'] = Data.obsm['X_pca']

Data.write(filename = 'temp.h5ad')")
  #source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  file.remove("temp.h5ad")

  if(params@return == "Seurat"){
    return(integrated)
  }
  else if(params@return == "SingleCellExperiment"){
    integrated = Seurat::as.SingleCellExperiment(integrated)
    return(integrated)
  }
  else{
    stop("Invalid return type, check params@return")
  }
}
