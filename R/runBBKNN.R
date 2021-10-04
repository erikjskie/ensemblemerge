#' Run Seurat V3 functions
#'
#' @import reticulate
#' @import Seurat
#' @import methods
#' @import SparseM
#' @import SingleCellExperiment
#'
#' @param data SingleCellExperiment object containing single cell counts matrix
#' @param params BBKNNParams object generated from setParams(method = "BBKNN") function
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_BBKNN <- function(params, data){
  ### load up python environment ###
  reticulate::py_config()
  sc <- reticulate::import("scanpy", delay_load = TRUE)
  ad <- reticulate::import("anndata", delay_load = TRUE, convert = FALSE)
  bbknn <- reticulate::import("bbknn", delay_load = TRUE)
  ### convert params to python ###
  py$min_genes = params@min_genes
  py$min_cells = params@min_cells
  py$svd_solver = params@svd_solver
  py$scale_factor = params@scale_factor
  py$npcs = params@npcs
  py$nhvg = params@nhvg
  py$batch = params@batch
  py$norm_data = params@norm_data
  py$scaling = params@scaling
  py$regressUMI = params@regressUMI
  py$norm_method = params@norm_method
  py$save_knn = params@save_knn
  py$copy = params@copy
  py$neighbors_within_batch = params@neighbors_within_batch
  py$approx = params@approx
  py$trim = params@trim
  py$n_neighbors = params@n_neighbors
  py$confounder_key = params@confounder_key
  py$ridge_regress = params@ridge_regress

  ### run BBKNN integration ###
  if(class(data) == "Seurat"){
    data <- Seurat::as.SingleCellExperiment(data, counts = "counts", data = NULL)
  }
  data = suppressWarnings(sceasy::convertFormat(data, from = "sce", to = "anndata", out = "temp.h5ad"))
  #py$adata = data
  #Sys.sleep(10)
  py_run_string("import numpy as np
import bbknn
import scanpy as sc

adata = sc.read('temp.h5ad')")

py_run_string("sc.pp.filter_cells(adata, min_genes=3)
sc.pp.filter_genes(adata, min_cells=500)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata)
sc.tl.pca(adata)")

py_run_string("import numpy as np
import bbknn
import scanpy as sc

#perform bbknn integration
if ridge_regress == False :
  bbknn.bbknn(adata, batch_key=batch)
else:
  if confounder_key == 'leiden':
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.4)
    bbknn.ridge_regression(adata, batch_key=[batch], confounder_key=['leiden'])
    sc.pp.pca(adata)
    bbknn.bbknn(adata, batch_key=batch)
  else:
    bbknn.ridge_regression(adata, batch_key=[batch], confounder_key=['CellType'])
    sc.pp.pca(adata)
    bbknn.bbknn(adata, batch_key=batch)")


py_run_string("adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
adata.write(filename = 'temp.h5ad')")

  #source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  file.remove("temp.h5ad")

  bbknn = py$adata$obsp["distances"]
  py$adata = NULL #remove adata object from python to reduce memory
  bbknn = as(bbknn, "matrix.csr")
  bbknn = as(bbknn, "dgCMatrix")
  rownames(bbknn) = colnames(integrated)
  colnames(bbknn) = colnames(integrated)
  bbknn = Seurat::as.Graph(bbknn)
  integrated[[params@graph_name]] <- bbknn
  bbknn = Seurat::as.Neighbor(bbknn)
  integrated[[params@nn_name]] <- bbknn

  if(params@return == "Seurat"){
    return(integrated)
  }
  else if(params@return == "SingleCellExperiment"){
    integrated = Seurat::as.SingleCellExperiment(integrated)
    S4Vectors::metadata(integrated)$bbknn = bbknn
    return(integrated)
  }
  else{
    stop("Invalid return type, check params@return")
  }
}
