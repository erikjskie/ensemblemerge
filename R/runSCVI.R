#' Run scVI functions
#'
#' @import Seurat
#' @import SingleCellExperiment
#' @import reticulate
#'
#' @param x SummarizedExpirement, SingleCellExperiment, or Seurat object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_scVI <- function(params, data){
  reticulate::py_config()
  scvi <- reticulate::import('scvi', delay_load = TRUE)
  anndata <- reticulate::import("anndata", delay_load = TRUE, convert = FALSE)
  bbknn <- reticulate::import("bbknn", delay_load = TRUE)
  sc <- reticulate::import('scanpy', delay_load = TRUE)
  np <- reticulate::import('numpy', delay_load = TRUE)
  #sce <- as.SingleCellExperiment(data)
  sce <- as.Seurat(data, data = NULL)
  sce = NormalizeData(sce)
  sce = FindVariableFeatures(sce)
  variablefeatures = VariableFeatures(sce)
  #sce <- sce[VariableFeatures(sce), ]
  sce <- data
  sce <- sce[variablefeatures, ]
  X <- t(assay(sce, 'counts'))
  adata <- anndata$AnnData(X = X)
  col_data <- colData(sce)
  adata$obs <- do.call('data.frame', c(             
    as.list(col_data),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  ))
  adata$layers <- list(counts = X)
  adata$obs$batch <- as.numeric(factor(colData(sce)[,params@batch]))
  sc$pp$normalize_total(adata, target_sum = 1e4)
  message("normalizing data")
  sc$pp$log1p(adata)
  adata$raw <- adata
  scvi$settings$seed <- 123L
  scvi$data$setup_anndata(adata, layer = "counts")
  message("Running scVI Model")
  vae <- scvi$model$SCVI(adata)
  vae$train()
  z <- vae$get_latent_representation()
  message("scVI model finished")
  data = Seurat::as.Seurat(data, data = NULL)
  rownames(z) <- colnames(data)
  data[['scvi']] <- Seurat::CreateDimReducObject(embeddings = z, key = "scVI_",  assay = DefaultAssay(data))
  if(params@return == "Seurat"){
    return(data)
  }
  else if(params@return == "SingleCellExperiment"){
    integrated = Seurat::as.SingleCellExperiment(data)
    S4Vectors::metadata(integrated)$scVI = z
    return(integrated)
  }
  else{
    stop("Invalid return type, check params@return")
  }
}