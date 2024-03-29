#' Run Harmony functions
#'
#' @import Seurat
#'
#' @param x SingleCellExperiment object containing single cell counts matrix
#' @param normData boolean check to normalize data
#' @param Datascaling boolean check to scale data
#' @param regressUMI boolean check to regress data
#' @param norm_method character vector specifying which normalization method to use
#' @param scale_factor numeric indicating what value to scale to
#' @param nfeatures numeric indicating how many genes to select as highly variable
#' @return returns a SingleCellExperiment object of the integrated data
#' @export
harmony_preprocess <- function(x,
                              normData = TRUE, Datascaling = TRUE, regressUMI = FALSE, 
                              norm_method = "LogNormalize", scale_factor = 10000, 
                              nfeatures = 300)
{
  b_seurat = x

  if(class(b_seurat) != "Seurat"){
    b_seurat <- Seurat::as.Seurat(b_seurat, counts = "counts", data = NULL)
  }

  if(normData){
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = norm_method, scale.factor = scale_factor)
  } else{
    b_seurat@data = b_seurat@raw.data
  }

  b_seurat <- Seurat::FindVariableFeatures(object = b_seurat, nfeatures = nfeatures)

  if(regressUMI && Datascaling) {
    b_seurat <- ScaleData(object = b_seurat, vars.to.regress = c("nUMI"))  # in case of read count data
  } else if (Datascaling) { # default option
    b_seurat <- ScaleData(object = b_seurat)
  }

  return(b_seurat)
}

#' @export
call_harmony_2 <- function(b_seurat, batch_label, npcs = 20, seed = 1, pca_name = "PCA")
{
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' needed for this function to work. Please install it.",
      call. = FALSE)
  }

  #Harmony settings
  theta_harmony = 2
  numclust = 50
  max_iter_cluster = 100


  b_seurat <- RunPCA(object = b_seurat, npcs = npcs, pc.genes = b_seurat@var.genes, verbose = FALSE, reduction.name = pca_name)

  b_seurat <- harmony::RunHarmony(object = b_seurat, batch_label, theta = theta_harmony, plot_convergence = FALSE, 
                          nclust = numclust, max.iter.cluster = max_iter_cluster)


  return(b_seurat)
}

#' @export
run_Harmony <- function(params, data){
  b_seurat = harmony_preprocess(x = data,
                              normData = params@norm_data, Datascaling = params@scaling, regressUMI = params@regressUMI, 
                              norm_method = params@norm_method, scale_factor = params@scale_factor, nfeatures = params@numHVG)
  integrated = call_harmony_2(b_seurat, batch_label = params@batch, npcs = params@npcs, seed = params@seed, pca_name = params@dimreduc_names[["PCA"]])

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
