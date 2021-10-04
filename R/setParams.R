#' set parameters for integration
#'
#' @importFrom methods new
#'
#' @param method character specifying the integration method i.e "Seurat"
#' @return returns a params S4 class depending on input
#' @export
setParams <- function(method = "Seurat", ...){
    availableMethods = c("Seurat", "Scanorama", "Harmony", 
                        "Liger", "BBKNN", "SMILE", "Uncorrected", "fastMNN", "scVI") # list available availableMethods

    ### checking valid parameters ###
    if(!(method %in% availableMethods)){
      stop(sprintf("method must be the following: %s", paste(availableMethods, collapse = ", ")))
    }

    ### running integration ###
    x = switch(method,
              "Seurat" = new("SeuratParams", ...),
              "Scanorama" = new("ScanoramaParams", ...),
              "Harmony" = new("HarmonyParams", ...),
              "Liger" = new("LigerParams", ...),
              "BBKNN" = new("BBKNNParams", ...),
              "SMILE" = new("SMILEParams", ...),
              "Uncorrected" = new("UncorrectedParams", ...),
              "fastMNN" = new("FastMNNParams", ...),
              "scVI" = new("scVIParams", ...))
    return(x)
}