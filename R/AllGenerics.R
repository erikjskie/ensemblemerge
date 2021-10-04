#' Generic merge function
#' @export
setGeneric("Merge", function(params, data) standardGeneric("Merge"))
setGeneric("Score", function(x, batch_label = NULL, method = "kBET", ...) standardGeneric("Score"))
setGeneric("getNeighborGraph", function(params, data, latent, ...) standardGeneric("getNeighborGraph"))