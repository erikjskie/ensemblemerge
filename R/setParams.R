#' Set parameters for integration
#'
#' @param methods The method names (default: 'Seurat')
#' @param npcs the size of latent dimension (default: 20L)
#' @param raw_assay Assay field in the Seurat object for the raw counts (default: 'RNA')
#' @param batch the name for the batch label (default: 'batch')
#' @param ... Additional arguments for new('SeuratPreprocess', ...)
#'
#' @importFrom methods new
#'
#' @return If only one method is used, the output is a parameter object for each method; If more than one methods are used, the output is a `EnsembleMergeParams` object
#'
#' @export
#'
setParams <- function(
	methods = "Seurat", 
	npcs = 20L,
	raw_assay = 'RNA',
	batch = 'batch',
	...
){

	.check_method(methods)

	ml <- new('MethodList', lapply(methods, function(method){
		m <- .new_object(method)
		m@npcs <- npcs
		m@raw_assay <- raw_assay
		m@batch <- batch
		m
	}))

	params <- new('Params')
	params@preprocess <- new('SeuratPreprocess', raw_assay = raw_assay, batch = batch, ...)
	params@constituent <- ml
	names(params@constituent) <- methods
	params

}

.new_object <- function(x, ...){
  switch(
		x,
		"Seurat" = new("SeuratMerge", ...),
		"Harmony" = new("HarmonyMerge", ...),
		"Scanorama" = new("ScanoramaMerge", ...),
		"Liger" = new("LigerMerge", ...),
		"BBKNN" = new("BBKNNMerge", ...),
		"Uncorrected" = new("UncorrectedMerge", ...),
		"fastMNN" = new("FastMNNMerge", ...),
		"scVI" = new("scVIMerge", ...)
	)
}
