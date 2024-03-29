#' Preprocess a Seurat objects by the Seurat pipeline
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
#'
#' @param data a Seurat object
#' @param params a SeuratPreprocess object 
#' @param merge a BaseMerge object
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a list of Seurat objects

#' @importFrom Seurat SplitObject NormalizeData FindVariableFeatures ScaleData SelectIntegrationFeatures
#'
setMethod(
	'Preprocess',
	signature(
		data  = 'Seurat',
		params = 'SeuratPreprocess',
		merge = 'BaseMerge'
	),
	function(
		data,
		params,
		merge,
		...
	){

		stopifnot(valid(data, params))

		rc <- rowSums(GetAssayData(data, 'data') > 0)
		cc <- colSums(GetAssayData(data, 'data') > 0)

		invalid <- rc < params@min_cells
		if (any(invalid)){
			data <- data[!invalid, ]
			sprintf('Preprocess | removing %d genes that are expressed in <%d (min_cells) cells', sum(invalid), params@min_cells) %>% message()
		}

		invalid <- cc < params@min_genes
		if (any(invalid)){
			data <- data[, !invalid]
			sprintf('Preprocess | removing %d cells that have <%d (min_genes) detected genes', sum(invalid), params@min_genes) %>% message()
		}

		cn <- colnames(data)

		if (params@batchwise){
			batch_list <- SplitObject(data, split.by = params@batch)
		}else{
			batch_list <- list(data)
		}

	  for(i in 1:length(batch_list)) {
			if(params@norm_data){
				batch_list[[i]] <- NormalizeData(
					batch_list[[i]],
					normalization.method = params@norm_method,
					scale.factor = params@scale_factor,
					verbose = FALSE
				)
	    }
	    batch_list[[i]] <- FindVariableFeatures(
				batch_list[[i]],
				selection.method = params@selection.method,
				nfeatures = params@numHVG,
				verbose = FALSE
			)
		}

		if (is(merge, 'SeuratMerge')){
			new('SeuratList', batch_list)
		}else{
			# select features that are repeatedly variable across datasets for integration
			features <- SelectIntegrationFeatures(
				object.list = batch_list,
				nfeatures = params@numHVG
			)
			data@assays[[data@active.assay]]@var.features <- features
			data <- ScaleData(object = data, verbose = FALSE)
			data
		}

	}
)


#' Preprocess a Seurat object by the Scanpy pipeline
#'
#' @param data  a Seurat object
#' @param params a ScanpyPreprocess object
#' @param merge a BaseMerge object
#' @param ... Additional arguments
#' @return returns a Seurat object
#'
#setMethod(
#	'Preprocess',
#	signature(
#		data  = 'Seurat',
#		params = 'ScanpyPreprocess',
#		merge = 'BaseMerge'
#	),
#	function(
#		data,
#		params,
#		merge,
#		...
#	){
#
#		data <- as.SingleCellExperiment(data)
#		Preprocess(params, data, ...)
#	}
#)

#' Preprocess a SingleCellExperiment object by the Scanpy pipeline
#'
#' @param data a SingleCellExperiment object
#' @param params a ScanpyPreprocess object
#' @return returns a SingleCellExperiment object
#' @param ... Additional arguments
#' @importFrom reticulate import
#'
#setMethod(
#	'Preprocess',
#	signature(
#		data  = 'SingleCellExperiment',
#		params = 'ScanpyPreprocess'
#	),
#	function(
#		data,
#		params,
#		...
#	){
#
#		h5ad_file <- tempfile(fileext = '.h5ad')
#		data = sceasy::convertFormat(data, from = "sce", to = "anndata", out = h5ad_file)
#		sc <- import("scanpy")
#		sc$pp$filter_cells(data, min_genes = params@min_genes)
#		sc$pp$filter_genes(data, min_cells = params@min_cells)
#		sc$pp$normalize_total(data, target_sum = params@scale_factor)
#		sc$pp$log1p(data)
#		sc$pp$highly_variable_genes(data, min_mean = params@min_mean, max_mean= params@max_mean, min_disp = params@min_disp)
#		sc$pp$scale(data)
#		sc$tl$pca(data)
#		data$write(filename = h5ad_file)
#		data <- zellkonverter::readH5AD(h5ad_file, reader = 'R')
#		rowData(data)$highly_variable <- as.logical(rowData(data)$highly_variable)
#		unlink(h5ad_file)
#		data
#	}
#)
