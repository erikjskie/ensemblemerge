#' parent

#' 
#'
#' @export
setClass(
	'parent', 
	representation(
		batch = "character",
    dimreduc_names = "character",
    return = "character",
    name = "character"
	),
  prototype(batch = "batch",
            dimreduc_names = c("PCA" = "pca",
                              "UMAP" = "umap",
                              "tSNE" = "tsne"),
            return = "SingleCellExperiment",
            name = "parent")
)

#' SeuratMerge
#'
#' @export
setClass(
	'SeuratMerge', 
	representation(
		npcs = "numeric",
    seed = "numeric",
    dims = "numeric",
    k.weight = "numeric"
	),
	contains = c('parent'),
  prototype(npcs = 20, 
            seed = 10,
            dims = 20,
            k.weight = 100)
)

#' SeuratNormalize
#'
#' @export
setClass(
	'SeuratNormalize', 
	representation(
		norm_data = "logical",
    scaling = "logical",
    regressUMI = "logical",
    min_cells = "numeric",
    min_genes = "numeric",
    norm_method = "character",
    scale_factor = "numeric"
	),
	contains = c('parent'),
  prototype(norm_data = TRUE,
            scaling = TRUE,
            regressUMI = FALSE,
            min_cells = 10,
            min_genes = 300,
            norm_method = "LogNormalize",
            scale_factor = 10000)
)

#' SeuratHVG
#'
#' @export
setClass(
	'SeuratHVG', 
	representation(
    numVG = "numeric",
    numHVG = "numeric"
	),
	contains = c('parent'),
  prototype(numVG = 300,
            numHVG = 2000)
)

#' SeuratParams
#'
#' @export
setClass(
	'SeuratParams', 
	representation(altExp_names = "character"),
	contains = c('SeuratNormalize',
                'SeuratHVG',
                'SeuratMerge'),
  prototype(altExp_names = "RNA",
  name = "Seurat")
)

#' HarmonyMerge
#'
#' @export
setClass(
	'HarmonyMerge', 
	representation(
    theta_harmony = "numeric",
    npcs = "numeric",
    seed = "numeric",
    num_clust = "numeric",
    max_iter_cluster = "numeric"
	),
	contains = c('parent'),
  prototype(theta_harmony = 2,
            num_clust = 50,
            max_iter_cluster = 100,
            npcs = 20,
            seed = 10)
)

#' HarmonyParams
#'
#' @export
setClass(
	'HarmonyParams', 
	representation(),
	contains = c('SeuratNormalize',
                'SeuratHVG',
                'HarmonyMerge'),
                prototype(name = "Harmony")
)

#' UncorrectedParams
#'
#' @export
setClass(
	'UncorrectedParams', 
	representation(vars_to_regress = "character",
                hvg = "logical"),
	contains = c('SeuratNormalize',
                'SeuratHVG',
                'SeuratMerge'),
  prototype(vars_to_regress = c("nUMI"),
                                hvg = TRUE,
                                name = "Uncorrected")
)

#' FastMNNParams
#'
#' @export
setClass(
	'FastMNNParams', 
	representation(),
	contains = c('UncorrectedParams'),
  prototype(name = "FastMNN")
)

#' LigerNormalize
#'
#' @export
setClass(
	'LigerNormalize', 
	representation(
	),
	contains = c('parent')
)

#' LigerHVG
#'
#' @export
setClass(
	'LigerHVG', 
	representation(
    var_threshold = "numeric"
	),
	contains = c('parent'),
  prototype(var_threshold = 0.1)
)

#' LigerMerge
#'
#' @export
setClass(
	'LigerMerge', 
	representation(
    k = "numeric",
    nrep = "numeric",
    lambda = "numeric"
	),
	contains = c('parent'),
  prototype(k = 20,
            nrep = 3,
            lambda = 5)
)

#' LigerParams
#'
#' @export
setClass(
	'LigerParams', 
	representation(),
	contains = c('LigerNormalize',
                'LigerHVG',
                'LigerMerge',
                'UncorrectedParams'),
  prototype(name = "Liger")
)

#' BBKNNNormalize
#'
#' @export
setClass(
	'BBKNNNormalize', 
	representation(
    min_genes = "numeric",
    min_cells = "numeric",
    svd_solver = "character",
    scale_factor = "numeric",
    npcs = "numeric",
    nhvg = "numeric",
    n_neighbors = "numeric",
    confounder_key = "character"
	),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_genes = 5,
            svd_solver = "arpack",
            scale_factor = 10000,
            npcs = 20,
            nhvg = 2000,
            n_neighbors = 10,
            confounder_key = "leiden")
)

#' BBKNNMerge
#'
#' @export
setClass(
	'BBKNNMerge', 
	representation(
    save_knn = "logical",
    copy = "logical",
    neighbors_within_batch = "numeric",
    approx = "logical",
    trim = "numeric",
    graph_name = "character",
    nn_name = "character",
    ridge_regress = "logical"
	),
	contains = c('parent'),
  prototype(save_knn = TRUE,
            copy = TRUE,
            neighbors_within_batch = 5,
            approx = FALSE,
            trim = 50,
            graph_name = "bbknn_graph",
            nn_name = "bbknn",
            ridge_regress = TRUE)
)

#' BBKNNParams
#'
#' @export
setClass(
	'BBKNNParams', 
	representation(),
	contains = c('BBKNNNormalize',
              'SeuratNormalize',
              'BBKNNMerge'),
  prototype(name = "bbknn")
)

#' ScanoramaParams
#'
#' @export
setClass(
	'ScanoramaParams', 
	representation(min_cells = "numeric",
          min_genes = "numeric",
          batch_size = "numeric",
          return_dense = "logical",
          knn = "numeric",
          svd_solver = "character",
          npcs = "numeric",
          nhvg = "numeric"),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_cells = 5,
            batch_size = 30,
            return_dense = TRUE,
            knn = 10,
            svd_solver = "arpack",
            npcs = 20,
            nhvg = 2000,
            name = "Scanorama")
)

#' SMILEParams
#'
#' @export
setClass(
	'SMILEParams', 
	representation(min_cells = "numeric",
          min_genes = "numeric",
          batch_size = "numeric",
          return_dense = "logical",
          knn = "numeric",
          svd_solver = "character",
          npcs = "numeric"),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_cells = 5,
            batch_size = 30,
            return_dense = TRUE,
            knn = 10,
            svd_solver = "arpack",
            npcs = 20,
            name = "SMILE")
)

#' scVIParams
#'
#' @export
setClass(
	'scVIParams', 
	representation(min_cells = "numeric",
          min_genes = "numeric",
          batch_size = "numeric",
          return_dense = "logical",
          knn = "numeric",
          svd_solver = "character",
          npcs = "numeric",
          nhvg = "numeric",
          use_cuda = "logical"),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_cells = 5,
            batch_size = 30,
            return_dense = TRUE,
            knn = 10,
            svd_solver = "arpack",
            npcs = 20,
            nhvg = 2000,
            name = "Scanorama",
            use_cuda = TRUE)
)