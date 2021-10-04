#' Get valid integration methods 
#'
#'
#' @return returns vector of valid names for integration methods
#' @export
getMethods <- function(){
  Methods = c("Seurat", "Harmony", "Liger", "Scanorama", "BBKNN", "fastMNN", "Uncorrected", "scVI")
  return(Methods)
}

#' Get valid parameter classes 
#'
#'
#' @return returns vector of valid names for integration methods
getParams <- function(){
  Params = c("SeuratParams", "HarmonyParams", "LigerParams", "ScanoramaParams", "BBKNNParams", "FastMNNParams", "UncorrectedParams")
  return(Params)
}


#' Run ensemblemerge functions
#'
#' @import Seurat
#' @import SingleCellExperiment
#' @import Matrix
#'
#' @param x SummarizedExpirement object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
EnsembleMerge <- function(data, methods = c("Seurat", "Harmony"), return = "Seurat", a = 1, b = 1, file = NA, normalization_method = "LogNormalize", latent = FALSE, ...){

  #* Add error if methods is not of class character or xParams
  if(class(methods) != "character" & class(methods) != "list"){
      stop("'methods' must be of class character or list of of Params objects set by setParams()")
  }

  #* Add error if methods length is not >= 2
  if(length(methods) < 2){
    stop("'methods' must be a list or vector of at least 2 valid integration methods")
  }

  if(class(methods[1][[1]]) %in% getParams()){
    #* Add error if invalid methods present
    for(i in 1:length(methods)){
      if(!(methods[i][[1]] %in% getParams())){
        stop(sprintf("'%s' is not an available method, check available methods with getParams()", methods[i]))
      }
    }
  }
  else{
    new_methods = list() #convert methods to params objects
    #* Add error if invalid methods present, also convert methods to params objects if no error
    for(i in 1:length(methods)){
      if(!(methods[i] %in% getMethods())){
        stop(sprintf("'%s' is not an available method, check available methods with getMethods()", methods[i]))
      }
      else{
        new_methods[[i]] = setParams(method = methods[i]) #assign params object to methods
      }
    }
    methods = new_methods #exchange variable name for downstream
    new_methods = NULL #remove temp methods holder
  }
  

  #* Add error if data is not a valid input type
  if(!(class(data) %in% c("Seurat", "SummarizedExperiment", "SingleCellExperiment"))){
      stop(sprintf("'%s' is not an valid input for data, available inputs are 'Seurat', 'SummarizedExperiment', 'SingleCellExperiment'", class(data)))
  }
  
  if(class(data) == 'SummarizedExperiment'){
      data = as(data, "SingleCellExperiment") # set data as SingleCellExperiment
    }

    if(class(data) == "Seurat"){
      data = Seurat::as.SingleCellExperiment(data)
    }

  ng = lapply(methods, function(ng){
    if(!is.na(normalization_method) & normalization_method == "CLR"){
      ng = getNeighborGraph(ng, data, normalization_method, latent = latent)
    }
    else{
      ng = getNeighborGraph(ng, data, normalization_method, latent)
    } 
  })

  return(ng)


  agreement = function(x){
    temp = base::Reduce(rbind, lapply(x, summary))
    temp = temp[which(duplicated(temp[,c("i", "j")])),]
    temp = unique(temp[,c("i", "j")])
    message("applying boolean values")
    ret = Matrix::Matrix(nrow = nrow(x[[1]]), ncol = ncol(x[[1]]), data = 0, sparse = TRUE)
    ret = as(ret, "dgCMatrix")
    ret = as(ret, "dgCMatrix")
    ret[as.matrix(temp[,c("i", "j")])] = 1
    return(ret)
  }

  #sigmoid = suppressWarnings(function(x){0.5 * (1 + sin((x*pi)-(pi/2)))})
  sigmoid = function(x, a = 1, b = 1){1/(1+a*(x/(1-x))^-(b*2))}

  agreement = agreement(ng)

  wt = list()
  for(i in 1:length(ng)){
    weight = ng[[i]]*agreement/sqrt((ng[[i]]*ng[[i]])%*%(agreement*agreement))
    weight[is.na(weight)] = 0
    weight = sum(sum(weight))
    message(sprintf("raw weight is %s", weight))
    wt = append(wt, weight)
  }

  tmpwt = as.vector(unlist(wt))
  for(i in 1:length(wt)){
    wt[i] = (wt[[i]] - nrow(data))/(max(tmpwt) - nrow(data)) #set 0 to min(tmpwt) to fix if needed
    message(sprintf("normalized weight is %f", wt[[i]]))
    wt[i] = sigmoid(wt[[i]], a , b)
    message(sprintf("sig-transformed weight is now %f", wt[i]))
  }

#tmpwt = as.vector(unlist(wt))
#for(i in 1:length(wt)){
  #wt[i] = (wt[[i]] - min(tmpwt))/(max(tmpwt) - min(tmpwt))
  #message(sprintf("renormalized weight is %f", wt[[i]]))
#}

kg = list()
for(i in 1:length(wt)){
  x = ng[[i]]*wt[[i]]
  x[is.na(x)] = 0
  x = as(x, "dgCMatrix")
  kg = append(kg, x)
}

  ng = Reduce("+", kg)/Reduce("+", wt)
  ng[is.na(ng)] = 0
  ng = as(ng, "dgCMatrix")

  weight = ng*agreement/sqrt((ng*ng)%*%(agreement*agreement))
  weight[is.na(weight)] = 0
  weight = sum(sum(weight))
  tmpwt = c(tmpwt,weight)
  weight = (weight - min(tmpwt))/(max(tmpwt) - min(tmpwt))
  message(sprintf("weight for ensemblemerge is %s", weight))
  wt = append(wt, weight)

  if(!is.na(file)){
    score = as.data.frame(c(method = "Method", weight = "Weight"), ncol = 2) %>% t()
    rownames(score) = NULL
    write.table(score, file = file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  }

  for(i in 1:length(methods)){
    message(sprintf("%s method scores %f out of 1", methods[[i]]@name, wt[[i]]))
    if(!is.na(file)){
    score = as.data.frame(c(method = methods[[i]]@name, weight = wt[[i]]), nrow = 2) %>% t()
    rownames(score) = NULL
    write.table(score, file = file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    }
  }
  message(sprintf("EnsembleMerge method scores %f out of 1", wt[[length(wt)]]))
  if(!is.na(file)){
    score = as.data.frame(c(method = "EnsembleMerge", weight = wt[[length(wt)]]), nrow = 2) %>% t()
    rownames(score) = NULL
    write.table(score, file = file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  }

  if(return == "SingleCellExperiment"){
    S4Vectors::metadata(data)$EnsembleMerge <- ng
  }
  else if(return == "Seurat"){
    data = as.Seurat(data, counts = "counts", data = NULL)
    data[["EnsembleMerge"]] = as.Graph(ng)
  }
  

  return(data)
}
