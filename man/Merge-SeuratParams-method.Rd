% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Merge.R
\name{Merge,SeuratParams-method}
\alias{Merge,SeuratParams-method}
\title{Merge SingleCellExperiment object}
\usage{
\S4method{Merge}{SeuratParams}(params, data)
}
\arguments{
\item{params}{SeuratParams object generated by setParams()}

\item{data}{SingleCellExperiment object containing single cell counts matrix}
}
\value{
returns a SingleCellExperiment object of the integrated data
}
\description{
Merge SingleCellExperiment object
}
