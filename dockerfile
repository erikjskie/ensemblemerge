FROM rocker/r-ver:4.0.4

ARG WHEN

# required
MAINTAINER Erik Skie <skiex003@umn.edu>

# copy the repo contents into the docker image at `/ensemblemerge`
COPY . /ensemblemerge

# install the dependencies of the R package located at `/portalDS`
RUN apt-get -y update -qq \ 
  && apt-get install -y --no-install-recommends \
    libgsl0-dev \
  && R -e "install.packages('devtools') \
	 install.packages('BiocManager', repos='http://cran.rstudio.com/') \
	 BiocManager::install(c('SummarizedExperiment', 'SingleCellExperiment', 'LoomExperiment')) \
	 remotes::install_github('immunogenomics/harmony') \
	 remotes::install_version('Seurat', version = '4.0.1') \
	 remotes::install_github('satijalab/seurat-wrappers') \
	 remotes::install_github('cellgeni/sceasy') \
	 remotes::install_github('erikjskie/ensemblemerge', force = FALSE, ref = 'development')"
