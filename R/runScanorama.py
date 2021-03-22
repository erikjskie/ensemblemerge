import numpy as np
import scanorama
import pandas as pd
import scanpy as sc

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)


groups = adata.obs.groupby(batch).indices

data = []
for group in groups:
  data.append(adata[groups[group]])

for adata in data:
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=int(nhvg), inplace=True)

for adata in data:
  adata.X = adata.X.tocsr()

adatas_cor = scanorama.correct_scanpy(data, return_dimred=True)

Data = adatas_cor[0].concatenate(adatas_cor[1:len(data)], index_unique=None)

sc.tl.pca(Data, svd_solver=svd_solver, n_comps=int(npcs))

Data.obsm["X_pca"] *= -1

Data.obsm["X_scanorama"] = Data.obsm["X_pca"]

Data.write(filename = "temp.h5ad")