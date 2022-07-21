# *** MUST RUN THESE LINES *** #
# conda activate cell2location
# export PYTHONPATH="/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/conda_envs/cell2loc_e/python3.9/site-packages/"
# **************************** #

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

results_folder = "/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/"
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# conda activate SeuratDisk
# bb@reductions$pca = NULL
# bb2 = Seurat::DietSeurat(bb, dimreducs = "umap")
# SaveH5Seurat(bb2, filename = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5seurat")
# Convert("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5seurat", dest = "h5ad")

adata_vis = sc.read(f'/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_070822.h5ad')
adata_ref = sc.read(f'/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5ad')

"""
Getting the bb Anndata in the format cell2location expects.
When exporting from Seurat, the raw counts matrix was put in the adata.raw location.
While cell2location (and other tools) expects the raw counts matrix to be in adata.X.
"""
adata_ref = adata_ref.raw.to_adata()
adata_ref.X = adata_ref.X.toarray()
cell2location.models.Cell2location.setup_anndata(adata_ref)


"""
Train the model
"""
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, batch_key='sample', labels_key='seuratclusters53')
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, use_gpu=False)

adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False})
mod.save(f"{ref_run_name}", overwrite=True)  # Save model

adata_file = f"{ref_run_name}/bb_with_trained_model.h5ad"
adata_ref.write(adata_file)

from matplotlib import pyplot as plt
mod.plot_history(20)
# plt.show()
plt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/mod_train.png")
plt.clf()

"""
Load the model:
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
"""
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")


mod_8_20 = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver, N_cells_per_location=8, detection_alpha=20)
mod_8_20 = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver, N_cells_per_location=8, detection_alpha=200)
mod.view_anndata_setup()