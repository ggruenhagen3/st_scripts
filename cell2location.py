# *** MUST RUN THESE LINES *** #
# conda activate cell2loc_env
# export PYTHONPATH="/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/conda_envs/cell2loc_env/python3.9/site-packages/"
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
# DefaultAssay(all_merge) = "Spatial"
# st = Seurat::DietSeurat(all_merge, dimreducs = "umap", assays = "Spatial")
# SaveH5Seurat(bb2, filename = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5seurat")
# Convert("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5seurat", dest = "h5ad")
# SaveH5Seurat(st, filename = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/st_diet_070822.h5seurat")
# Convert("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/st_diet_070822.h5seurat", dest = "h5ad")

adata_vis = sc.read(f'/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/st_diet_070822.h5ad')
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
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, batch_key='sample', labels_key='seuratclusters15')
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, use_gpu=False)

adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False})
mod.save(f"{ref_run_name}", overwrite=True)  # Save model

adata_file = f"{ref_run_name}/bb_with_trained_model15.h5ad"
adata_ref.write(adata_file)

# adata_ref = sc.read_h5ad(adata_file)
# mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

from matplotlib import pyplot as plt
mod.plot_history(20)
# plt.show()
plt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/mod_train.png")
plt.clf()

mod.plot_QC()
plt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/mod_train_qc.png")
plt.clf()

inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                      for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]
inf_aver.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/bb_reference_signatures15.csv")
# inf_aver = pd.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/bb_reference_signatures.csv", index_col=0)

"""
Use the model on the Spatial data
"""
# prepare anndata for cell2location model
adata_vis = adata_vis.raw.to_adata()
adata_vis.X = adata_vis.X.toarray()
adata_vis.var_names = inf_aver.index
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
del adata_vis.var

# mod_8_20  = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver, N_cells_per_location=8, detection_alpha=20)
mod_8_200 = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver, N_cells_per_location=8, detection_alpha=200) # I think this the one we want
mod_8_200.view_anndata_setup()
mod_8_200.train(max_epochs=30000, batch_size=None, train_size=1, use_gpu=False)
adata_vis = mod_8_200.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod_8_200.adata.n_obs, 'use_gpu': False})

adata_vis.obsm['q05_cell_abundance_w_sf'].to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/cell2location_spatial15_output_q05.csv")
adata_vis.obsm['q95_cell_abundance_w_sf'].to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/cell2location_spatial15_output_q95.csv")
adata_vis.obsm['means_cell_abundance_w_sf'].to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/cell2location_spatial15_output_means.csv")
adata_vis.obsm['stds_cell_abundance_w_sf'].to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/cell2location_spatial15_output_stds.csv")

mod_8_200.save(f"{run_name}", overwrite=True)
adata_file = f"{run_name}/sp_trained15.h5ad"
adata_vis.write(adata_file)

# model = torch.load("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/cell2location_map/model.pt", map_location="cpu")
# attr_dict = model["attr_dict"]
# attr_dict["history_"]  # elbo_train