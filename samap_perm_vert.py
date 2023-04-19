from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles,
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import multiprocessing
import itertools
from functools import reduce
from sklearn.preprocessing import LabelBinarizer

# Helper Functions
def my_mapper(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
    meta = sm.samap.adata.obs
    all_mz_cluster = meta[mz_col].unique()
    all_mm_cluster = meta[mm_col].unique()
    all_mz_cluster.sort()
    all_mm_cluster.sort()
    lb = LabelBinarizer(sparse_output=True)
    knn_mz  = lb.fit_transform(sm.samap.adata.obs[mz_col]).T.dot(knn)
    knn_sum = pd.DataFrame(lb.fit_transform(sm.samap.adata.obs[mm_col]).T.dot(knn_mz.T).todense())
    mz_counts = meta[mz_col].value_counts().sort_index()
    mm_counts = meta[mm_col].value_counts().sort_index()
    mz_mm_counts = pd.DataFrame(np.array(mz_counts) * np.array(mm_counts).reshape(len(mm_counts), 1))
    knn_mean = knn_sum / mz_mm_counts
    knn_mean.index = all_mm_cluster
    knn_mean.columns = all_mz_cluster
    knn_mean = knn_mean.drop("unassigned", axis=0)
    knn_mean = knn_mean.drop("unassigned", axis=1)
    return(knn_mean)

def perm_mapper(seed_num, mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
    np.random.seed(seed_num)
    sm.samap.adata.obs['mz_perm'] = sm.samap.adata.obs[mz_col]
    sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'])
    sm.samap.adata.obs['mm_perm'] = sm.samap.adata.obs[mm_col]
    sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'])
    perm_map = my_mapper('mz_perm', 'mm_perm')
    perm_map = perm_map.values.tolist()
    perm_map = list(itertools.chain(*perm_map))
    return(perm_map)

def my_p(x):
    return( 1 - (np.count_nonzero(perm_dist<x) / perm_dist.size))

with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/vert2_pair_sm_ran.pkl', 'rb') as out_file:
    sm = pickle.load(out_file)

turtle_meta = pd.read_csv("/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/bcs/data/turtle_meta.csv", index_col = 0)
sm.samap.adata.obs['cp_areaident'] = 'unassigned'
sm.samap.adata.obs.loc[turtle_meta.index, 'cp_areaident'] = turtle_meta['areaident']

sm.samap.adata.obs['tg_region_cluster'] = sm.samap.adata.obs['tg_region'] + "_" + sm.samap.adata.obs['tg_cluster_orig2']
sm.samap.adata.obs.loc[sm.samap.adata.obs['tg_region_cluster'] == "unassigned_unassigned", 'tg_region_cluster'] = "unassigned"

am_meta = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/axolotl_metadata.csv', index_col=0)
am_meta['cellclusters'] = am_meta['cellclusters'].replace(["glut_SUBSET_", "GABA_SUBSET_", "epen_clus_", "npc_SUBSET_", "oligodendrocyte_", "endothelial_", "microglia_"], ["GLUT", "GABA", "EPEN", "NB", "OLIG", "ENDO", "MG"], regex=True)
sm.samap.adata.obs['am_cluster'] = 'unassigned'
sm.samap.adata.obs.loc[am_meta.index,'am_cluster'] = am_meta['cellclusters']
sm.samap.adata.obs['am_region_cluster'] = 'unassigned'
sm.samap.adata.obs.loc[am_meta.index,'am_region_cluster'] = sm.samap.adata.obs.loc[am_meta.index,'am_region'] + "_" + am_meta['cellclusters']

sm.samap.leiden_clustering(res=3)
sm.samap.adata.obs.to_csv("~/scratch/bcs/samc/vert2.csv")
sm.samap.adata.obs.to_csv("~/scratch/bcs/samc/vert2_turtle.csv")
sm.samap.adata.obs.to_csv("~/scratch/bcs/samc/vert2_mouse.csv")
sm.samap.adata.obs.to_csv("~/scratch/bcs/samc/vert2_bird.csv")
sm.samap.adata.obs.to_csv("~/scratch/bcs/samc/vert2_axolotl.csv")

# KNN Graph
knn = sm.samap.adata.obsp['knn']
nonzero_mask = np.array(knn[knn.nonzero()] > 0)[0]
rows = knn.nonzero()[0][nonzero_mask]
cols = knn.nonzero()[1][nonzero_mask]
knn[rows, cols] = 1

# MOUSE
real = my_mapper(mz_col = 'mz_good_names', mm_col = 'mm_ClusterName')
real.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_zeisel_mapping_mine3.csv")

# Permutations: returns a list of dataframes that contain 0/1 depending on if the permutation was greater than the real value
perm_nums = 1000
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_good_names'] * perm_nums, ['mm_ClusterName'] * perm_nums))

# Find permutation p-values
perm_dist = list(itertools.chain(*perm_list))
np.quantile(perm_dist, 0.95)
max(perm_dist)
perm_dist = np.array(perm_dist)
perm_p = real.applymap(my_p)

perm_p.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_zeisel_mapping_mine_p3.csv")

# TURTLE
real = my_mapper(mz_col = 'mz_good_names', mm_col = 'cp_detail')
real.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_turtle_mapping_mine3.csv")

# Permutations: returns a list of dataframes that contain 0/1 depending on if the permutation was greater than the real value
perm_nums = 1000
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_good_names'] * perm_nums, ['cp_detail'] * perm_nums))

# Find permutation p-values
perm_dist = list(itertools.chain(*perm_list))
np.quantile(perm_dist, 0.95)
max(perm_dist)
perm_dist = np.array(perm_dist)
perm_p = real.applymap(my_p)

perm_p.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_turtle_mapping_mine_p3.csv")

# BIRD
real = my_mapper(mz_col = 'mz_good_names', mm_col = 'tg_cluster_orig2')
real.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_bird_mapping_mine3.csv")

# Permutations: returns a list of dataframes that contain 0/1 depending on if the permutation was greater than the real value
perm_nums = 1000
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_good_names'] * perm_nums, ['tg_cluster_orig2'] * perm_nums))

# Find permutation p-values
perm_dist = list(itertools.chain(*perm_list))
np.quantile(perm_dist, 0.95)
max(perm_dist)
perm_dist = np.array(perm_dist)
perm_p = real.applymap(my_p)

perm_p.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_bird_mapping_mine_p3.csv")

# AXOLOTL
real = my_mapper(mz_col = 'mz_good_names', mm_col = 'am_cluster')
real.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_axolotl_mapping_mine3.csv")

# Permutations: returns a list of dataframes that contain 0/1 depending on if the permutation was greater than the real value
perm_nums = 1000
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_good_names'] * perm_nums, ['am_cluster'] * perm_nums))

# Find permutation p-values
perm_dist = list(itertools.chain(*perm_list))
np.quantile(perm_dist, 0.95)
max(perm_dist)
perm_dist = np.array(perm_dist)
perm_p = real.applymap(my_p)

perm_p.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/vert2_axolotl_mapping_mine_p3.csv")

# import scipy
# dist_vect = scipy.spatial.distance.pdist(sm.samap.adata.obsm['X_umap'])
# dist = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(sm.samap.adata.obsm['X_umap']), checks = False) # overloads memory even at 192 GB of memory
#
# mz_idx = np.where(sm.samap.adata.obs['species'] == 'mz')[0]
# mm_idx = np.where(sm.samap.adata.obs['species'] == 'mm')[0]
# tg_idx = np.where(sm.samap.adata.obs['species'] == 'tg')[0]
# am_idx = np.where(sm.samap.adata.obs['species'] == 'am')[0]
# cp_idx = np.where(sm.samap.adata.obs['species'] == 'cp')[0]
# # dist_all = dok_matrix((len(mz_idx), sm.samap.adata.obsm['X_umap'].shape[0]-len(mz_idx)))
# dist_all = scipy.sparse.csr_matrix((len(mz_idx), sm.samap.adata.obsm['X_umap'].shape[0]-len(mz_idx)), dtype=np.int8)
#
# dist = euclidean_distance2(sm.samap.adata.obsm['X_umap'][mz_idx], sm.samap.adata.obsm['X_umap'][mm_idx])
# dist_all[mz_idx[:, None], mm_idx] = dist
# knn = dist
# real = my_mapper2(mz_col = 'mz_good_names', mm_col = 'mm_ClusterName')
# # dist = np.pad(dist, ((0,dist.shape[1]-dist.shape[0]),(0,0)), mode='constant', constant_values=0)
# # out = dist.T + dist
# # np.fill_diagonal(out,np.diag(dist))
#
#
# from scipy.sparse import dok_matrix
# from scipy.spatial import distance
# import pandas as pd
# from numba import njit, prange
# import numpy as np
#
# @njit(parallel=True)
# def euclidean_distance(coords1, coords2):
#     # allocate output array
#     c1_length, c2_length = len(coords1), len(coords2)
#     out = np.empty(shape=(c1_length, c2_length), dtype=np.float64)
#     for lat_ix in prange(c1_length):
#         for lon_ix in prange(c2_length):
#             if lat_ix >= lon_ix: # do the reverse for the upper triangle
#                 out[lat_ix, lon_ix] = (
#                     (coords1[lat_ix, 0] - coords2[lon_ix, 0]) ** 2
#                     + (coords1[lat_ix, 1] - coords2[lon_ix, 1]) ** 2
#                 ) ** 0.5
#             else:
#                 out[lat_ix, lon_ix] = 0
#     return out
#
#
# @njit(parallel=True)
# def euclidean_distance2(coords1, coords2):
#     # allocate output array
#     c1_length, c2_length = len(coords1), len(coords2)
#     out = np.empty(shape=(c1_length, c2_length), dtype=np.float64)
#     for lat_ix in prange(c1_length):
#         for lon_ix in prange(c2_length):
#             out[lat_ix, lon_ix] = (
#                 (coords1[lat_ix, 0] - coords2[lon_ix, 0]) ** 2
#                 + (coords1[lat_ix, 1] - coords2[lon_ix, 1]) ** 2
#             ) ** 0.5
#     return out
#
# def square_to_condensed(i, j, n):
#     assert i != j, "no diagonal elements in condensed matrix"
#     if i < j:
#         i, j = j, i
#     return n*j - j*(j+1)//2 + i - 1 - j
#
# c = list(itertools.product(np.where(sm.samap.adata.obs['mz_good_names'] == "8.1_Glut")[0], np.where(sm.samap.adata.obs['mm_ClusterName'] == "TEGLU6")[0]))
# d, e = zip(*c)
# this_idx = square_to_condensed(d, e, sm.samap.adata.obsm['X_umap'].shape[0])
# with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
#     tmp = mp_pool.starmap(square_to_condensed, zip(d, e, sm.samap.adata.obsm['X_umap'].shape[0] * len(d)))
#
# this_idx = square_to_condensed(np.where(sm.samap.adata.obs['mz_good_names'] == "8.1_Glut"), np.where(sm.samap.adata.obs['mm_ClusterName'] == "TEGLU6"), sm.samap.adata.obsm['X_umap'].shape[0])
#
# def my_mapper2(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
#     meta = sm.samap.adata.obs
#     all_mz_cluster = meta[mz_col].unique()
#     all_mm_cluster = meta[mm_col].unique()
#     all_mz_cluster.sort()
#     all_mm_cluster.sort()
#     lb = LabelBinarizer(sparse_output=True)
#     knn_mz  = lb.fit_transform(sm.samap.adata.obs.loc[sm.samap.adata.obs['species'] == 'mz', mz_col]).T.dot(knn)
#     knn_sum = pd.DataFrame(lb.fit_transform(sm.samap.adata.obs.loc[sm.samap.adata.obs['species'] == 'mm', mm_col]).T.dot(knn_mz.T))
#     mz_counts = meta[mz_col].value_counts().sort_index()
#     mm_counts = meta[mm_col].value_counts().sort_index()
#     mz_mm_counts = pd.DataFrame(np.array(mz_counts) * np.array(mm_counts).reshape(len(mm_counts), 1))
#     knn_mean = knn_sum / mz_mm_counts
#     knn_mean.index = all_mm_cluster
#     knn_mean.columns = all_mz_cluster
#     knn_mean = knn_mean.drop("unassigned", axis=0)
#     knn_mean = knn_mean.drop("unassigned", axis=1)
#     return(knn_mean)

# def run_euc(point):
#     list_a = sm.samap.adata.obsm['X_umap'][point]
#     list_b = sm.samap.adata.obsm['X_umap']
#     return np.array([[ np.linalg.norm(i-j) for j in list_b] for i in list_a])
# 
# with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
#     dist_list = mp_pool.map(run_euc2, range(0, sm.samap.adata.obsm['X_umap'].shape[0]))
# 
# sm.samap.adata.obsm['X_umap']
# 
# def run_euc2(point):
#     return(scipy.spatial.distance.cdist([sm.samap.adata.obsm['X_umap'][point]], sm.samap.adata.obsm['X_umap']))
# 
# import scipy
# dist = scipy.spatial.distance.cdist(sm.samap.adata.obsm['X_umap'][1:5], sm.samap.adata.obsm['X_umap'][1:5])
# dist = scipy.spatial.distance.pdist(sm.samap.adata.obsm['X_umap'][1:5])
# 
# my_list = list(range(0, sm.samap.adata.obsm['X_umap'].shape[0]))
# n=10000
# chunks = [my_list[i * n:(i + 1) * n] for i in range((len(my_list) + n - 1) // n )]
# 
# def run_euc3(chunk_idx):
#     return(scipy.spatial.distance.cdist(sm.samap.adata.obsm['X_umap'][chunks[chunk_idx]], sm.samap.adata.obsm['X_umap']))
# 
# 
# with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
#     dist_list = mp_pool.map(run_euc3, range(0, len(chunks)))
# 
# def dist_1(mat,vec):
#     res=np.empty(mat.shape[0],dtype=mat.dtype)
#     for i in nb.prange(mat.shape[0]):
#         acc=0
#         for j in range(mat.shape[1]):
#             acc+=(mat[i,j]-vec[j])**2
#         res[i]=np.sqrt(acc)
#     return res