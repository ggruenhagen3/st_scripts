# conda activate SAMap
# Imports
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from sklearn.preprocessing import LabelBinarizer

# BLAST Results
map_dir = '/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/samap_directory/mouse_mz/maps/'
gene_info = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/gene_info_3.csv')
mm_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/mm_gene_prot_table.tsv', sep=None, header=None)
mz_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/mz_gene_prot_table.tsv', sep=None, header=None)
mz_prot["loc"] = "LOC" + mz_prot.iloc[:,1].astype(str)
mz_prot = mz_prot.merge(gene_info, how='inner', on="loc")
mm_prot = mm_prot.loc[mm_prot.iloc[:,6] != "-",:]
mz_prot = mz_prot.loc[mz_prot.iloc[:,6] != "-",:]
mm_prot_tup = list(zip(mm_prot.iloc[:,5], mm_prot.iloc[:,15]))
mz_prot_tup = list(zip(mz_prot.iloc[:,5], mz_prot["seurat_name"]))
protein_to_gene_names = {'mz': mz_prot_tup, 'mm':mm_prot_tup}

# Objects containing gene x counts
# fn1 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_pr.h5ad'
# fn2 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/zeisel_tel_norm_pr.h5ad'
# fn2 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/saunders_norm_pr.h5ad'
# fn2 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/oritiz_norm_pr.h5ad'
fn1 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5ad'
fn2 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/oritiz.h5ad'
filenames = {'mz':fn1,'mm':fn2}
sam1=SAM()
sam1.load_data(fn1)
sam2=SAM()
sam2.load_data(fn2)
# sams = {'mz':sam1,'mm':sam2}
sams = {'mz':fn1,'mm':fn2}
sm = SAMAP(sams, names=protein_to_gene_names, f_maps = map_dir)
# sm.run(neigh_from_keys={'mz':'seuratclusters53', 'mm':'ClusterName'})
sm.run()
samap = sm.samap

# Get celltype-celltype mapping #
# Cluster-wise
# keys_cluster = {'mz':'seuratclusters53','mm':'ClusterName'}
# D,MappingTable = get_mapping_scores(sm,keys_cluster,n_top = 0)
# MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_zei_mapping_table.csv")

# Zeisel
# Mouse broad clusters
# samap.adata.obs["mm_broad"] = samap.adata.obs["mm_Region"] + "_" + samap.adata.obs["mm_ClusterName"]
# sm.samap.adata.obs["mm_broad"] = sm.samap.adata.obs["mm_Region"] + "_" + sm.samap.adata.obs["mm_ClusterName"]
# sm.sams['mm'].adata.obs["broad"] = sm.sams['mm'].adata.obs['Region'].astype('str') + "_" + sm.sams['mm'].adata.obs['TaxonomyRank4'].astype('str')
# keys_broad = {'mz':'seuratclusters53','mm':'broad'}
# D,MappingTable = get_mapping_scores(sm,keys_broad,n_top = 0)
# MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_zei_broad_mapping_table.csv")

# Gene pairs that contribute positively to the correlation
# gpf_cluster = GenePairFinder(sm,keys=keys_cluster)
# gene_pairs_cluster = gpf_cluster.find_all(align_thr=0.2)
# gene_pairs_cluster.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_zei_gene_pair.csv")

# gpf_broad = GenePairFinder(sm,keys=keys_broad)
# gene_pairs_broad = gpf_broad.find_all(align_thr=0.2)
# gene_pairs_broad.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_zei_broad_gene_pair.csv")

# Saunders
# Mouse subclusters
# sm.sams['mm'].adata.obs["region_subcluster"] = sm.sams['mm'].adata.obs['region'].astype('str') + "_" + sm.sams['mm'].adata.obs['subcluster'].astype('str')
# keys_rs = {'mz':'seuratclusters53','mm':'region_suplt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/tmp.png")bcluster'}
# D,MappingTable = get_mapping_scores(sm,keys_rs,n_top = 0)
# MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_saunders_subcluster_mapping_table.csv")

# sm.sams['mm'].adata.obs["region_cluster"] = sm.sams['mm'].adata.obs['region'].astype('str') + "_" + sm.sams['mm'].adata.obs['cluster'].astype('str')
# keys_rc = {'mz':'seuratclusters53','mm':'region_cluster'}
# D,MappingTable = get_mapping_scores(sm,keys_rc,n_top = 0)
# MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_saunders_cluster_mapping_table.csv")

# Oritiz
keys_cluster = {'mz':'seuratclusters53','mm':'ABA_parent'}
D,MappingTable = get_mapping_scores(sm,keys_cluster,n_top = 0)
MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_oritz_mapping_table.csv")

keys_cluster = {'mz':'structure','mm':'ABA_parent'}
D,MappingTable = get_mapping_scores(sm,keys_cluster,n_top = 0)
MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_st_oritz_mapping_table.csv")

# sm.samap.adata.obs['mz_struct'] = sm.samap.adata.obs['mz_struct_b2_vdc']
# sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_struct'] == "Dl-d", 'mz_struct'] = "Dl-v"
# sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_struct'] == "SP-u", 'mz_struct'] = "Vx"
# sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_struct'] == "tract", 'mz_struct'] = "ON"
my_map = my_mapper(mz_col = 'mz_good_names', mm_col = 'mm_b_parent', mode = "mz_to_mm")
my_map = my_map[my_map.columns].astype(float)
sns.clustermap(my_map, annot=False, cmap = "viridis", figsize = [my_map.shape[1]/2, my_map.shape[0]/2], z_score=1)
plt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_oritz_b_mapping_mztomm_z.png")
my_map.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_oritz_b_mapping_mine.csv")

knn = sm.samap.adata.obsp['knn']
knn = sm.samap.adata.obsp['connectivities']
nonzero_mask = np.array(knn[knn.nonzero()] > 0)[0]
rows = knn.nonzero()[0][nonzero_mask]
cols = knn.nonzero()[1][nonzero_mask]
knn[rows, cols] = 1

real = my_mapper(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_b_parent')

perm_nums = 1000
import multiprocessing
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_struct_b2_vdc'] * perm_nums, ['mm_b_parent'] * perm_nums))

perm_dist = list(itertools.chain(*perm_list))
np.quantile(perm_dist, 0.95)
max(perm_dist)
perm_dist = np.array(perm_dist)
perm_p = real.applymap(my_p)
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


def p_col(x):
    all_col_values = list()
    for i in range(0, len(perm_list)):
        all_col_values.extend(perm_list[i].iloc[:,x].tolist())
    ge_mat = np.greater(np.array([all_col_values]), np.array([real.iloc[:,x]]).T)
    p = ge_mat.sum(axis=1) / (len(perm_list)*real.shape[0])
    return(p)

with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    p_list = mp_pool.map(p_col, range(0, real.shape[1]))

p_perm = pd.DataFrame(p_list, columns = real.index, index = real.columns).T


# Chord celltype-celltype mapping visulizations
# from matplotlib import pyplot as plt
# import holoviews as hv
# skp = sankey_plot(MappingTable, align_thr=0.05, species_order = ['mz', 'mm'])
# hv.save(skp, '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_zei_mapping.html')
# chp = chord_plot(MappingTable, align_thr=0.05)
# hv.save(chp, '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_zei_mapping_chord.html')



# Save samap
import pickle
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/mz_mm_oritiz_samap.pkl', 'wb') as out_file:
    pickle.dump(samap, out_file)


with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_saunders_sm_ran.pkl', 'wb') as out_file:
    pickle.dump(sm, out_file)

# Read samap
import pickle
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_zeisel_sm_ran.pkl', 'rb') as out_file:
    sm = pickle.load(out_file)

with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/st_.pkl', 'rb') as out_file:
    sm = pickle.load(out_file)

# Export as h5ad so it can be loaded into Seurat
samap.adata.write_h5ad("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/mz_mm_oritiz.h5ad")

turtle_neurons_meta = pd.read_csv("~/scratch/bcs/data/turtle_neurons_meta.csv", index_col = 0)
turtle_neurons_meta = turtle_neurons_meta.loc[turtle_neurons_meta.index.isin(sm.samap.adata.obs.index),]
sm.samap.adata.obs['cp_detail'] = sm.samap.adata.obs['cp_cluster']
sm.samap.adata.obs.loc[turtle_neurons_meta.index, 'cp_detail'] = turtle_neurons_meta['clusters']

#
sm.samap.leiden_clustering(res=3)
sm.samap.adata.obs.to_csv("~/scratch/bcs/samc/tmp.csv")
sm.samap.adata.obs['value'] = 1
mz_both = pd.pivot_table(sm.samap.adata.obs, index=sm.samap.adata.obs['leiden_clusters'], columns=['mz_good_names'], aggfunc=sum)
mz_both.columns = mz_both.columns.droplevel(0)
mz_both = mz_both.drop("unassigned", axis=1)
mz_both = mz_both / mz_both.sum(axis=0)

mm_both = pd.pivot_table(sm.samap.adata.obs, index=sm.samap.adata.obs['leiden_clusters'], columns=['mm_ClusterName'], aggfunc=sum)
mm_both.columns = mm_both.columns.droplevel(0)
mm_both = mm_both.drop("unassigned", axis=1)
mm_both = mm_both / mm_both.sum(axis=0)

mz_mm_cor = pd.concat([mz_both, mm_both], axis=1).corr()
mz_mm_cor = mz_mm_cor.loc[mz_mm_cor.index.isin(sm.samap.adata.obs['mm_ClusterName']), mz_mm_cor.columns.isin(sm.samap.adata.obs['mz_good_names'])]
mz_mm_cor.to_csv("~/scratch/bcs/results/bb_zeisel_cluster_cor2.csv")

#
# gene_map = sm.gnnm_refined
# thresh_mask = np.array(gene_map[gene_map.nonzero()] > 0)[0]
# rows = gene_map.nonzero()[0][thresh_mask]
# cols = gene_map.nonzero()[1][thresh_mask]
gene_map = sm.gnnm_refined
rows, cols = gene_map.nonzero()
gene_map_df = pd.DataFrame({'value': gene_map[rows,cols].tolist()[0], 'mz_gene':sm.gns[rows], 'mm_gene':sm.gns[cols] })
gene_map_df = gene_map_df.loc[gene_map_df['mz_gene'].str.startswith('mz_') & ~gene_map_df['mm_gene'].str.startswith('mz_'), :]
# gene_map_df = gene_map_df.loc[gene_map_df['value'] > 0.25,]
gene_map_df.to_csv("~/scratch/bcs/data/bb_zeisel_gene_ortholog.csv")

# Experimenting with mapping and permutations
# perm_nums = 1000
# import multiprocessing
# with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
#     perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_good_names'] * perm_nums, ['mm_b_parent'] * perm_nums, ['mz_to_mm'] * perm_nums))
#
# from functools import reduce
# perm_p = reduce(lambda x, y: x.add(y, fill_value=0), perm_list)
# perm_p = perm_p / perm_nums
#
# def my_mapper_slow(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent', mode = "mutual"):
#     meta = sm.samap.adata.obs
#     all_species = meta['species'].unique()
#     all_mz_cluster = meta.loc[meta['species'] == all_species[0], mz_col].unique()
#     all_mm_cluster = meta.loc[meta['species'] == all_species[1], mm_col].unique()
#     xsim_mean = pd.DataFrame(columns=all_mz_cluster, index=all_mm_cluster)
#     for mz_cluster in all_mz_cluster:
#         for mm_cluster in all_mm_cluster:
#             if mode == "mutual":
#                 this_mat = sm.samap.adata.obsp['xsim'][np.where(meta[mz_col] == mz_cluster)[0],:]
#                 xsim_mean.loc[mm_cluster, mz_cluster] = this_mat[:,np.where(meta[mm_col] == mm_cluster)[0]].mean().astype(float)
#             elif mode == "mz_to_mm":
#                 this_mat = sm.samap.adata.obsp['knn'][np.where(meta[mz_col] == mz_cluster)[0],:]
#                 xsim_mean.loc[mm_cluster, mz_cluster] = this_mat[:,np.where(meta[mm_col] == mm_cluster)[0]].mean().astype(float)
#             elif mode == "mm_to_mz":
#                 this_mat = sm.samap.adata.obsp['knn'][:,np.where(meta[mz_col] == mz_cluster)[0]]
#                 xsim_mean.loc[mm_cluster, mz_cluster] = this_mat[np.where(meta[mm_col] == mm_cluster)[0],].mean().astype(float)
#     return(xsim_mean)
#
#
# my_map = my_mapper(mz_col = 'mz_good_names', mm_col = 'mm_ClusterName', mode = "mz_to_mm")
# my_map2 = my_mapper2(mz_col = 'mz_good_names', mm_col = 'mm_ClusterName')
# my_small =  my_map.loc[my_map2.index, my_map2.columns].astype("float")
# my_map2.iloc[0:5, 0:5], my_small.iloc[0:5, 0:5]
# my_map3 = my_mapper3(mz_col = 'mz_good_names', mm_col = 'mm_ClusterName')
# my_map4 = my_mapper3(mz_col = 'mz_good_names', mm_col = 'mm_ClusterName')
#
# from sklearn.preprocessing import LabelBinarizer
# knn = sm.samap.adata.obsp['knn']
# nonzero_mask = np.array(knn[knn.nonzero()] > 0)[0]
# rows = knn.nonzero()[0][nonzero_mask]
# cols = knn.nonzero()[1][nonzero_mask]
# knn[rows, cols] = 1
# knn3 = knn[np.where(sm.samap.adata.obs['species'] == 'mz')[0],:]
# knn3 = knn3[:,np.where(sm.samap.adata.obs['species'] == 'mm')[0]]
# def my_mapper2(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
#     meta = sm.samap.adata.obs
#     all_mz_cluster = meta[mz_col].unique()
#     all_mm_cluster = meta[mm_col].unique()
#     all_mz_cluster.sort()
#     all_mm_cluster.sort()
#     lb = LabelBinarizer(sparse_output=True)
#     knn_mz  = lb.fit_transform(sm.samap.adata.obs[mz_col]).T.dot(knn)
#     knn_sum = pd.DataFrame(lb.fit_transform(sm.samap.adata.obs[mm_col]).T.dot(knn_mz.T).todense())
#     mz_counts = meta[mz_col].value_counts().sort_index()
#     mm_counts = meta[mm_col].value_counts().sort_index()
#     mz_mm_counts = pd.DataFrame(np.array(mz_counts) * np.array(mm_counts).reshape(len(mm_counts), 1))
#     knn_mean = knn_sum / mz_mm_counts
#     knn_mean.index = all_mm_cluster
#     knn_mean.columns = all_mz_cluster
#     knn_mean = knn_mean.drop("unassigned", axis=0)
#     knn_mean = knn_mean.drop("unassigned", axis=1)
#     return(knn_mean)
#
# def my_mapper3(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
#     """
#     This doesn't work because of differences in cluster size
#     """
#     meta = sm.samap.adata.obs
#     all_mz_cluster = meta.loc[meta['species'] == 'mz', mz_col].unique()
#     all_mm_cluster = meta.loc[meta['species'] == 'mm', mm_col].unique()
#     all_mz_cluster.sort()
#     all_mm_cluster.sort()
#     lb = LabelBinarizer(sparse_output=True)
#     knn_mz  = lb.fit_transform(meta.loc[meta['species'] == 'mm', mm_col]).T.dot(knn3.T).T
#     knn_sum = pd.DataFrame(lb.fit_transform(meta.loc[meta['species'] == 'mz', mz_col]).T.dot(knn_mz).todense())
#     knn_sum.columns = all_mm_cluster
#     knn_sum.index = all_mz_cluster
#     mz_counts = meta.loc[meta['species'] == 'mz', mz_col].value_counts().sort_index()
#     mz_counts = mz_counts * 20
#     knn_mean = (knn_sum.T / mz_counts).T
#     # mm_counts = meta.loc[meta['species'] == 'mm', mm_col].value_counts().sort_index()
#     # knn_mean = knn_mean / mm_counts
#     return(knn_mean)
#
# def perm_mapper_original(seed_num, mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent', mode = "mutual"):
#     np.random.seed(seed_num)
#     sm.samap.adata.obs['mz_perm'] = sm.samap.adata.obs[mz_col]
#     sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'])
#     sm.samap.adata.obs['mm_perm'] = sm.samap.adata.obs[mm_col]
#     sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'])
#     perm_map = my_mapper('mz_perm', 'mm_perm', mode)
#     perm_map = perm_map.ge(real) * 1
#     return(perm_map)
#
# def perm_mapper2(seed_num, mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent', mode = "mutual"):
#     np.random.seed(seed_num)
#     sm.samap.adata.obs['mz_perm'] = sm.samap.adata.obs[mz_col]
#     sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'])
#     sm.samap.adata.obs['mm_perm'] = sm.samap.adata.obs[mm_col]
#     sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'])
#     perm_map = my_mapper('mz_perm', 'mm_perm', mode)
#     perm_map = perm_map.values.tolist()
#     perm_map = list(itertools.chain(*perm_map))
#     return(perm_map)
#
# perm_dist = list(itertools.chain(*perm_list))
# np.quantile(perm_dist, 0.95)
# perm_dist_df = pd.DataFrame(perm_dist, columns=['value'])
# perm_dist_df['real'] = False
# real_dist = real.values.tolist()
# real_dist = list(itertools.chain(*real_dist))
# real_dist_df = pd.DataFrame(real_dist, columns=['value'])
# real_dist_df['real'] = True
# dist_df = real_dist_df.append(perm_dist_df.iloc[1:1000,], ignore_index=True)
#
# plt.figure(figsize=[6, 6])
# sns.displot(dist_df, x="value", hue="real", kind="kde", fill=True)
# plt.savefig("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/bb_zeisel_samap_perm.svg")
