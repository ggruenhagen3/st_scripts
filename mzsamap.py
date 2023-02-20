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
# keys_rs = {'mz':'seuratclusters53','mm':'region_subcluster'}
# D,MappingTable = get_mapping_scores(sm,keys_rs,n_top = 0)
# MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_saunders_subcluster_mapping_table.csv")

# sm.sams['mm'].adata.obs["region_cluster"] = sm.sams['mm'].adata.obs['region'].astype('str') + "_" + sm.sams['mm'].adata.obs['cluster'].astype('str')
# keys_rc = {'mz':'seuratclusters53','mm':'region_cluster'}
# D,MappingTable = get_mapping_scores(sm,keys_rc,n_top = 0)
# MappingTable.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/mz_mm_saunders_cluster_mapping_table.csv")

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


with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/mz_mm_oritiz_sm_ran.pkl', 'wb') as out_file:
    pickle.dump(sm, out_file)

# Read samap
# import pickle
# with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/mz_mm_saunders_sm_ran.pkl', 'rb') as out_file:
#     sm = pickle.load(out_file)

# Export as h5ad so it can be loaded into Seurat
samap.adata.write_h5ad("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/mz_mm_oritiz.h5ad")
