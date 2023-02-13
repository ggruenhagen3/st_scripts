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

# TUTORIAL
# BLAST results
A=pd.read_csv('example_data/maps/plsc/pl_to_sc.txt',sep='\t',index_col=0,header=None)
B=pd.read_csv('example_data/maps/plsc/sc_to_pl.txt',sep='\t',index_col=0,header=None)
A.head()

# Raw Data
fn1 = 'example_data/planarian.h5ad'
fn2 = 'example_data/schistosome.h5ad'
fn3 = 'example_data/hydra.h5ad'

# Load Objects
filenames = {'pl':fn1,'sc':fn2,'hy':fn3}
sam1=SAM()
sam1.load_data(fn1)
sam2=SAM()
sam2.load_data(fn2)
sam3=SAM()
sam3.load_data(fn3)
sams = {'pl':sam1,'sc':sam2,'hy':sam3}
sm = SAMAP(
        sams,
        f_maps = 'example_data/maps/',
    )

# Mine
map_dir = '/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/samap_directory/mouse_mz/maps/'
# A=pd.read_csv(map_dir+'mmmz/mz_to_mm.txt',sep='\t',index_col=0,header=None)
# B=pd.read_csv(map_dir+'mmmz/mm_to_mz.txt',sep='\t',index_col=0,header=None)
# '/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_pr.h5ad' made from SAMAP
# mm_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/samap_directory/mouse_mz/mouse_protein_table.tsv', sep=None)
# mz_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/samap_directory/mouse_mz/mz_protein_table.csv')
mm_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/mm_gene_prot_table.tsv', sep=None, header=None)
mz_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/mz_gene_prot_table.tsv', sep=None, header=None)
mm_prot = mm_prot.loc[mm_prot.iloc[:,6] != "-",:]
mz_prot = mz_prot.loc[mz_prot.iloc[:,6] != "-",:]
mm_prot_tup = list(zip(mm_prot.iloc[:,5], mm_prot.iloc[:,15]))
mz_prot_tup = list(zip(mz_prot.iloc[:,5], mz_prot.iloc[:,15]))
protein_to_gene_names = {'mz': mz_prot_tup, 'mm':mm_prot_tup}
fn1 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_pr.h5ad'
fn2 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/zeisel_tel_norm_pr.h5ad'
filenames = {'mz':fn1,'mm':fn2}
sam1=SAM()
sam1.load_data(fn1)
sam2=SAM()
sam2.load_data(fn2)
sams = {'mz':sam1,'mm':sam2}
# sams = {'mz':fn1,'mm':fn2}
sm = SAMAP(sams, names=protein_to_gene_names, f_maps = map_dir)

