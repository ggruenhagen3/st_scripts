# Load packages
import pickle
import scipy
import scipy.sparse as sparse
import pandas
import numpy as np
import time
import random
import argparse
import h5py
import multiprocessing
from itertools import repeat

global data_mat
global gene_labels
global cond_labels

def parseArgs():
    parser = argparse.ArgumentParser(description='Find the correlation')
    parser.add_argument("-o", "--output_folder", help="Output Folder", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/st/results/",
                        const="/storage/home/hcoda1/6/ggruenhagen3/scratch/st/results/")
    parser.add_argument("-a", "--do_abs", help="Take the absolute value of the correlations?", action="store_true")
    args = parser.parse_args()
    return args.output_folder, args.do_abs

def corOnlyAndWrite(this_idx, output_path):
    """
    Given idexes of cells, create a matrix and find correlations only
    :param this_idx: Indexes of columns
    :param output_path: Output path of h5 correlation matrix file
    :return success: Function completed? True/False
    """
    cor = pandas.DataFrame(data=sparse_corrcoef(data_mat[:, this_idx].todense()), index=gene_labels, columns=gene_labels)
    if do_abs:
        print("Taking absolute value of correlations")
        cor = cor.abs()
    else:
        print("NOT taking absolute value of correlations. Using raw values.")
    h5f = h5py.File(output_path, 'w')
    h5f.create_dataset('name', data=cor)
    h5f.close()
    return True

def sparse_corrcoef(A, B=None):
    """
    Find correlations in sparse matrix
    """
    if B is not None:
        A = sparse.vstack((A, B), format='csr')
    A = A.astype(np.float64)
    n = A.shape[1]
    # Compute the covariance matrix
    rowsum = A.sum(1)
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)
    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    coeffs = C / np.sqrt(np.outer(d, d))
    return coeffs

def main():
    # Start the timer
    start_time = time.perf_counter()

    # Read Inputs
    global do_abs
    output_folder, do_abs = parseArgs()

    # Make Sparse Matrix
    # import pandas
    # import scipy.sparse as sparse
    # test = pandas.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_sct_data.csv', index_col = 0)
    # # test.iloc[0::, 1::]
    # test_sparse = sparse.csr_matrix(test)
    # sparse.save_npz('/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_sct_data.npz', test_sparse)

    # Read BB data
    global data_mat
    global gene_labels
    global cond_labels
    data_mat = sparse.load_npz("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_sct_data.npz")
    gene_labels = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_sct_data_names.csv").iloc[:, 1].to_numpy()

    # Change file name based on input
    base_name = "real"
    if do_abs:
        base_name = base_name + "_abs"

    # Find Correlations
    print("Finding Correlations")
    base_name = base_name + "_cor"
    output_path = output_folder + "/" + base_name + ".h5"
    all_cor_success = corOnlyAndWrite(range(0, data_mat.shape[1]), output_path)

if __name__ == '__main__':
    main()
