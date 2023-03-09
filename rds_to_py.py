import scanpy as sc
import pandas as pd
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='Find the correlation')
    parser.add_argument("-i", "--input_path", help="Input File Path")
    args = parser.parse_args()
    return args.input_path

def main():
    input_path = parseArgs()
    base_name = input_path.split(".")[0]
    adata = sc.read(base_name + ".h5ad")

    adata.var = pd.DataFrame({"gene":list(adata.var['features'])}, index = list(adata.var['features']))
    adata.var_names = adata.var['gene']
    tempAdata = adata.raw.to_adata()
    tempAdata.var = adata.var
    tempAdata.var_names = adata.var_names
    adata.raw = tempAdata

    adata.write_h5ad(base_name + ".h5ad")

if __name__ == '__main__':
    main()
