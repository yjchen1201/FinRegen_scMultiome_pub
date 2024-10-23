#!/usr/bin/env python3
import os
import sys
import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
import arboreto
import dask
from dask.distributed import Client
from distributed import LocalCluster, Client
import numpy as np


parser = argparse.ArgumentParser(description='GRNboost2 process gene expression matrices.')
parser.add_argument('celltype', type=str, help='Cell type to process')

# Parse the arguments
args = parser.parse_args()

celltype = args.celltype

stages = ["0dpa", "1dpa", "2dpa", "4dpa", "6dpa"]

for stage in stages:
    #ex_matrix = pd.read_csv(f'{celltype}_{stage}_gene_exp_matrix.csv')
    ex_matrix = pd.read_csv(f'{celltype}_{stage}_gene_exp_matrix_transposed.csv')

    # Filter out columns with all zeros
    col_sums = ex_matrix.sum(axis=0)
    ex_matrix = ex_matrix.loc[:, col_sums != 0]

    # Load TF names and filter them based on presence in ex_matrix
    tf_names = load_tf_names("GRNs_TF.txt")
    tf_names = [tf for tf in tf_names if tf in ex_matrix.columns]

    # Setup Dask cluster and client
    local_cluster = LocalCluster(n_workers=12, threads_per_worker=1)
    custom_client = Client(local_cluster)

    # Run GRNBoost2
    network = grnboost2(expression_data=ex_matrix, tf_names=tf_names, seed=1234, client_or_address=custom_client)

    # Save the network to a TSV file
    network.to_csv(f'{celltype}_{stage}_grnboost2_network.tsv', sep='\t', header=False, index=False)
    print(f'{celltype} {stage} network inference finished!')
custom_client.close()
local_cluster.close()
