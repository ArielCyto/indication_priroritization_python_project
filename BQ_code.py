import pandas as pd
import numpy as np
from scipy.stats import bootstrap
from algorithm_data.BQ_code_utils import generate_expression_matrix
BQ_ndarray = pd.read_csv("algorithm_data/BQ_try.csv",index_col=[0])
# ask for effector genes
effector_genes = ['NAT2', 'ABCA2', 'ABCA3']


# check if the genes are available:
submit = False
for gene in effector_genes:
    gene_status = (gene in (BQ_ndarray.index.unique().values))
    if gene_status:
        submit = True
    else:
       print("gene", gene, "is an invalid input, please remove it and try again")
       submit = False
       # check for alternative gene names
       all_genes = pd.Series(BQ_ndarray.index.unique().values)
       all_genes.index = all_genes
       result = all_genes.str.contains(pat=gene.upper())
       if result.any():
           alternative_gene_names = result[np.where(result)[0]].index
           print("optional alternatives: ")
           print(alternative_gene_names.to_list())
       break
#
#
# # run analysis
if (submit == True):
    expression_matrix = generate_expression_matrix(effector_genes)  # assume all the genes are there, need to add message for missing genes
    expression_matrix['effector_genes_average'] = expression_matrix[effector_genes].mean(axis=1)



