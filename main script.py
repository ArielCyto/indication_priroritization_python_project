import numpy as np
from pycytocc.client import get_outputs
import pybigquery
import pandas as pd


# load input data
expression_matrix = pd.read_csv("algorithm_data/small_expression_matrix.csv",index_col=[0])
annotation_table = pd.read_csv("algorithm_data/small_annotation_table.csv",index_col=[0])

# define effector genes by the user
effector_genes = ['ABCA1',	'ABCA2', 'ABCA3', 'NAT2']

# check if the genes are available:
submit = False
# for gene in effector_genes:
   # gene_status = (gene in name(expression_matrix))



