import numpy as np
from pycytocc.client import get_files

# load input data
expression_matrix = get_files("wf-219a0c80f1")
annotation_table = get_files("wf-c546d82645")

# define effector genes by the user
effector_genes = ['BRCA1', 'BRCA2', 'HOXB13']

# check if the genes are available:
submit = False
# for gene in effector_genes:
   # gene_status = (gene in name(expression_matrix))



