import pandas as pd
import numpy as np
import random


# load input data
expression_matrix = pd.read_csv("algorithm_data/small_expression_matrix.csv",index_col=[0])
annotation_table = pd.read_csv("algorithm_data/small_annotation_table.csv",index_col=[0])
expression_matrix.head()
annotation_table.head()
# define effector genes by the user
effector_genes = ['ABCA1',	'ABCA2', 'ABCA3', 'NAT2']

# check if the genes are available:
submit = False
for gene in effector_genes:
    gene_status = (gene in expression_matrix.columns)
    if gene_status:
        submit = True
    else:
       print("gene ", gene, " is an invalid input, please remove it and try again")
       submit = False
       # I need to add here an option to check for alternative names in the data
       break

# run analysis
if (submit == True):
    # per sample (row), calculate the selected effector_genes average expression
    mini_expression_matrix = expression_matrix[expression_matrix.columns.intersection(effector_genes)]
    mini_expression_matrix['effector_genes_average'] = mini_expression_matrix.mean(axis=1)

    # work with the annotation_table object of gene_exp_data with all the necessary columns
    annotation_table['effector_genes_average'] = mini_expression_matrix['effector_genes_average']

    # calculate the group median
    groups = annotation_table['group'].unique()
    groups_df = pd.DataFrame(index = [groups], columns= ['median'])
    for group in groups:
       tmp = annotation_table.loc[annotation_table['group'] == group]
       tmp_median = np.median(tmp['effector_genes_average'])
       groups_df.loc[group, 'median'] = tmp_median

    print(groups_df)


    # bootstrapping per group

   # group_median = []
   # x = np.random.normal(loc=300.0, size=1000)
   # print(np.mean(x))
   # for i in range(50):
   #     y = random.sample(x.tolist(), 4)
   #     median = np.median(y)
   #     group_median.append(median)
   # print("hi")









