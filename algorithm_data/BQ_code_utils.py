import pandas as pd
BQ_ndarray = pd.read_csv("algorithm_data/BQ_try.csv",index_col=[0])

def generate_expression_matrix(effector_genes):
    # # subset the table that I need
    effector_genes_data = BQ_ndarray.loc[BQ_ndarray.index.isin(effector_genes)]

    # extract samples
    sample_ids = effector_genes_data['sample_ids']
    sample_ids = sample_ids.to_list()
    sample_final_list = []
    for cell in sample_ids:
        samples = cell.split("' ")
        new_list = [x for x in samples if x != '']
        new_list = [s.replace("'", "") for s in new_list]
        new_list = [s.replace("]", "") for s in new_list]
        new_list = [s.replace("[", "") for s in new_list]
        new_list = [x.strip() for x in new_list]
        sample_final_list.extend(new_list)
    sample_final_list = list(dict.fromkeys(sample_final_list))

    new_expression_matrix = pd.DataFrame(index=[sample_final_list], columns=[effector_genes])

    # extract values
    for gene in effector_genes:
        gene_values = effector_genes_data.loc[effector_genes_data.index == gene]['values']
        gene_values = gene_values.to_list()
        values_final_list = []
        for cell in gene_values:
            gene_vals = cell.split(" ")
            new_list = [x for x in gene_vals if x != '']
            new_list = [x.strip() for x in new_list]
            new_list = [x for x in new_list if x != ']']
            new_list = [x for x in new_list if x != '[']
            new_list = [s.replace("'", "") for s in new_list]
            new_list = [s.replace("]", "") for s in new_list]
            new_list = [s.replace("[", "") for s in new_list]
            values_final_list.extend(new_list)
        new_expression_matrix[gene] = [float(i) for i in values_final_list]
    # fix the structure of the new_expression_matrix
    new_expression_matrix.index.name = 'sample_id'
    new_expression_matrix.reset_index(inplace=True)
    new_expression_matrix['sample_id'] = new_expression_matrix.loc[:, 'level_0']
    new_expression_matrix.set_index(['level_0'], inplace = True)

    return(new_expression_matrix)


