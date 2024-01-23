import pandas as pd
import numpy as np
from scipy.stats import bootstrap

# load input data - *** will be changed ***
expression_matrix = pd.read_csv("algorithm_data/small_expression_matrix.csv",index_col=[0])
annotation_table = pd.read_csv("algorithm_data/small_annotation_table.csv",index_col=[0])

# check if the requested genes are available, return submit = True if they are - *** will be changed ***
def check_input_genes(gene_list):
    # check if the genes are available:
    submit = False
    for gene in gene_list:
        gene_status = (gene in expression_matrix.columns)
        if gene_status:
            submit = True
        else:
            print("gene", gene, "is an invalid input, please remove it and try again")
            submit = False
            # check for alternative gene names
            all_genes = pd.Series(expression_matrix.columns.values)
            all_genes.index = expression_matrix.columns.values
            result = all_genes.str.contains(pat=gene.upper())
            if result.any():
                alternative_gene_names = result[np.where(result)[0]].index
                print("optional alternatives: ")
                print(alternative_gene_names.to_list())
            break
    return submit

# calculate average expression per sample by a given gene_list
# def calc_avg(gene_list):
#     effector_genes_expression_matrix = expression_matrix[expression_matrix.columns.intersection(gene_list)]
#     return(effector_genes_expression_matrix.mean(axis=1))

def calc_median(gene_list):
    effector_genes_expression_matrix = expression_matrix[expression_matrix.columns.intersection(gene_list)]
    return(effector_genes_expression_matrix.median(axis=1))

# calculate median, CI_low and CI_high for given and return the results in df
def calc_statistics(annotation_table):
    groups = annotation_table['group'].unique()
    groups_df = pd.DataFrame(index=[groups], columns=['median', 'CI_low', 'CI_high'])

    # calculate the median average expression per group
    for group in groups:
        group_median = np.median((annotation_table.loc[annotation_table['group'] == group])['effector_genes_median'])
        groups_df.loc[group, 'median'] = group_median

    # calculate CI_low, CI_high per group
    for group in groups:
        # collect all the group's sample_id into group_samples object
        group_samples = annotation_table.loc[annotation_table['group'] == group]['sample_id']
        # collect the average expression of these samples into group_median_vector
        group_median_vector = annotation_table.loc[annotation_table['sample_id'].isin(group_samples)][
            'effector_genes_median']
        group_median_vector = (group_median_vector.to_numpy(),)  # samples must be in a sequence
        # run bootstrapping on group_median_vector and extract the CI_low, CI_high results
        res = bootstrap(group_median_vector, np.median, confidence_level=0.95)
        CI_low, CI_high = res.confidence_interval
        groups_df.loc[group, 'CI_low'] = CI_low
        groups_df.loc[group, 'CI_high'] = CI_high
    return(groups_df)

# main function
def rank_groups_by_genes(gene_list, rank_by = 'CI_low'):
    # check if the genes are available:
    submit = check_input_genes(gene_list)

    # run analysis
    if (submit == True):
        # per sample (row), calculate the selected effector_genes median expression
        annotation_table['effector_genes_median'] = calc_median(gene_list)

        # get a df with the calculated 'median', 'CI_low', 'CI_high' per group
        stat_df = calc_statistics(annotation_table)

        # Generate final table to display for the customer
        final_table = ((annotation_table.drop(['sample_id', 'effector_genes_median'], axis=1)).drop_duplicates())
        final_table.index = [annotation_table['group'].unique()]
        final_table['median'] = stat_df['median']
        final_table['CI_low'] = stat_df['CI_low']
        final_table['CI_high'] = stat_df['CI_high']
        # rank by rank_by
        final_table = final_table.sort_values(by=rank_by, ascending=False)
        final_table.pop("n_per_indication")
        # generate the ranking column
        final_table['rank'] = np.arange(1, len(final_table) + 1, 1)
        final_table.insert(0, "rank", final_table.pop("rank"))
        # print results
        print(final_table.to_string())
