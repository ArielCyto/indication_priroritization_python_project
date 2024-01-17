import script_utils # where the functions are located

# Check different cases
script_utils.rank_groups_by_genes(['ABCA1', 'ABCA2', 'ABCA3', 'NAT2'])
script_utils.rank_groups_by_genes(['ABCA1', 'ABCA2', 'NAT2'], 'CI_high')
script_utils.rank_groups_by_genes(['ABCA1', 'ABCA2', 'NAT']) # should fail for no NAT gene
script_utils.rank_groups_by_genes(['ABAT', 'ADCY6', 'NAT2'], 'median')


