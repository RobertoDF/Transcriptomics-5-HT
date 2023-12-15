import os
os.chdir("..")
import pandas as pd
from Utils.Settings import output_folder_calculations, order_broad_division, class_to_division, class_to_broad_division, genes_families, neuron_cluster_groups_order, manifest, download_base, output_folder, output_folder_supplementary, family_name, threshold_expression
from Utils.Utils import percentage_above_threshold
import numpy as np

metadata = manifest['file_listing']['WMB-10X']['metadata']

rpath = metadata['cell_metadata_with_cluster_annotation']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath )
cell = pd.read_csv(file, keep_default_na=False)
cell.set_index('cell_label',inplace=True)

matrices = cell.groupby(['dataset_label','feature_matrix_label'])[['library_label']].count()
matrices.columns  = ['cell_count']


rpath = metadata['example_genes_all_cells_expression']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath.replace('example_genes_all_cells_expression.csv', f'{family_name}_genes_all_cells_expression.csv'))
exp = pd.read_csv(file)
exp.set_index('cell_label',inplace=True)
exp = exp.sort_index(axis=1)
expression_total = round((((exp>threshold_expression).sum(axis=0)/exp.shape[0])*100),2)

cell["division"] = cell['class'].map(class_to_division)
cell["broad_division"] = cell['class'].map(class_to_broad_division)

metadata = manifest['file_listing']['WMB-neighborhoods']['metadata']
rpath = metadata['cluster_group_membership']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
group_membership = pd.read_csv(file) # cluster can belong to two different groups
metadata = manifest['file_listing']['WMB-taxonomy']['metadata']
rpath = metadata['cluster_annotation_term_set']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
term_set = pd.read_csv(file)
metadata = manifest['file_listing']['WMB-neighborhoods']['metadata']
rpath = metadata['dimension_reduction']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
maps = pd.read_csv(file)
maps.set_index('name',inplace=True)

cell_with_membership = cell.reset_index().merge(group_membership[["cluster_group_label", "cluster_alias", "cluster_group_name"]], on='cluster_alias').set_index("cell_label")

selected_genes = exp.columns.sort_values()

joined = cell.join(exp)
joined_boolean =  cell.join( exp>threshold_expression  )
subsampled = joined.loc[::30]

joined_with_membership = cell_with_membership.join(exp)

joined_boolean_with_membership =  cell_with_membership.join( exp>threshold_expression  )
subsampled_with_membership = joined_with_membership.loc[::30]

total_cells_by_neurotransmitter = round(joined.groupby(["neurotransmitter"]).size()/joined.shape[0]*100,2)

expression_by_neurotransmitter = round(joined.groupby("neurotransmitter")[selected_genes].apply(percentage_above_threshold), 2)

expression_by_group = round(joined_with_membership.groupby("cluster_group_name")[selected_genes].apply(percentage_above_threshold), 2)

expression_by_class = round(joined_with_membership.groupby("class")[selected_genes].apply(percentage_above_threshold), 2)
# Assuming 'expression_by_neurotransmitter' is a DataFrame
correlation_matrix = expression_by_neurotransmitter.corr()

# Mask to get the lower triangle, excluding the diagonal
mask = np.tril(np.ones_like(correlation_matrix, dtype=bool), k=-1)

lower_triangle_neurotransmitter = correlation_matrix.where(mask)

# Assuming 'expression_by_neurotransmitter' is a DataFrame
correlation_matrix = expression_by_group.corr()

# Mask to get the lower triangle, excluding the diagonal
mask = np.tril(np.ones_like(correlation_matrix, dtype=bool), k=-1)

lower_triangle_groups = correlation_matrix.where(mask)

correlation_matrix = expression_by_class .corr()

# Mask to get the lower triangle, excluding the diagonal
mask = np.tril(np.ones_like(correlation_matrix, dtype=bool), k=-1)

lower_triangle_class = correlation_matrix.where(mask)

mean_expression_by_class = round(joined.groupby("class")[selected_genes].apply(percentage_above_threshold).mean(), 2)


correlation_matrix = exp.corr()

# Mask to get the lower triangle, excluding the diagonal
mask = np.tril(np.ones_like(correlation_matrix, dtype=bool), k=-1)

lower_triangle_corr = round(correlation_matrix.where(mask), 2).stack().sort_values(ascending=False)

coloc = pd.read_pickle(f"{output_folder_calculations}/total_colocalization_{family_name}.pkl")
coloc.rename(columns={"Value":"Co-localization"}, inplace=True)
coloc['Gene1'] = pd.Categorical(coloc['Gene1'], categories=selected_genes, ordered=True)

coloc['Gene2'] = pd.Categorical(coloc['Gene2'], categories=selected_genes, ordered=True)
coloc = coloc[coloc["Co-localization"]<100]
mean_coloc = round(coloc.groupby("Gene2")["Co-localization"].mean(), 2)

_ = joined_boolean[selected_genes]
at_least_2_receptors = {}
for gene in selected_genes:
    at_least_2_receptors[gene] = (np.sum(_[_[gene]==True].sum(axis=1)>=2)/_[_[gene]==True].shape[0])*100

at_least_2_receptors = pd.DataFrame.from_dict(at_least_2_receptors, orient="index", columns=["Percentage co-localized (%)"])

at_least_3_receptors = {}
for gene in selected_genes:
    at_least_3_receptors[gene] = (np.sum(_[_[gene]==True].sum(axis=1)>=3)/_[_[gene]==True].shape[0])*100

at_least_4_receptors = {}
for gene in selected_genes:
    at_least_4_receptors[gene] = (np.sum(_[_[gene]==True].sum(axis=1)>=4)/_[_[gene]==True].shape[0])*100

at_least_5_receptors = {}
for gene in selected_genes:
    at_least_5_receptors[gene] = (np.sum(_[_[gene]==True].sum(axis=1)>=5)/_[_[gene]==True].shape[0])*100

at_least_5_receptors = pd.DataFrame.from_dict(at_least_5_receptors, orient="index", columns=["Percentage co-localized (%)"])

at_least_2_receptors_per_group = {}
for area in _['cluster_group_name'].unique():
    at_least_2_receptors_per_group[area] = (np.sum(_[_['cluster_group_name']==area][selected_genes].sum(axis=1)>=2)/_[_['cluster_group_name']==area].shape[0])*100

at_least_2_receptors_per_group = pd.DataFrame.from_dict(at_least_2_receptors_per_group, orient="index", columns=["Percentage co-localized (%)"])


