import os
os.chdir("..")
from Utils.Utils import percentage_above_threshold, percentage_above_threshold_MER
import time
import anndata
from pathlib import Path
from Utils.Settings import class_to_division, class_to_broad_division, output_folder_calculations, manifest, download_base, \
    family_name, threshold_expression, threshold_enriched_clusters
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
import pickle
import shap
from sklearn.metrics import confusion_matrix

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

cell_with_membership = cell_with_membership[cell_with_membership['cluster_group_name']!="WholeBrain"]

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


# Convert the DataFrame to a NumPy array
corr_matrix_np = expression_by_neurotransmitter.T.corr().to_numpy()

# Replace the diagonal elements with NaN
np.fill_diagonal(corr_matrix_np, np.nan)

# Convert the NumPy array back to a DataFrame if needed
corr_by_neurotransmitter = pd.DataFrame(corr_matrix_np, index=expression_by_neurotransmitter.index, columns=expression_by_neurotransmitter.index)


# Convert the DataFrame to a NumPy array
corr_matrix_np = expression_by_group.T.corr().to_numpy()

# Replace the diagonal elements with NaN
np.fill_diagonal(corr_matrix_np, np.nan)

# Convert the NumPy array back to a DataFrame if needed
corr_by_group = pd.DataFrame(corr_matrix_np, index=expression_by_group.index, columns=expression_by_group.index)



# Convert the DataFrame to a NumPy array
corr_matrix_np = expression_by_class.T.corr().to_numpy()

# Replace the diagonal elements with NaN
np.fill_diagonal(corr_matrix_np, np.nan)

# Convert the NumPy array back to a DataFrame if needed
corr_by_class = pd.DataFrame(corr_matrix_np, index=expression_by_class.index, columns=expression_by_class.index)


# Convert the DataFrame to a NumPy array
corr_matrix_np = expression_by_class.T.corr().to_numpy()

# Replace the diagonal elements with NaN
np.fill_diagonal(corr_matrix_np, np.nan)

# Convert the NumPy array back to a DataFrame if needed
corr_by_class = pd.DataFrame(corr_matrix_np, index=expression_by_class.index, columns=expression_by_class.index)

# Convert the DataFrame to a NumPy array
corr_matrix_np = exp.corr().to_numpy()

# Replace the diagonal elements with NaN
np.fill_diagonal(corr_matrix_np, np.nan)

# Convert the NumPy array back to a DataFrame if needed
corr_by_cell = pd.DataFrame(corr_matrix_np, index=exp.corr().index, columns=exp.corr().index)



df_melted = corr_by_cell.reset_index().melt(id_vars='index', var_name='Gene2', value_name='Value')
df_melted['Gene_Pair'] = df_melted['index'] + '-' + df_melted['Gene2']
corr_by_cell = df_melted[['Gene_Pair', 'Value']]
corr_by_cell["Value"] = round(corr_by_cell["Value"], 3)
corr_by_cell.set_index("Gene_Pair", inplace=True)
corr_by_cell = corr_by_cell.sort_values("Value", ascending=False)
corr_by_cell = corr_by_cell.iloc[::2]


mean_expression_by_class = round(joined.groupby("class")[selected_genes].apply(percentage_above_threshold).mean(), 2)

sem_expression_by_class = round(joined.groupby("class")[selected_genes].apply(percentage_above_threshold).sem(), 2)


coloc = pd.read_pickle(f"{output_folder_calculations}/total_colocalization_{family_name}.pkl")

with open(f"{output_folder_calculations}/total_colocalization_{family_name}_by_neigh.pkl", "rb") as pkl_file:
    coloc_by_neigh = pickle.load(pkl_file)

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

_ = joined_boolean_with_membership[joined_boolean_with_membership['cluster_group_name']!="WholeBrain"]

at_least_2_receptors_per_group = {}
for area in _['cluster_group_name'].unique():
    at_least_2_receptors_per_group[area] = (np.sum(_[_['cluster_group_name']==area][selected_genes].sum(axis=1)>=2)/_[_['cluster_group_name']==area].shape[0])*100

at_least_2_receptors_per_group = pd.DataFrame.from_dict(at_least_2_receptors_per_group, orient="index", columns=["Percentage co-localized (%)"])


## Load MERFISH

datasets = ['Zhuang-ABCA-1']

metadata = {}
for d in datasets :
    metadata[d] = manifest['file_listing'][d]['metadata']

cell = {}

for d in datasets:
    rpath = metadata[d]['cell_metadata']['files']['csv']['relative_path']
    file = os.path.join(download_base, rpath)
    cell[d] = pd.read_csv(file, dtype={"cell_label": str})
    cell[d].set_index('cell_label', inplace=True)

    sdf = cell[d].groupby('brain_section_label')

    print(d, ":", "Number of cells = ", len(cell[d]), ", ", "Number of sections =", len(sdf))
taxonomy_metadata = manifest['file_listing']['WMB-taxonomy']['metadata']

rpath = taxonomy_metadata['cluster_to_cluster_annotation_membership_pivoted']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
cluster_details = pd.read_csv(file,keep_default_na=False)
cluster_details.set_index('cluster_alias', inplace=True)

rpath = taxonomy_metadata['cluster_to_cluster_annotation_membership_color']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
cluster_colors = pd.read_csv(file)
cluster_colors.set_index('cluster_alias', inplace=True)

cell_extended = {}

for d in datasets :
    cell_extended[d] = cell[d].join(cluster_details,on='cluster_alias')
    cell_extended[d] = cell_extended[d].join(cluster_colors,on='cluster_alias')
rpath = metadata[datasets[0]]['gene']['files']['csv']['relative_path'] # first dataset is big
file = os.path.join( download_base, rpath)
gene = pd.read_csv(file)
gene.set_index('gene_identifier',inplace=True)
print("Number of genes = ", len(gene))


expression_matrices = {}

for d in datasets :
    expression_matrices[d] = manifest['file_listing'][d]['expression_matrices']


gene_filtered  = pd.read_csv(Path(output_folder_calculations, "selected_genes_MERFISH.csv")).set_index("gene_identifier")
selected_genes = np.sort(gene_filtered["gene_symbol"].values)
selected_genes

cell_expression = {}
cell_expression_bool = {}

for d in datasets:
    expression_matrices[d]
    rpath = expression_matrices[d][d]['log2']['files']['h5ad']['relative_path']
    file = os.path.join(download_base, rpath)

    adata = anndata.read_h5ad(file, backed='r')

    start = time.process_time()
    gdata = adata[:, gene_filtered.index].to_df()
    gdata.columns = gene_filtered.gene_symbol
    cell_expression[d] = cell_extended[d].join(gdata)
    cell_expression_bool[d] = cell_extended[d].join(gdata.astype("bool"))

    print(d, "-", "time taken: ", time.process_time() - start)

    adata.file.close()
    del adata
ccf_coordinates = {}

for d in datasets:
    rpath = manifest['file_listing'][d + '-CCF']['metadata']['ccf_coordinates']['files']['csv']['relative_path']
    file = os.path.join(download_base, rpath)
    ccf_coordinates[d] = pd.read_csv(file)
    ccf_coordinates[d].set_index('cell_label', inplace=True)
    ccf_coordinates[d].rename(columns={'x': 'x_ccf', 'y': 'y_ccf', 'z': 'z_ccf'}, inplace=True)
    cell_expression[d] = cell_expression[d].join(ccf_coordinates[d], how='inner')
    cell_expression_bool[d] = cell_expression_bool[d].join(ccf_coordinates[d], how='inner')

metadata_par = manifest['file_listing']['Allen-CCF-2020']['metadata']
rpath = metadata_par['parcellation_to_parcellation_term_membership_acronym']['files']['csv']['relative_path']
file = os.path.join(download_base, rpath)
parcellation_annotation = pd.read_csv(file)
parcellation_annotation.set_index('parcellation_index', inplace=True)
parcellation_annotation.columns = ['parcellation_%s' % x for x in parcellation_annotation.columns]

rpath = metadata_par['parcellation_to_parcellation_term_membership_color']['files']['csv']['relative_path']
file = os.path.join(download_base, rpath)
parcellation_color = pd.read_csv(file)
parcellation_color.set_index('parcellation_index', inplace=True)
parcellation_color.columns = ['parcellation_%s' % x for x in parcellation_color.columns]

for d in datasets:
    cell_expression[d] = cell_expression[d].join(parcellation_annotation, on='parcellation_index')
    cell_expression[d] = cell_expression[d].join(parcellation_color, on='parcellation_index')
    cell_expression_bool[d] = cell_expression_bool[d].join(parcellation_annotation, on='parcellation_index')
    cell_expression_bool[d] = cell_expression_bool[d].join(parcellation_color, on='parcellation_index')

data_merfish = pd.concat(cell_expression).sort_values(by='brain_section_label')
data_merfish = data_merfish[data_merfish['parcellation_division'] != "unassigned"]

color_dict = data_merfish[['parcellation_division_color', 'parcellation_division']].drop_duplicates().set_index(
    'parcellation_division').to_dict()['parcellation_division_color']
color_dict.update(data_merfish[['parcellation_structure_color', 'parcellation_structure']].drop_duplicates().set_index(
    'parcellation_structure').to_dict()['parcellation_structure_color'])

sel = "neurotransmitter"  # "cluster_group_name"#"neurotransmitter"


def decoddddddd(joined_boolean, sel):
    df = joined_boolean[[sel] + list(selected_genes)]

    df.set_index(sel, inplace=True)

    # Assuming 'df' is your DataFrame
    X = df  # Features (Htr expression levels)
    y = df.index  # Target (neurotransmitter type)

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, stratify=y)

    # Initialize the Random Forest classifier
    # You can adjust 'n_estimators' and 'max_depth' based on your dataset and computational resources
    rf_classifier = RandomForestClassifier(n_estimators=10, max_depth=10, n_jobs=30, class_weight='balanced')
    # rf_classifier = HistGradientBoostingClassifier( class_weight='balanced')

    # Train the model
    rf_classifier.fit(X_train, y_train)

    # Make predictions
    y_pred = rf_classifier.predict(X_test)

    # Evaluate the model
    accuracy = accuracy_score(y_test, y_pred) * 100
    print("Accuracy:", accuracy)
    report = classification_report(y_test, y_pred, output_dict=True)
    print("\nClassification Report:\n", classification_report(y_test, y_pred))

    cm = confusion_matrix(y_test, y_pred, normalize="true")

    cm = pd.DataFrame(cm, columns=y_test.astype("category").categories,
                      index=y_test.astype("category").categories) * 100

    X_sample = X_test.sample(n=10000, weights=y_test.map(10000 / y_test.value_counts()), random_state=4)

    # Calculate SHAP values
    explainer = shap.TreeExplainer(rf_classifier, n_jobs=40)
    shap_values = explainer.shap_values(X_sample)

    out = []
    for i in range(len(shap_values)):
        out.append(pd.Series(np.abs(shap_values[i]).mean(0), name=sorted(joined_boolean[sel].unique())[i],
                             index=X_test.columns))

    shap_matrix = pd.concat(out, axis=1).T

    return cm, shap_matrix, accuracy, report


cm_neurotransmitter, shap_matrix_neurotransmitter, accuracy_neurotransmitter, report_neurotransmitter  = decoddddddd(joined_boolean, sel)

sel = "cluster_group_name"#"class"#"cluster_group_name"#"neurotransmitter"

cm_neighborood, shap_matrix_neighborood, accuracy_neighborood, report_neighborood = decoddddddd(joined_boolean_with_membership, sel);

sel = "class"#"class"#"cluster_group_name"#"neurotransmitter"

cm_class, shap_matrix_class, accuracy_class, report_class = decoddddddd(joined_boolean, sel);

report_class = pd.DataFrame(report_class)#.loc["recall"]

recall_values = report_class.loc['recall']
selected_indices = recall_values[recall_values > 0.4].index.values

# Join the values with commas, and "and" before the last value
formatted_class_string = ', '.join(selected_indices[:-1]) + ', and ' + selected_indices[-1] if len(selected_indices) > 1 else selected_indices[0]


subset = joined_with_membership[joined_with_membership["cluster_group_name"] == 'TH-EPI-Glut']
correlation_TH_EPI = subset[exp.columns].corr()
correlation_TH_EPI = round(correlation_TH_EPI[correlation_TH_EPI<1], 3)


out = {}
for i, neigh in enumerate(['Pallium-Glut',
     'Subpallium-GABA',
     'MB-HB-CB-GABA',
     'MB-HB-Glut-Sero-Dopa',
     'HY-EA-Glut-GABA',
     'TH-EPI-Glut']):
            coloc = coloc_by_neigh[neigh]
            out[neigh] = coloc[coloc["Colocalized (%)"]<100].groupby("Gene2")["Colocalized (%)"].mean()

coloc_matrix = round(pd.DataFrame.from_dict(out), 2)

gene="Htr1a"
clu_by_expr = joined.groupby("cluster")[gene].apply(percentage_above_threshold).sort_values(ascending=False)
perc_enriched_htr1a = round((joined_boolean[joined_boolean["cluster"].isin(clu_by_expr[clu_by_expr > threshold_enriched_clusters].index)][gene].sum() /joined_boolean[gene].sum())*100, 2)
gene="Htr1b"
clu_by_expr = joined.groupby("cluster")[gene].apply(percentage_above_threshold).sort_values(ascending=False)
perc_enriched_htr1b = round((joined_boolean[joined_boolean["cluster"].isin(clu_by_expr[clu_by_expr > threshold_enriched_clusters].index)][gene].sum() /joined_boolean[gene].sum())*100, 2)
gene="Htr1d"
clu_by_expr = joined.groupby("cluster")[gene].apply(percentage_above_threshold).sort_values(ascending=False)
perc_enriched_htr1d = round((joined_boolean[joined_boolean["cluster"].isin(clu_by_expr[clu_by_expr > threshold_enriched_clusters].index)][gene].sum() /joined_boolean[gene].sum())*100, 2)
