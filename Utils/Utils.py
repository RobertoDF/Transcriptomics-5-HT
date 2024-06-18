import os
print("Current working directory:", os.getcwd())
from pathlib import Path
import matplotlib.colors as plt_colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from allensdk.core.reference_space_cache import ReferenceSpaceCache
from Utils.Settings import download_base, genes_families, class_to_broad_division, output_folder_calculations, manifest, \
    threshold_expression, threshold_expression_MERFISH
from sklearn.model_selection import cross_val_score, StratifiedKFold,  cross_val_predict
from sklearn.metrics import make_scorer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import  balanced_accuracy_score, classification_report
import shap
from sklearn.metrics import confusion_matrix


metadata = manifest['file_listing']['WMB-neighborhoods']['metadata']
rpath = metadata['cluster_group_membership']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
group_membership = pd.read_csv(file) # cluster can belong to two different groups

neuron_cluster_groups = group_membership["cluster_group_name"].unique()

colors = ['#6FADCF', '#F6AE99', '#73628A', '#4F5D75',
          '#D6EFFF', '#9B5094', '#BCD39C', '#8FC0A9']


# Create a dictionary mapping each value to a color
cluster_groups_cmap = dict(zip(neuron_cluster_groups, colors))

########

base_colors = sns.color_palette("husl", n_colors=len(genes_families))
genes_cmap = {}

for idx, (family, members) in enumerate(genes_families.items()):
    shades = sns.light_palette(base_colors[idx], n_colors=len(members) + 1)[1:]
    for receptor, shade in zip(members, shades):
        genes_cmap[receptor] = shade

genes_cmap ["Any Htr"] = "#494949"

# Convert RGB to HEX
genes_cmap_rgb = {k: plt_colors.rgb2hex(v) for k, v in genes_cmap.items()}

#########
metadata = manifest['file_listing']['Allen-CCF-2020']['metadata']
rpath = metadata['parcellation_to_parcellation_term_membership_acronym']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
parcellation_annotation = pd.read_csv(file)
parcellation_annotation.set_index('parcellation_index',inplace=True)
parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]

output_dir = output_folder_calculations
reference_space_key = os.path.join('annotation', 'ccf_2017')
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest=Path(output_dir) / 'manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1)

broad_division_color_map = {}
for area in set(class_to_broad_division.values()):
    if area == "Non-neuronal":
        broad_division_color_map[area] = [.8, .8, .8]
    else:
        broad_division_color_map[area] = np.array(tree.get_structures_by_acronym([area])[0]["rgb_triplet"]) / 255


def Naturize():
    """
    Change figure panels lettering  to lowercase
    """
    for label in plt.figure(1).texts:
        if len(label.get_text())==1:
            label.set_text(label.get_text().lower())

def Naturize_text(legends_supplementary):
    """
    Change figure references to lowercase
    """
    for k, v in legends_supplementary.items():
        v_list = list(v)
        for n, char in enumerate(v_list):
            try:
                if (char.isupper()) & (v[n-1] == "(") & (v[n+1] == ")") & (char.isalpha() is True):
                    v_list[n] = char.lower()
                elif (char.isupper()) & (v[n-1] != ".") & (v[n-2] != ".") & (char.isalpha() is True) &\
                        (v[n-1].isalpha() is False) & (v[n+1].isalpha() is False)& (v[n+1]!="²"):
                    v_list[n] = char.lower()
                elif (char.isupper()) & (v[n-1].isdigit() is True):
                    v_list[n] = char.lower()
            except:
                pass
        legends_supplementary[k] = "".join(v_list)
    return legends_supplementary



def percentage_non_zero(series):
    return (series != 0).sum() / len(series) * 100

def percentage_above_threshold(series):
    return (series > threshold_expression).sum() / len(series) * 100

def percentage_above_threshold_MER(series):
    return (series > threshold_expression_MERFISH).sum() / len(series) * 100

def optimize_df(df):
    float_cols = df.select_dtypes(include='float64').columns
    df[float_cols] = df[float_cols].astype('float16')

    # Convert object columns to category if they have a limited number of unique values
    obj_cols = df.select_dtypes(include='object').columns
    for col in obj_cols:
        if df[col].nunique() <100:  # Adjust this threshold as needed
            df[col] = df[col].astype('category')
    return df


def decoddddddd(joined_boolean, sel, selected_genes, n_splits):

    df = joined_boolean[[sel] + list(selected_genes)]

    df.set_index(sel, inplace=True)

    # Assuming 'df' is your DataFrame
    X = df  # Features (Htr expression levels)
    y = df.index  # Target (neurotransmitter type)

    # Initialize the Random Forest classifier
    # You can adjust 'n_estimators' and 'max_depth' based on your dataset and computational resources
    model = RandomForestClassifier(n_estimators=200, max_depth=10, class_weight='balanced', n_jobs=20)

    skf = StratifiedKFold(n_splits=n_splits)

    scorer = make_scorer(balanced_accuracy_score)
    scores = cross_val_score(model, X, y, cv=skf, scoring=scorer, n_jobs=n_splits)

    accuracy = np.mean(scores) * 100

    print(scores)
    print("Accuracy:", accuracy)

    y_pred = cross_val_predict(model, X, y, cv=n_splits, n_jobs=n_splits)

    report = classification_report(y, y_pred, output_dict=True)
    print("\nClassification Report:\n", classification_report(y, y_pred))

    cm = confusion_matrix(y, y_pred, normalize="true")

    cm = pd.DataFrame(cm, columns=y.astype("category").categories,
                      index=y.astype("category").categories) * 100

    X_sample = X.sample(n=10000, weights=y.map(10000 / y.value_counts()), random_state=4)

    model.fit(X, y)
    # Calculate SHAP values
    explainer = shap.TreeExplainer(model, n_jobs=40)
    shap_values = explainer.shap_values(X_sample)

    out = []
    for i in range(len(shap_values)):
        out.append(pd.Series(np.abs(shap_values[i]).mean(0), name=sorted(joined_boolean[sel].unique())[i],
                             index=X.columns))

    shap_matrix = pd.concat(out, axis=1).T

    return cm, shap_matrix, accuracy, report