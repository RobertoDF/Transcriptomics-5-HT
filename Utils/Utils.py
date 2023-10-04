import os
import pandas as pd
from Utils.Settings import  download_base, url, htr_families, class_to_broad_division, class_to_division
from pathlib import Path
import matplotlib.colors as plt_colors
import seaborn as sns
import json
import requests
import numpy as np
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import os
import matplotlib.pyplot as plt

manifest = json.loads(requests.get(url).text)

output_dir = download_base
reference_space_key = os.path.join('annotation', 'ccf_2017')
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest=Path(output_dir) / 'manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1)


## color maps

base_colors = sns.color_palette("husl", n_colors=len(htr_families))
htr_cmap = {}

for idx, (family, members) in enumerate(htr_families.items()):
    shades = sns.light_palette(base_colors[idx], n_colors=len(members) + 1)[1:]
    for receptor, shade in zip(members, shades):
        htr_cmap[receptor] = shade

htr_cmap_rgb = {k: plt_colors.rgb2hex(v) for k, v in htr_cmap.items()}



broad_division_color_map = {}
for area in set(class_to_broad_division.values()):
    if area == "Non-neuronal":
        broad_division_color_map[area] = [.8,.8,.8]
    else:
        broad_division_color_map[area] = np.array(tree.get_structures_by_acronym([area])[0]["rgb_triplet"])/255



metadata = manifest['file_listing']['WMB-10X']['metadata']
rpath = metadata['example_genes_all_cells_expression']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
exp = pd.read_csv(file)
exp.set_index('cell_label',inplace=True)
exp = exp.sort_index(axis=1)

rpath = metadata['cell_metadata_with_cluster_annotation']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath )
cell = pd.read_csv(file, keep_default_na=False)
cell.set_index('cell_label',inplace=True)

cell["division"] = cell['class'].map(class_to_division)
cell["broad_division"] = cell['class'].map(class_to_broad_division)

cell['broad_division_color'] = cell['broad_division'].map(broad_division_color_map)

joined = cell.join(exp)
subsampled = joined.loc[::30]

joined_boolean =  cell.join( exp.astype("bool") )

htrgenes = exp.columns

htrgenes = htrgenes.sort_values()





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
                        (v[n-1].isalpha() is False) & (v[n+1].isalpha() is False)& (v[n+1]!="�"):
                    v_list[n] = char.lower()
                elif (char.isupper()) & (v[n-1].isdigit() is True):
                    v_list[n] = char.lower()
            except:
                pass
        legends_supplementary[k] = "".join(v_list)
    return legends_supplementary







