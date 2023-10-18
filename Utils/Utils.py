import os
import pandas as pd
from Utils.Settings import  download_base, url, htr_families, class_to_broad_division, class_to_division, output_folder_calculations, manifest
from pathlib import Path
import matplotlib.colors as plt_colors
import seaborn as sns
import json
import requests
import numpy as np
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import os
import matplotlib.pyplot as plt
import param
import numpy as np
import holoviews as hv
from holoviews.util.transform import dim
from holoviews.plotting.bokeh.element import ColorbarPlot, LegendPlot
from holoviews.plotting.bokeh.styles import line_properties, fill_properties

from holoviews.core.dimension import Dimension
from holoviews.core.data import Dataset

# cmaps

metadata = manifest['file_listing']['WMB-10X']['metadata']

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

base_colors = sns.color_palette("husl", n_colors=len(htr_families))
htr_cmap = {}

for idx, (family, members) in enumerate(htr_families.items()):
    shades = sns.light_palette(base_colors[idx], n_colors=len(members) + 1)[1:]
    for receptor, shade in zip(members, shades):
        htr_cmap[receptor] = shade

# Convert RGB to HEX
htr_cmap_rgb = {k: plt_colors.rgb2hex(v) for k, v in htr_cmap.items()}

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
                        (v[n-1].isalpha() is False) & (v[n+1].isalpha() is False)& (v[n+1]!="�"):
                    v_list[n] = char.lower()
                elif (char.isupper()) & (v[n-1].isdigit() is True):
                    v_list[n] = char.lower()
            except:
                pass
        legends_supplementary[k] = "".join(v_list)
    return legends_supplementary



def percentage_non_zero(series):
    return (series != 0).sum() / len(series) * 100


