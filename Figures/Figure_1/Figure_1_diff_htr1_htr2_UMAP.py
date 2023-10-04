import numpy as np
from Utils.Utils import broad_division_color_map
from Utils.Settings import htr_families
import matplotlib.patches as mpatches
from Utils.Utils import subsampled
import matplotlib.pyplot as plt

fig,axs = plt.subplots(1,3, figsize=(15,5))
_ = subsampled.dropna()

legend_handles = []

# Create a Patch (a colored box) for each name and color, and add it to legend_handles
for name, color in broad_division_color_map.items():
    legend_handles.append(mpatches.Patch(color=color, label=name))

axs[0].legend(handles=legend_handles, title='Broad Division', ncol=int(np.ceil(len(legend_handles)/2)), loc='upper left', bbox_to_anchor=(0, 1.12),
             handlelength=.5, handletextpad=0.1, borderaxespad=0, framealpha=.2)

axs[0].scatter(_['x'], _['y'], c=_["broad_division_color"], s=0.5, marker='.')
axs[0].axis('off')

for gene in htr_families["Htr1"]:
   axs[1].scatter(_['x'], _['y'], c=_[gene], s=0.5, marker='.', cmap="Reds")
axs[1].axis('off')
axs[1].set_title('Htr1 family')

for gene in htr_families["Htr2"]:
   axs[2].scatter(_['x'], _['y'], c=_[gene], s=0.5, marker='.', cmap='YlOrBr')

axs[2].axis('off')
axs[2].set_title('Htr2 family')

plt.show()