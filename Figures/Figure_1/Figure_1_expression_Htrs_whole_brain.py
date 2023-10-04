import seaborn as sns
import numpy as np
from Utils.Utils import htr_cmap_rgb
import pandas as pd
from Utils.Utils import exp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt


ax = sns.heatmap(pd.DataFrame(np.sort(exp.values, axis=0)[::-200], columns=exp.columns,
                              index=exp.reset_index().index[::200]).T, xticklabels=5000, cbar_kws={'label': 'log(CPM)'}, cmap="Greys")
ax.set_xlabel("Number of cells")
ax.set_ylabel("5-HT receptors")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)

colors_list = [htr_cmap_rgb[receptor] for receptor in exp.columns]

for ytick, color in zip(ax.get_yticklabels(), colors_list):
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

axins = inset_axes(ax, width="50%", height="40%", loc=1)
((exp.astype("bool").sum(axis=0)/exp.shape[0])*100).plot.bar(ax=axins, color=colors_list)
axins.set_xlabel("5-HT receptors")
axins.set_ylabel("Percentage (%)")

for xtick, color in zip(axins.get_xticklabels(), colors_list):
    xtick.set_color(htr_cmap_rgb[xtick.get_text()])

plt.show()