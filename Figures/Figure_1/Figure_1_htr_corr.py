from Utils.Utils import exp, htr_cmap_rgb
import seaborn as sns
import matplotlib.pyplot as plt

colors_list = [htr_cmap_rgb[receptor] for receptor in exp.columns]

_ = exp.corr()
ax = sns.heatmap(_[_<1], vmax=.5)
for ytick in ax.get_xticklabels():
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])
for ytick, color in zip(ax.get_yticklabels(), colors_list):
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

plt.show()