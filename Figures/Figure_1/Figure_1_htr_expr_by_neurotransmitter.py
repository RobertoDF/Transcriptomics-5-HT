from Utils.Utils import joined_boolean, htrgenes, htr_cmap_rgb
import seaborn as sns
import matplotlib.pyplot as plt

ax = sns.heatmap(joined_boolean.groupby("neurotransmitter")[htrgenes].sum().div(joined_boolean.groupby("neurotransmitter").size(), axis="index") * 100,
            yticklabels=True, cbar_kws={'label': 'Percentage cells (%)'}, cmap="rocket")
ax.set_xlabel("5-HT receptors")
for ytick in ax.get_xticklabels():
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

plt.show()