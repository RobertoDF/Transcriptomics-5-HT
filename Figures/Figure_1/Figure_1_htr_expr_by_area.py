from Utils.Utils import joined_boolean, htrgenes, htr_cmap_rgb
import seaborn as sns
import matplotlib.pyplot as plt

fig,axs = plt.subplots(1,2,figsize=(15,5))
sns.heatmap(joined_boolean.groupby("broad_division")[htrgenes].sum().div(joined_boolean.groupby("broad_division").size(), axis="index") * 100,
            yticklabels=True, cbar_kws={'label': 'Percentage cells (%)'}, ax=axs[0])
axs[0].set_xlabel("5-HT receptors")
for ytick in axs[0].get_xticklabels():
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

sns.heatmap(joined_boolean.groupby("division")[htrgenes].sum().div(joined_boolean.groupby("division").size(), axis="index") * 100,
            yticklabels=True, cbar_kws={'label': 'Percentage cells (%)'}, ax=axs[1])
axs[1].set_xlabel("5-HT receptors")
for ytick in axs[1].get_xticklabels():
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

plt.show()