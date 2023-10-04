from Utils.Utils import joined_boolean, htrgenes, htr_cmap_rgb
import seaborn as sns
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5,10))
sns.heatmap(joined_boolean.groupby("class")[htrgenes].sum().div(joined_boolean.groupby("class").size(), axis="index") * 100,
            yticklabels=True, cbar_kws={'label': 'Percentage cells (%)'}, ax=ax)
ax.set_xlabel("5-HT receptors")
for ytick in ax.get_xticklabels():
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

plt.show()