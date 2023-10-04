from Utils.Utils import htr_cmap, htr_cmap_rgb
from Utils.Settings import output_folder_calculations
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

coloc = pd.read_pickle(f"{output_folder_calculations}/colocalization_broad.pkl")

ax = sns.barplot(data=coloc.groupby( "Gene1")["Value"].mean().reset_index(), x="Gene1", y="Value", palette=htr_cmap)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
ax.set_xlabel("5-HT receptors")
ax.set_ylabel("Percentage co-localized (%)")
ax.tick_params(axis='x', rotation=45)
for ytick in ax.get_xticklabels():
    ytick.set_color(htr_cmap_rgb[ytick.get_text()])

plt.show()