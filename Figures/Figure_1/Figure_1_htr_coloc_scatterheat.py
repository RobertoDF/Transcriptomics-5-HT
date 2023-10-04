from Utils.Utils import htr_cmap, htr_cmap_rgb
from Utils.Settings import output_folder_calculations
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

coloc = pd.read_pickle(f"{output_folder_calculations}/colocalization_broad.pkl")
with sns.axes_style("whitegrid"):
    g = sns.relplot(
        data=coloc[coloc["Value"]<100],
        y="Gene1", x="Gene2", hue="Value", size="Value",
        palette="vlag", hue_norm=(0, 100), edgecolor=".7",
        height=10, sizes=(50, 250), size_norm=(coloc[coloc["Value"]<100]["Value"].min(), coloc[coloc["Value"]<100]["Value"].max()),
    )
    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(.02)
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    for artist in g.legend.legendHandles:
        artist.set_edgecolor(".7")

    for ax in g.axes.flat:  # Loop through the individual axes
        for label in ax.get_xticklabels():
            # Your code to modify the x-tick labels goes here
            label.set_color(htr_cmap_rgb[label.get_text()])
        for label in ax.get_yticklabels():
            # Your code to modify the x-tick labels goes here
            label.set_color(htr_cmap_rgb[label.get_text()])

plt.show()