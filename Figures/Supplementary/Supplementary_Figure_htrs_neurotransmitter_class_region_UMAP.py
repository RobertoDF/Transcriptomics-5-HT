import numpy as np
from Utils.Utils import htrgenes
from Utils.Settings import htr_families
import matplotlib.patches as mpatches
from Utils.Utils import subsampled
import matplotlib.pyplot as plt



k = len(htrgenes) + 6

# Calculate number of rows for the grid
ncols = 3
nrows = int(np.ceil(k / ncols))

# Create a figure with the desired overall width and height
fig, axs = plt.subplots(nrows,ncols,figsize=(25,45))

for n, item in enumerate(["neurotransmitter", "class"]):
    _ = subsampled[[item, item + "_color"]].drop_duplicates()
    _ = _.sort_values(item)
    # Create custom legend elements
    legend_elements = [mpatches.Patch(color=row[item + '_color'], label=row[item]) for x, row in _.iterrows()]
    axs[0,n].legend(handles=legend_elements, ncol=3, loc="center")

_ = subsampled[["region_of_interest_acronym", "region_of_interest_color"]].drop_duplicates()
_ = _.sort_values("region_of_interest_acronym")
# Create custom legend elements
legend_elements = [mpatches.Patch(color=row["region_of_interest_color"], label=row["region_of_interest_acronym"]) for x, row in _.iterrows()]
axs[0,2].legend(handles=legend_elements, ncol=3, loc="center")

axs[1,0].scatter(subsampled['x'], subsampled['y'], c=subsampled['neurotransmitter_color'], s=0.5, marker='.')
axs[1,0].set_title("Neurotransmitter")
axs[1,2].scatter(subsampled['x'], subsampled['y'], c=subsampled['region_of_interest_color'], s=0.5, marker='.')
axs[1,2].set_title("Region")
axs[1,1].scatter(subsampled['x'], subsampled['y'], c=subsampled['class_color'], s=0.5, marker='.')
axs[1,1].set_title("Class")

for n, gene in enumerate(htrgenes):
    row = n // ncols
    col = n % ncols
    sc = axs[row + 2, col ].scatter(subsampled['x'], subsampled['y'], c=subsampled[gene], s=0.5, marker='.', cmap=plt.cm.magma_r)
    axs[row + 2, col].set_title(gene)

for ax_row in axs:
    for ax in ax_row:
        ax.axis('off')

plt.tight_layout()

plt.show()