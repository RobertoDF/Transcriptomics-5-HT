import pylustrator
import matplotlib.pyplot as plt
from Utils.Settings import output_folder, output_folder_figures_calculations, Adapt_for_Nature_style
from Utils.Utils import Naturize
from Utils.Settings import root_github_repo

pylustrator.start()

pylustrator.load(f"{root_github_repo}/Figures/Figure_1/Figure_1_expression_Htrs_whole_brain.py", offset=[0, 0])
pylustrator.load(f"{root_github_repo}/Figures/Figure_1/Figure_1_diff_htr1_htr2_UMAP.py", offset=[1, 0])
pylustrator.load(f"{root_github_repo}/Figures/Figure_1/Figure_1_htr_expr_by_neurotransmitter.py", offset=[0, 1])

#if Adapt_for_Nature_style is True:
#    Naturize()

plt.show()
#plt.savefig(f"{output_folder}/Figure_1", dpi=300)