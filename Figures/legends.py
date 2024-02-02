from pathlib import Path
import numpy as np
import pandas as pd
from Utils.Settings import output_folder_calculations

gene_filtered  = pd.read_csv(Path(output_folder_calculations, "selected_genes_RNAseq.csv"))
selected_genes_cl = np.sort(gene_filtered["gene_symbol"].values)


legends =  {"Figure 1. Overview of Htrs translation in the RNA-seq dataset.":
                "(A) Heatmap showing absolutwe number of cells expressing each Htrs. Inset shows the same information in percentage of the total. "
                "(B) UMAP representation color-coded by neighborhood metadata (left), Htr1 (middle) and Htr2 (right) expression. "
                "(C) Htr expression prevalence in cells grouped by neurotransmitter release (top). Confusion matrix of the multi-label random forest classifier showing "
                "true label on y axis and predicted label on x axis (middle). Matrix of absolute SHAP values for each group and receptor (bottom). "
                "(D) Htr expression prevalence in cells grouped byclass. "
                "(E) Htrs expression correlation matrix. "
                "(F) Htrs colocalization matrix. Each dot represents the percentage of colocalization of gene on x axis in cells expressing gene on y axis. "
                "(G) Percentage of cells expressing the gene on x axis expressing at least another Htrs gene (top) or at least other 4 Htrs (bottom). "
                "(H) Pie charts representing the proportion of principal Htrs grouped by intracellular effector for each neighborhood. Principal effector was "
                "identified by summing the expression of Htrs. Each number represents the "
                "number of cells in thousands. "}

legends.update({f"{gene} transcription": f"(A) On the left, {gene} prevalence across neighborhoods with squared Pearson correlation coefficient (R²) between RNA-seq and MERFISH dataset.. On the right, amount of {gene} RNA detected using "
                                              f"RNA-seq (top) and MERFISH (bottom). "
                                              f"(B) Amount of colocalization with each Htrs by cells expressing {gene} RNA (left). Number of Htrs RNA detected in cells "
                                              f"expressing {gene} RNA (right). "
                                              f"(C) Prevalence of {gene} RNA across all classes of cells in RNA-seq and MERFISH dataset. Inset represents the linear regression between the two datasets. "
                                              f"On te right, absolute number of cells expressing {gene} RNA by class ranked in descending order (top ten). "
                                              f"(D) Prevalence of {gene} RNA across divisions (left) and structures (right). Inset represents the proportion of cells expressing {gene} RNA that belongs to enriched clusters. "
                                              f"(E) Prevalence (top) and average amount of RNA expression across the antero-posteroir axis. "
                                              f"(F) Expression of {gene} RNA detected by MERFISH in 4 representative slices. Border color represents the position on the antero-posterior axis. "
                                              f"  " for gene in selected_genes_cl})