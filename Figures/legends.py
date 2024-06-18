from pathlib import Path
import numpy as np
import pandas as pd
from Utils.Settings import output_folder_calculations, threshold_expression

gene_filtered  = pd.read_csv(Path(output_folder_calculations, "selected_genes_RNAseq.csv"))
selected_genes_for_main_fig = ['Htr1a', 'Htr1b', 'Htr1f', 'Htr2a', 'Htr2c',
      'Htr4', 'Htr7',  'Htr3a',]


legends =  {"Figure 1. Overview of Htrs translation.":
                "(A) Barplot showing absolute number of cells transcribing each Htrs, amount of expression is represented in greyscale, no threshold is applied. "
                f"Inset shows the prevalence of each Htr using  a threshold (log(CPM)>{threshold_expression}). "
                "(B) UMAP representation color-coded by neighborhood metadata (left), Htr1 (middle) and Htr2 (right) transcription. "
                "(C) Htr expression prevalence in cells grouped by neurotransmitter release (top). Confusion matrix of the multi-label random forest classifier showing "
                "true label on y axis and predicted label on x axis (middle). Matrix of absolute SHAP values for each group and receptor (bottom). Glut = Glutamate, GABA = Gamma-Aminobutyric Acid, Glut-GABA = Glutamate and Gamma-Aminobutyric Acid, Dopa = Dopamine, None = No specific neurotransmitter, GABA-Glyc = Gamma-Aminobutyric Acid and Glycine, Chol = Acetylcholine (Cholinergic), Hist = Histamine, Sero = Serotonin, Nora = Norepinephrine. "
                "(D) Htr expression prevalence in cells grouped byclass. "
                "(E) Htrs expression correlation matrix. "
                "(F) Htrs colocalization matrix. Each dot represents the percentage of colocalization of gene on x axis in cells transcribing gene on y axis. "
                "(G) Top: Percentage of cells transcribing the number of Htrs on the x axis. Percentage of cells transcribing the gene on x axis transcribing at least another Htrs gene (middle) or at least other 4 Htrs (bottom). "
                "(H) Pie charts representing the proportion of principal Htrs grouped by intracellular effector for each neighborhood. Principal effector was "
                "identified by summing the expression of Htrs. Each number represents the "
                "number of cells in thousands. "}

legends.update({f"Figure {n+2}. {gene} transcription": f"(A) On the left, dotplot representing {gene} prevalence across neighborhoods with squared Pearson correlation coefficient (R²) between RNA-seq "
                                                       f"and MERFISH dataset. On the right, violinplots representing the amount of {gene} RNA detected using "
                                              f"RNA-seq (top) and MERFISH (bottom). "
                                              f"(B) Amount of colocalization with each Htrs by cells expressing {gene} RNA in the scRNAseq dataset (left). Number of Htrs RNA detected in cells "
                                              f"expressing {gene} RNA in the scRNAseq dataset  (right). "
                                              f"(C) Prevalence of {gene} RNA across all classes of cells in RNA-seq and MERFISH dataset. Inset represents the linear regression between the two datasets. "
                                              f"On te right, absolute number of cells expressing {gene} RNA in the scRNAseq by class, ranked in descending order (top ten). "
                                              f"(D) Ranked prevalence of {gene} RNA across divisions (left) and structures in the MERFISH dataset(right). "
                                                       f"Inset represents the proportion of cells expressing {gene} RNA that belongs to enriched clusters. "
                                              f"(E) Prevalence (top) and average amount of RNA expression in the MERFISH dataset across the antero-posterior axis of cells belonging to enriched clusters. "
                                              f"(F) Expression of {gene} RNA detected by MERFISH in 4 representative slices. Border color represents the position on the antero-posterior axis. "
                                              f"  " for n, gene in enumerate(selected_genes_for_main_fig)})