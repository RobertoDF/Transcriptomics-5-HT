from Utils.Utils import htrgenes, joined_boolean
import pandas as pd
from Utils.Settings import output_folder_calculations

coexp = {}
for gene_target in htrgenes:
    for gene_to_check in htrgenes:
        coexp[gene_target, gene_to_check] = joined_boolean[[gene_target, gene_to_check]].all(axis=1).sum()/(joined_boolean[gene_target]>0).sum()


coloc = pd.DataFrame.from_dict(coexp, orient='index', columns=['Value'])
coloc[['Gene1', 'Gene2']] = pd.DataFrame(coloc.index.tolist(), index=coloc.index)
coloc = coloc.reset_index(drop=True)

# Reorder the columns for clarity
coloc = coloc[['Gene1', 'Gene2', 'Value']]
coloc["Value"] = coloc["Value"] * 100

coloc.to_pickle(f"{output_folder_calculations}/colocalization_broad.pkl")

coexp = {}
for division in joined_boolean["division"].unique():
    for gene_target in htrgenes:
        for gene_to_check in htrgenes:
            coexp[gene_target, gene_to_check, division] = joined_boolean[joined_boolean["division"]==division][[gene_target, gene_to_check]].all(axis=1).sum()/(joined_boolean[joined_boolean["division"]==division][gene_target]>0).sum()

coloc = pd.DataFrame.from_dict(coexp, orient='index', columns=['Value'])
coloc[['Gene1', 'Gene2', "Division"]] = pd.DataFrame(coloc.index.tolist(), index=coloc.index)
coloc = coloc.reset_index(drop=True)

# Reorder the columns for clarity
coloc = coloc[["Division",'Gene1', 'Gene2', 'Value']]

coloc.to_pickle(f"{output_folder_calculations}/colocalization.pkl")