from pathlib import Path
import json
import requests


# gene Family breakdown

family_name = "htr"

genes_families = {
    'Htr1': ['Htr1a', 'Htr1b', 'Htr1d', 'Htr1f'],
    'Htr2': ['Htr2a', 'Htr2b', 'Htr2c'],
    'Htr3': ['Htr3a', 'Htr3b'],
    'Htr4': ['Htr4'],
    'Htr5': ['Htr5a', 'Htr5b'],
    'Htr6': ['Htr6'],
    'Htr7': ['Htr7']
}

## DIRECTORIES

root_data = Path("/alznew/roberto/Allen_Institute/")
root_github_repo = r"/alznew/roberto/Github/Transcriptomics-5-HT"
output_folder_calculations = Path(f"{root_github_repo}/Processed")
output_folder_processed_lfps = Path(f"{root_data}/Processed_lfps")
output_folder = Path(f"{root_github_repo}/Output_figures")
output_folder_supplementary = Path(f"{root_github_repo}/Output_figures/Supplementary")
manuscript_folder = Path(f"{root_github_repo}/Manuscript")
utils_folder = Path(f"{root_github_repo}/Utils")
root_visualizer_data = Path('/alznew/roberto/Github/Transcriptomics-5-HT-huggingface/Data/')

version = '20230830'
download_base = Path(f"{root_data}/abc_download_root")

manifest_path = 'releases/%s/manifest.json' % version

url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/' + manifest_path

manifest = json.loads(requests.get(url).text)

# thresholds expression

threshold_expression = 3.5 # for RNAseq
threshold_expression_MERFISH = .3
threshold_enriched_clusters = 70 # percentage cells expressing receptor to be considered an enriched cluster



## custom divisions

class_to_division = {'01 IT-ET Glut': "CTX Glut",
'02 NP-CT-L6b Glut': "CTX Glut",
'03 OB-CR Glut':"CTX Glut",
'04 DG-IMN Glut':"CTX Glut",
'05 OB-IMN GABA':"CTX Gaba",
'06 CTX-CGE GABA': "CTX Gaba",
'07 CTX-MGE GABA': "CTX Gaba",
'08 CNU-MGE GABA': "CNU Gaba",
'09 CNU-LGE GABA': "CNU Gaba",
'10 LSX GABA': "CNU Gaba",
'11 CNU-HYa GABA': "CNU Gaba",
'12 HY GABA': "HY Gaba",
'13 CNU-HYa Glut': "HY Glut",
'14 HY Glut': "HY Glut",
'15 HY Gnrh1 Glut': "HY Glut",
'16 HY MM Glut': "HY Glut",
'17 MH-LH Glut': "TH Glut",
'18 TH Glut': "TH Glut",
'19 MB Glut': "MB Glut",
'20 MB GABA': "MB Gaba",
'21 MB Dopa': "MB Dopa",
'22 MB-HB Sero': "MB Sero",
'23 P Glut': "HB Glut",
'24 MY Glut': "HB Glut",
'25 Pineal Glut': "TH Glut",
'26 P GABA': "HB Gaba",
'27 MY GABA': "HB Gaba",
'28 CB GABA': "CB Gaba",
'29 CB Glut': "CB Glut",
'30 Astro-Epen': "Non-neuronal",
'31 OPC-Oligo': "Non-neuronal",
'32 OEC': "Non-neuronal",
'33 Vascular': "Non-neuronal",
'34 Immune': "Non-neuronal"}

class_to_broad_division = {
    '01 IT-ET Glut': "CTX",
    '02 NP-CT-L6b Glut': "CTX",
    '03 OB-CR Glut': "CTX",
    '04 DG-IMN Glut': "CTX",
    '05 OB-IMN GABA': "CTX",
    '06 CTX-CGE GABA': "CTX",
    '07 CTX-MGE GABA': "CTX",
    '08 CNU-MGE GABA': "CNU",
    '09 CNU-LGE GABA': "CNU",
    '10 LSX GABA': "CNU",
    '11 CNU-HYa GABA': "CNU",
    '12 HY GABA': "IB",
    '13 CNU-HYa Glut': "IB",
    '14 HY Glut': "IB",
    '15 HY Gnrh1 Glut': "IB",
    '16 HY MM Glut': "IB",
    '17 MH-LH Glut': "IB",
    '18 TH Glut': "IB",
    '19 MB Glut': "MB",
    '20 MB GABA': "MB",
    '21 MB Dopa': "MB",
    '22 MB-HB Sero': "MB",
    '23 P Glut': "HB",
    '24 MY Glut': "HB",
    '25 Pineal Glut': "IB",
    '26 P GABA': "HB",
    '27 MY GABA': "HB",
    '28 CB GABA': "CB",
    '29 CB Glut': "CB",
    '30 Astro-Epen': "Non-neuronal",
    '31 OPC-Oligo': "Non-neuronal",
    '32 OEC': "Non-neuronal",
    '33 Vascular': "Non-neuronal",
    '34 Immune': "Non-neuronal"
}

order_division = ['Isocortex', 'OLF', 'HPF', 'CTXsp', 'STR', 'PAL', 'TH', 'HY', 'MB',
       'P',  'MY', 'CB',]
order_broad_division = ["CTX", "CNU", "MB", "HB", "IB", "CB", "Non-neuronal"]
neuron_cluster_groups_order = ['Pallium-Glut', 'Subpallium-GABA', 'MB-HB-CB-GABA','MB-HB-Glut-Sero-Dopa', 'HY-EA-Glut-GABA','TH-EPI-Glut','NN-IMN-GC']

## ANALYSIS PARAMETERS

#splits of the stratified cross-validation
n_splits = 5

# Nature style figure letters

Adapt_for_Nature_style = False

