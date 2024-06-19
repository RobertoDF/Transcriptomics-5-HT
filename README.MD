# 🧬Transcriptomic Mapping of the 5-HT Receptor Landscape🧬

The complete manuscript (figures + text) can be reproduced using Python code present in this repository. Each number present in the text is an f-String that can be traced back to a data structure (see [End2EndPaper](https://github.com/RobertoDF/End2EndPaper)).

## 🟧[Online visualizer](https://rdef654875678597657-5-ht-transcriptomics.hf.space) 🧠

*The 'Spatial MERFISH' and 'Overview genes by brain structure' are associated with the MERFISH dataset, remaining tabs are associated with the RNA-seq dataset. 
Data source for each dashboard is annotated in the title between parentheses*

**All colorbars represent prevalence (% of cells transcribing the selected gene)**

### 'Spatial MERFISH'
Five interactive controls enable the selections of different datasets from (Zhang et al., 2023), brain section, gene, class and subclass. The datasets available are 2 coronal (Zhuang-ABCA-1/2) and 2 sagittal (Zhuang-ABCA-3/4). The controls allow visualization of different slices, specific genes, and selected groups. The dashboard includes six panels: 1. Lineplot representing the proportion of cells selected across the spatial axis associated to each dataset, 2. Lineplot representing the amount of transcription across space of the selected gene, 3. Lineplot representing the percentage of cells across space in which RNA of the selected gene was detected (threshold set at 0.3), 4. Barplot representing the percentage of Htr positive cells in the selected slice grouped by brain structure (number in each bar is the absolute number of cells), 5-6. Slice selected with gene transcription (left) and atlas metadata (right). 
###  'Gene by class/subclass/supertype/cluster'
This dashboard has two interactive controls for selecting neighborhood group and gene. For each class of neurons, three levels of visualization are provided: 1. Violinplots: Gene prevalence by subclass, 2.Violinplots: Prevalence by supertype, 3. Barplots: Prevalence by cluster. 
### 'Overview genes by class/subclass'
 This dashboard includes four interactive controls for selecting class, subclass, type of grouping, and sorting. The plot can be grouped at different clustering depths: classes, subclasses, supertypes and even individual clusters (the number of groups that can visualized at the same time is limited by the maximum recursion depth of Holoviews). The plot can be sorted by the group�s alphabetical name or gene transcription. Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits of the current selection. Y axis is populated by the name of the groups selected by the "Group by" selector. X axis shows each Htrs. 
### 'Overview genes by brain structure'  
This dashboard includes four interactive controls for selecting data source, division, neurotransmitter, and sorting. Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits of the current selection. Gene prevalence is limited to cluster enriched in the according gene (prevalence within cluster of the gene >70%). The y axis is populated by the brain structures belonging to the currently selected brain division. For each division we can refine our selection by isolating neurons releasing a specific neurotransmitter. X axis shows each Htrs. 
First, enriched clusters in the scRNA-seq dataset are identified, then the proportion of cells belonging to enriched clusters over the total number of cells per region is analyzed. To handle cases where most cells do not belong to enriched clusters and are ignored, a 'Data Source Selector' is used to bypass scRNA-seq data and use MERRFISH data directly.In this case we look directly at the ratio of cells transcribing each gene over the total number of cells per region.  

## ⬛️Description

The entire analysis is contained in 2 Jupyter notebooks located in the folder 'Figures': 'Figure_1.ipynb' and 'Figure_2.ipynb'. To adapt the code for the visualization of different genes, it is sufficient to change the'family_name' and 'genes_families' variables in the Utils.Settings.py file. Data is downloaded following the [instructions](https://alleninstitute.github.io/abc_atlas_access/intro.html) provided by the Allen Institute. Notebooks to download the RNA-seq and MERFISH datasets are contained in the 'Load_Data' folder. To explore the expression of different genes, it is necessary to download the associated expression matrices by changing the selected genes in the 'Download_RNAseq_data.ipynb' notebook. This can be achieved by modifying the cells underneath the headings 'Select genes RNA-seq' and 'Select genes MERFISH'.

## ⬜️Installation
```
conda create --name Transciptomics-5-HT --file requirements.txt
```
To use the interactive visualizer locally by running the Figures/Interactive_vizs.ipynb notebook 
you need to download [this](https://huggingface.co/spaces/RDeF654875678597657/5-HT-Transcriptomics/tree/main/Data) data folder and change the 'root_visualizer_data' variable in Utils.Settings 
to reflect its position on your machine.

## 🟩License

cc-by-nc-nd-4.0 🪪
