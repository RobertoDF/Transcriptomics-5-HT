# 🧬Transcriptomic Mapping of the 5-HT Receptor Landscape🧬

The complete manuscript (figures + text) can be reproduced using Python code present in this repository. Each number present in the text is an f-String that can be traced back to a data structure (see End2EndPaper).

## 🟧[Online visualizer](https://huggingface.co/spaces/RDeF654875678597657/5-HT-Transcriptomics) 🧠

*The 'Spatial MERFISH' and 'Overview genes by brain structure' are associated with the MERFISH dataset, remaining tabs are associated with the RNA-seq dataset. Each tab is associated to different interactive controls and panels. 
More details in the associated manuscript.*

### 'Spatial MERFISH': 
5 interactive controls enable the selections of different datasets from {Zhang, 2023 #2887}, brain section, gene, class and subclass. The datasets available are 2 coronal (Zhuang-ABCA-1/2) and 2 sagittals (Zhuang-ABCA-3/4). The brain section selector enables the visualization of different slices. The gene selector enables the selection of a specific gene. Class and subclass selector restrict the visualization to selected groups. 6 panels are provided. From top to bottom: lineplot representing the proportion of cells selected out the cells available across the spatial axis associated to each dataset, lineplot representing the amount of transcription across space of the selected gene, lineplot representing the percentage of cells across space in which RNA of the selected gene was detected (threshold set at 0.3), barplot representing the percentage of Htr positive cells in the selected slice grouped by brain structure (number in each bar is the absolute number of cells) and two panels representing the slice selected with the gene transcription on the left and atlas metadata on the right. 

### 'Gene by class/subclass/supertype/cluster': 
2 interactive controls enable the selections of neighborhood group and gene. The neighborhood selector enables the selection of a specific neighborhood. The gene selector enables the selection of a specific gene. For each class of neurons we provide 3 levels of visualization. On top, violinplots representing the gene prevalence by subclass; in the middle, violinplots representing prevalence by supertype and on the bottom barplots representing prevalence by cluster. Each subclass is color-coded according to the panel available for each class. 

### 'Overview genes by class': 
4 interactive controls enable the selections of class, subclass, type of grouping and sorting. The class and subclass selectors enable the selection of a specific class and subclass, respectively. 
The plot can be grouped at different levels of detail: classes, subclasses, supertypes and even individual clusters (the number of groups that can visualized at the same time is limited by the maximum recursion depth of Holoviews). 
The plot can be sorted by the group´s alphabetical name or gene expression. Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits of the current selection. Y axis is populated by the name of the groups selected by the "Group by" selector. X axis shows each Htrs. 

### 'Overview genes by brain structure': 
3 interactive controls enable the selections of which data source to use, division and neurotransmitter. The division and neurotransmitter selectors enable the selection of a specific brain division and neurotransmitter, respectively. Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits of the current selection. Gene prevalence is limited to cluster enriched in the according gene (prevalence within cluster of the gene >70%). The y axis is populated by the brain structures belonging to the currently selected brain division. For each division we can refine our selection by isolating neurons releasing a specific neurotransmitter. X axis shows each Htrs. 

First we identify enriched clusters in the RNA-seq dataset then we look at the proportion of cells belonging to enriched clusters over the total number of cells per region. 
This creates a problem in the cases where a 5-HT receptor was not deemed important by the clustering algo, 
in these cases there might be a really small amount of enriched clusters and most cells will be ignored. 
To solve this we create a data source selector to enable the possibility to bypass RNA-seq and look at MERRFISH data only. In this case we look directly at the ratio of cells expressing transcribing each gene over the total number of cells per region. 


## ⬛️Description

The entire analysis is contained in 2 Jupyter notebooks located in the folder 'Figures': 'Figure_1.ipynb' and 'Figure_2.ipynb'. To adapt the code for the visualization of different genes, it is sufficient to change the'family_name' and 'genes_families' variables in the Utils.Settings.py file. Data is downloaded following the [instructions](https://alleninstitute.github.io/abc_atlas_access/intro.html) provided by the Allen Institute. Notebooks to download the RNA-seq and MERFISH datasets are contained in the 'Load_Data' folder. To explore the expression of different genes, it is necessary to download the associated expression matrices by changing the selected genes in the 'Download_RNAseq_data.ipynb' notebook. This can be achieved by modifying the cells underneath the headings 'Select genes RNA-seq' and 'Select genes MERFISH'.

## ⬜️Installation
```
conda create --name Transciptomics-5-HT --file requirements.txt
```

## 🟩License

cc-by-nc-nd-4.0 🪪
