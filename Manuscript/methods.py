from Utils.Settings import threshold_expression_MERFISH, threshold_enriched_clusters


methods = {"Jupyter notebooks structure": "The entire analysis is contained in 2 jupyter notebooks hosted on Github at https://github.com/RobertoDF/Transcriptomics-5-HT. "
                                          "Within the 'Figures' folder, 'Figure_1.ipynb' and 'Figure_2.ipynb' notebooks reproduce all figures contained in the paper. "
                                            "To adapt the code for the visualization of different genes it "
                                          "is sufficient to change the 'family_name' and 'genes_families' variables in Utils.Settings.py file. Data is downloaded "
                                          "following the instructions provided by the Allen Institute (https://alleninstitute.github.io/abc_atlas_access/intro.html). "
                                          "Notebooks to download the RNA-seq and MERFISH datasets are contained in the 'Load_Data' folder. To explore the expression of different genes, "
                                          "it is necessary to download the associated expression matrices by changing the selected genes in the 'Download_RNAseq_data.ipynb' notebook, "
                                          "this can be achieved by modifying the cells underneath the headings 'Select genes RNA-seq' and 'Select genes MERFISH'. ",

            "Data preparation":             "We loaded the metadata and the precomputed expression matrices for the RNA-seq dataset relative to all Htr genes "
                                                      "(see 'Load_data/Download_RNAseq_data.ipynb'). We also loaded the metadata relative to the "
                                                      "'cluster_group_name' residing in originally in a different .csv file ('Find membership df' in 'Figure_1.ipynb'). "
                                                      "This information is referred to as 'membership'. "
                                                      "Each of this data structure is a pandas dataframe that can be easily joined together according to the unique "
                                                      "cell label index ('joined' dataframe). A different dataframe containing membership information is created ('joined_with_membershiop'), "
                                                      "this is necessary because some cells belong to multiple 'cluster_group_name' or 'neighborhood' and therefore causes a doubling of some rows. "
                                                      "We will use the dataframe contianing 'membership information' only in to visualize information relative to 'cluster_group_name'. \n"
                                                        "The MERFISH dataset was loaded in a similar fashion (see 'Load data MERFISH' in 'Figure_2.ipynb'). "
                                            "This dataset is split in 4 different dataframes ('Zhuang-ABCA-1', "
                                            "'Zhuang-ABCA-2', 'Zhuang-ABCA-3' and 'Zhuang-ABCA-4') stored in a dictionary ('cell_expression')."
                                            " We concatenated the 4 dataframe in one data structure called 'data_merfish' "
                                            "using the '.concat()' pandas method. Additionally, we used the spatial informations of each cell belonging to the MERFISH dataset for "
                                                             "the registration to the Allen Mouse Brain Common Coordinate Framework (CCF) and , subsequently, "
                                                             "we assigned parcellations labels ('CCF registration and parcellation annotation' in 'Figure_2.ipynb'). ",

            "Overview figure visualization": "This figure relies uniquely on the RNA-seq dataset. "
                                                      "In panel A we use a heatmap to visualize both the amount of transcription per cell and the number of cells transcribing each Htr "
                                                      "contained in the dataset. In panel B we used the precomputed UMAP coordinates and plot on the color axis information about the "
                                                      "most transcribed gene per selected family. \n"
                                                      "In panel C and D we plot the percentage of cells transcribing each Htr grouped by neurotransmitter release. We take advantage of the pandas 'Group by' "
                                                      "function to concisely perform this computation: after grouping by the selected variable (in this case 'neurotransmitter') we apply a function called "
                                                      "'percentage_above_threshold' to compute the percentage of cells within a group expressing a gene above a threshold. "
                                                      "The 'percentage_above_threshold' function is defined within the "
                                                      "'Utils.Utils.py' file. "
                                                      "The threshold is stored in the 'Utils.Settings.py' file ('threshold_expression'). "
                                                      "The confusion matrix is computed within the 'decodddit' function in 'Figure_1.ipynb'. Here we use a boolean version of the 'joined' dataframe created using "
                                                      "the same threshold ('threshold_expression').  The dataset containing boolean values for gene expression (joined_boolean) was "
                                                      "filtered to include columns of interest, specifically a selector column (sel) and a list of selected genes (selected_genes)."
                                                      " The resulting DataFrame was indexed by the selector column, which represented the target variable (neurotransmitter type), "
                                                      "while the remaining columns contained features corresponding to the expression levels of various serotonin receptor genes (Htr). "
                                                      "The features for classification were defined as the expression levels of the serotonin receptor genes, and the target variable was "
                                                      "the neurotransmitter type. The dataset was divided into training and testing sets using a stratified sampling approach to ensure that "
                                                      "the class distribution of the target variable was maintained in both subsets. The test set comprised 5% of the total data. A Random "
                                                      "Forest classifier was initialized with parameters set to 10 decision trees, a maximum tree depth of 10, and 30 parallel jobs to leverage "
                                                      "computational resources effectively. The classifier was configured to handle class imbalances by adjusting the class weights. The model "
                                                      "was then trained using the training data. The trained Random Forest model was used to predict the neurotransmitter types on the test "
                                                      "dataset. The performance of the model was evaluated by comparing the predicted labels with the actual labels. The accuracy of the model "
                                                      "was calculated and transcribed as a percentage. Additionally, a comprehensive classification report was generated, providing metrics such as "
                                                      "precision, recall, and F1-score for each class. A confusion matrix, normalized by the true labels, was also produced to visualize the model's"
                                                      " classification performance across different neurotransmitter types. The data manipulation and analysis were conducted using the pandas "
                                                      "library. The machine learning model was implemented using the scikit-learn library, specifically the RandomForestClassifier for"
                                                      " classification tasks. The evaluation of the model's performance was performed using scikit-learn's accuracy score, classification report, "
                                                      "and confusion matrix functions. SHAP (SHapley Additive exPlanations) values were calculated to interpret the feature "
                                                      "importance of the Random Forest classifier. An explainer object was created using SHAP's TreeExplainer, which was specifically "
                                                      "designed for tree-based models. The explainer was initialized with the trained Random Forest classifier, and the number of parallel "
                                                      "jobs was set to 40 to leverage computational resources effectively. The SHAP values were computed for a sample of the feature set "
                                                      "(X_sample). These values indicate the contribution of each feature to the model's predictions. \n"
                                                      "In panel E we plot the correlation between transcription of different Htr genes by using the pandas 'corr()' method. \n "
                                                      "To plot the co-localization data of panel F a dictionary named 'coexp' was initialized to store the co-localization results. "
                                                      "This dictionary would eventually hold the percentage of co-localization for each pair of genes. A nested "
                                                      "loop was employed to iterate through each pair of selected genes, excluding a placeholder category labeled "
                                                      "'Any Htr'. For each target gene and gene to check, the following computations were performed: "
                                                      "Co-localization Calculation: For each gene pair, the boolean DataFrame joined_boolean was used to"
                                                      " check whether both genes were transcribed (True) in each sample. This was done using the .all(axis=1) "
                                                      "method, which returned True for rows where both genes were transcribed. The sum of these True values "
                                                      "indicated the total number of samples where both genes were co-transcribed. Normalization: This sum was "
                                                      "then normalized by dividing it by the total number of samples where the target gene was transcribed. This "
                                                      "provided the percentage of samples where the gene pair was co-transcribed relative to the expression of the "
                                                      "target gene. Storing Results: The computed co-localization percentage for each gene pair was stored in the "
                                                      "coexp dictionary with the gene pair as the key. After computing the co-localization percentages for all gene "
                                                      "pairs, the results were converted into a pandas DataFrame for further analysis and visualization. The same colocalization "
                                                      "was used in the barplots of panel G. \n"
                                                      "For panel H we aggregated Htr expression by family by aggregating the expression levels of specific serotonin "
                                                      "receptor genes. These genes were grouped into four primary families: Htr1/5: Summing the expression levels of "
                                                      "genes Htr1a, Htr1b, Htr1d, Htr1f, Htr5a, and Htr5b. Htr2: Summing the expression levels of genes Htr2a, Htr2b, "
                                                      "and Htr2c. Htr4/6/7: Summing the expression levels of genes Htr4, Htr6, and Htr7. Htr3: Summing the expression "
                                                      "levels of genes Htr3a and Htr3b. These aggregated values were combined with additional columns representing "
                                                      "neuronal classifications (class, subclass, supertype, and cluster_group_name). The columns of the resulting "
                                                      "DataFrame were labeled accordingly, and a new column ('Primary Htr family') was added. "
                                                      "This column identified the primary serotonin receptor family for each entry by determining the "
                                                      "family with the highest aggregated expression. ",

            "Receptor figure preparation and visualization":  "This figure relies on both the RNA-seq and MERFISH datasets. "
                                                              "In panel A we plot both the prevalence and the average amount of transcription of the selected gene in the two datasets. "
                                                             "We excluded from the analysis the 'NN-IMN-GC' neighborhood because of consistently low transcription across all Htr genes. For the visualization "
                                                              "of gene expression patterns across different 'neighborhoods', we utilized the Seaborn library in Python to create point plots. "
                                                              "Specifically, we employed the 'sns.pointplot' function to illustrate the expression levels of a given gene across various groups. The 'sns.violinplot' "
                                                              "function was used to plot violin plots of amount of transcription per group. \n"
                                                              "In panel B we used the same co-localization data used in Figure 1 panel F (RNA-seq dataset), This barplot is a 'sliced' version of that panel focusing on one receptor at the time. "
                                                              "To visualize the number of colocalized genes (barplot on the right), we utilized a boolean DataFrame ('joined_boolean') to"
                                                              " filter for selected genes and focus on the expression status of a particular gene. "
                                                              "We then calculated the sum of true values (indicating gene transcription) across each row "
                                                              "where the specific gene was transcribed. The distribution of these sums was normalized to obtain "
                                                              "the percentage of samples exhibiting co-expression of the genes. \n"
                                                              "In panel C on the left we repeat the same computation of panel A but using 'class' as grouping variable. On the right, "
                                                              "we plotted the raw number of cells transcribing the selected gene across different classes. "
                                                              "We first filtered the 'joined' DataFrame to include only rows where the expression level of a specific gene exceeded a defined threshold ('threshold_expression'). "
                                                              "We then counted the occurrences of each class in this filtered dataset. The top 10 classes with the highest counts were "
                                                              "selected for visualization. Using Seaborn's barplot function, we created a bar plot to display the distribution of these classes. "
                                                              "The y-axis represented the count of occurrences, while the x-axis denoted the different classes. \n"
                                                              "In panel D we plotted the prevalence of the selected gene in brain regions at two different hyerarchical levels, 'division' and 'structure'. "
                                                              "Here we take advantage of the high-confidence label integration between the 'RNA-seq' and 'MERFISH dataset' {Zhang, 2023 #2887}. "
                                                              "Each cell of the 'MERFISH' "
                                                              "dataset is assigned a cell-type label ('class', 'subclass', 'supertype' and 'cluster') from the clustering of the 'RNA-seq' {Yao, 2023 #2886}."
                                                              "To analyze the expression of specific genes across different brain regions and neuronal clusters, we utilized a multi-step data processing approach. "
                                                              "First, we calculated the percentage of cells within each cluster expressing the target gene above a defined threshold, grouping the data by cluster. "
                                                              "This allowed us to identify clusters with high gene expression levels in the RNA-seq. Next, we focused on clusters with significant gene expression, "
                                                              "filtering the 'MERFISH' dataset "
                                                              "to include only these enriched clusters. We then computed the prevalence of gene expression within these clusters across different parcellation divisions and "
                                                              "structures. This was done by normalizing the number of cells expressing the gene in each division or structure by the total number of cells in that division or "
                                                              "structure, expressed as a percentage. The results were visualized using bar plots to illustrate the top 10 parcellation divisions and structures with the highest "
                                                              "gene expression prevalence. Additionally, "
                                                              "we included an inset pie chart to show the proportion of gene expression attributable to the enriched clusters relative to the total gene expression. This "
                                                              "pie chart highlighted the contribution of these enriched clusters to the overall expression of the target gene. This pie plot shows that often only a minority of cells "
                                                              "transcribing a selected gene belong to enriched clusters. "
                                                              "It is possible to visualize the proportion of cells transcribing a gene grouped by area bypassing the 'enriched cluster' computation using the online visualizer "
                                                              "('Overview genes by brain structure', Data source selector='MERFISH'). \n"
                                                              "Panel E: We calculated the percentage of cells within each cluster that expressed the target gene above a defined threshold, "
                                                              "allowing us to identify clusters with enriched gene expression. "
                                                              "The number of such enriched clusters was printed for reference. Next, we focused on cells within these"
                                                              " enriched clusters and calculated the prevalence of the target gene's expression across different "
                                                              "brain sections. This was done by normalizing the number of cells expressing the gene in each section "
                                                              "by the total number of cells in that section, expressed as a percentage. The results were plotted using a "
                                                              "line plot to illustrate the gene's prevalence across brain sections. \n"
                                                              "Panel F: To visualize the expression of a specific gene in the top four brain sections, we implemented a function plot_4_best_sections. "
                                                              "This function aimed to identify and plot the sections with the highest gene expression levels. The steps are as follows: "
                                                                "Data Preparation: We first prepared the dataset by selecting the relevant brain sections and ensuring that unassigned "
                                                              "parcellation divisions were excluded. We merged this dataset with cluster membership information to provide context for the gene expression data."
                                                                "Gene Expression Calculation: The percentage of cells within each cluster expressing the target gene above a defined threshold was calculated.^"
                                                              "This allowed us to identify clusters with enriched gene expression."
                                                                "Section Identification: We calculated the prevalence of the target gene's expression in each brain section. "
                                                              "Using these prevalence values, we identified the top four sections with the highest gene expression. "
                                                              "Peaks in the expression data, spaced adequately apart, were determined "
                                                                 "using the find_peaks function. The top four peaks were selected for visualization. "
                                                                "Color Assignment: A specific color was assigned to each of the "
                                                              "top four sections to differentiate them in the plots. A predefined list of colors was used "
                                                                  "to ensure consistency and clarity. "
                                                                "Plotting: For each of the top four sections, the gene expression data was plotted. "
                                                              "The plot_slice function was used to generate the plots "
                                                                   "for each section, and the border color of each subplot was set to match the "
                                                                   "assigned color for the respective section. This helped in visually distinguishing each section. "
                                                                "The final figure comprised four subplots, each representing one of the top four brain sections with the highest gene expression levels, "
                                                                                                                              "providing a clear and comparative visualization of the gene expression "
                                                                                                                              "patterns across these key sections.",

           "Data Visualizer":   "The visualizer was built in Python using Matplotlib, Holoviews and Panel libraries. It is available as a jupyter notebook ('Figures/INteractive_vizs.ipynb') and online. "
                                "The jupyter notebook can be used by following the installation instructions available in https://github.com/RobertoDF/Transcriptomics-5-HT. "
                                "The visualizer is deployed and accessible online on the Hugging Face portal. "
                                "It is organized in 4 different tabs: 'Spatial MERFISH', 'Gene by class/subclass/supertype/cluster', "
                                "'Overview genes by class' and 'Overview genes by brain structure'. "
                                "The 'Spatial MERFISH' and 'Overview genes by brain structure' are associated with the MERFISH dataset, "
                                "remaining tabs are associated with the RNA-seq dataset. Each tab is associated to different interactive controls and panels. "
                                "'Spatial MERFISH': 5 interactive controls enable the selections of different datasets from {Zhang, 2023 #2887}, brain section, gene, class and subclass. "
                                "The datasets available are 2 coronal (Zhuang-ABCA-1/2) and 2 sagittals (Zhuang-ABCA-3/4). The brain section selector enables the visualization of different slices. "
                                "The gene selector enables the selection of a specific gene. Class and subclass selector restrict the visualization to selected groups. 6 panels are provided. "
                                "From top to bottom: lineplot representing the proportion of cells selected out the cells available across the spatial axis associated to each dataset, "
                                "lineplot representing the amount of transcription across space of the selected gene, lineplot representing the percentage of "
                                "cells across space in which RNA of "
                                f"the selected gene was detected (threshold set at {threshold_expression_MERFISH }), barplot representing the percentage of Htr positive cells in the selected slice "
                                f"grouped by brain structure (number in each bar is the absolute number of cells) and two panels representing the slice "
                                f"selected with the gene transcription on the left and atlas metadata on the right. \n"
                                
                                "'Gene by class/subclass/supertype/cluster': 2 interactive controls enable the selections of neighborhood group and gene. "
                                "The neighborhood selector enables the selection of a specific neighborhood. "
                                "The gene selector enables the selection of a specific gene. "
                                "For each class of neurons we provide 3 levels of visualization. On top, violinplots representing the gene prevalence by subclass; in the middle, "
                                "violinplots representing prevalence by supertype and on the bottom "
                                "barplots representing prevalence by cluster. Each subclass is color-coded according to the panel available for each class. \n"
                                
                                "'Overview genes by class': 4 interactive controls enable the selections of class, subclass, type of grouping and sorting. "
                                "The class and subclass selectors enable the selection of a specific class and subclass, respectively. The plot can begrouped at different levels of detail: classes, "
                                "subclasses, supertypes and even individual clusters (the number of groups that can visualized at the same time is limited by the maximum recursion depth of Holoviews). "
                                "The plot can be sorted by the group´s alphabetical name or gene expression. Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits "
                                'of the current selection. Y axis is populated by the name of the groups selected by the "Group by" selector. X axis shows each Htrs. \n'
                                
                                "'Overview genes by brain structure': 4 interactive controls enable the selections of dazta soirce, division, neurotransmitter and sorting. "
                                "The division and neurotransmitter selectors enable the selection of a specific brain division and neurotransmitter, respectively. "
                                "Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits "
                                f"of the current selection. Gene prevalence is limited to cluster enriched in the according gene (prevalence within cluster of the gene >{threshold_enriched_clusters}%). "
                                f"The y axis is populated by the brain structures belonging to the currently selected brain division. For each division we can "
                                "refine our selection by isolating neurons releasing a specific neurotransmitter. X axis shows each Htrs. \n"
                                "First we identify enriched clusters in the RNA-seq dataset then we look at the proportion of cells belonging to enriched clusters over the total number of cells per region. "
                                "This creates a problem in the cases where a 5-HT receptor was not deemed important by the clustering algo, "
                                "in these cases there might be a really small amount of enriched clusters and most cells will be ignored. "
                                "To solve this we create a data source selector to enable the possibility to bypass RNA-seq and look at MERRFISH data only. In this case we look directly at the ratio of cells expressing transcribing each gene over the total number of cells per region. "

           }
