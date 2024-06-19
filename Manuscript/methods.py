from Utils.Settings import threshold_expression_MERFISH, threshold_enriched_clusters


methods = {"Jupyter notebooks structure": "The entire analysis is contained in 2 jupyter notebooks hosted on Github at https://github.com/RobertoDF/Transcriptomics-5-HT. "
                                          "For data analysis and visualization we employed mainly pandas, numpy, matplotlib, scikit-learn and seaborn python libraries. "
                                          "Within the 'Figures' folder, 'Figure_1.ipynb' and 'Figure_2.ipynb' notebooks reproduce all figures contained in the paper. "
                                          "All parameters relative to the analysis are contained in Utils.Settings.py. "
                                          "Data is downloaded "
                                          "following the instructions provided by the Allen Institute (https://alleninstitute.github.io/abc_atlas_access/intro.html). "
                                          "Notebooks to download the scRNA-seq and MERFISH datasets are contained in the 'Load_Data' folder. "
                                          "To explore the transcription of different genes, "
                                          "it is necessary to download the associated transcription matrices by changing the selected genes in the 'Download_RNAseq_data.ipynb' notebook, "
                                          "this can be achieved by modifying the cells underneath the headings 'Select genes scRNA-seq' and 'Select genes MERFISH'. It is also necessary "
                                          "to change the 'family_name' and 'genes_families' variables in Utils.Settings.py file. ",

            "Data preparation":             "We loaded the metadata and the precomputed transcription matrices ('exp' pandas dataframe) for the scRNA-seq dataset relative to all Htr genes "
                                                      "(see 'Load_data/Download_RNAseq_data.ipynb'). We also loaded the metadata relative to the "
                                                      "'cluster_group_name' (or 'neighborhood' in the text) residing originally in a different .csv file ('Find membership df' in 'Figure_1.ipynb'). "
                                                      "This information is referred to as 'membership'. Additionally we loaded cell metadata information ('cell' dataframe). "
                                                      "Each of these data structures are pandas dataframes that can be easily joined together according to the unique "
                                                      "cell label index ('joined' dataframe). A different dataframe containing membership information is created ('joined_with_membership'), "
                                                      "this is necessary because some cells are assigned to multiple 'cluster_group_name' and therefore cause the duplication of some dataframe´s rows. "
                                                      "We used the dataframe containing 'membership information' only to visualize information relative to 'cluster_group_name'. \n"
                                                        "The MERFISH dataset was loaded in a similar fashion (see 'Load data MERFISH' in 'Figure_2.ipynb'). "
                                            "This dataset is split in 4 different dataframes ('Zhuang-ABCA-1', "
                                            "'Zhuang-ABCA-2', 'Zhuang-ABCA-3' and 'Zhuang-ABCA-4') stored in a dictionary ('cell_expression'). "
                                            "We concatenated the 4 dataframe in one data structure called 'data_merfish' "
                                            "using the '.concat()' pandas method. Additionally, we used the spatial information of each cell belonging to the MERFISH dataset for "
                                                             "the registration to the Allen Mouse Brain Common Coordinate Framework (CCF) and, subsequently, "
                                                             "we assigned parcellations labels ('CCF registration and parcellation annotation' in 'Figure_2.ipynb'). The dataframes loaded by both datasets "
                                            "already included all the clustering labels (class, subclass, supertype and cluster). For details about the clustering see 'Clustering scRNA-seq data' section in {Yao, 2023 #2886}. ",

            "Overview figure": "This figure relies uniquely on the scRNA-seq dataset. "
                                                      "In panel A we used a heatmap to visualize both the amount of transcription per cell and the number of cells transcribing each Htr "
                                                      "contained in the dataset using the 'exp' dataframe. "
                                                     "\nIn panel B we used the precomputed UMAP coordinates available in the 'joined' dataframe to create a scatterplot and "
                                                    "plotted on the color axis information about the "
                                                      "most transcribed gene per selected family (either Ht1 or Ht2). \n"
                                                      "In panel C we plotted the percentage of cells transcribing each Htr grouped by neurotransmitter release. "
                                                        "We take advantage of pandas 'Group by' "
                                                      "function to concisely perform this computation: after grouping by the selected variable (in this case 'neurotransmitter') we apply a function called "
                                                      "'percentage_above_threshold' to compute the percentage of cells within a group transcribing a gene above a threshold. "
                                                      "The 'percentage_above_threshold' function is defined within the "
                                                      "'Utils.Utils.py' file. "
                                                      "The threshold is stored in the 'Utils.Settings.py' file ('threshold_expression'). "
                                                      "The confusion matrix is computed within the 'decoddddddd' function defined in Utils.Utils.py. "
                                                     "This function uses a boolean version of the 'joined' dataframe created using "
                                                      "the same threshold ('threshold_expression'). The dataset containing boolean values for gene transcription ('joined_boolean') was "
                                                      "filtered to include columns of interest, specifically a selector column ('sel') and a list of selected genes ('selected_genes')."
                                                      " The resulting dataframe was indexed by the selector column, which represented the target variable, "
                                                      "while the remaining columns contained features corresponding to the transcription levels of various serotonin receptor genes (Htr). "
                                                      "In this particular case, the features for classification were defined as the boolean transcription of the various 5-HT receptor genes, and the target variable was "
                                                      "the neurotransmitter type. A Random Forest classifier ('RandomForestClassifier' from scokit-learn) was initialized with 200 estimators, a maximum depth of 10, balanced class weights, and parallel processing across 20 jobs. "
                                                        "Linear models such as 'LogisticRegression' and 'LinearDiscriminantAnalysis' were found to underperfom the Random Forest classifier "
                               "(respectevely, 0.3768 and 0.249 accuracy vs 0.385 for the Random Forest classifier, see 'Test linear models' in Figure_1.ipynb). "
                               "Using Stratified K-Fold cross-validation with 5 ('n_splits' set in Utils.Settings.py) folds, balanced accuracy scores were computed, and the mean accuracy was reported. "
                                                      "Predictions were generated with cross-validation ('cross_val_predict' function in scikit-learn). "
                               "The performance of the model was evaluated by comparing the predicted labels with the actual labels. "
                               " Additionally, a comprehensive classification report was generated, providing metrics such as "
                                                      "precision, recall, and F1-score for each class. A confusion matrix, normalized by the true labels, was also produced to visualize the model's"
                                                      " classification performance across different neurotransmitter types. The evaluation of the model's performance was performed using scikit-learn's 'balanced_accuracy_score', 'classification_report', "
                                                      "and 'confusion_matrix' functions. SHAP (SHapley Additive exPlanations) values were calculated to interpret the feature "
                                                      "importance of the Random Forest classifier. An explainer object was created using SHAP's 'TreeExplainer', which was specifically "
                                                      "designed for tree-based models. The explainer was initialized with the trained Random Forest classifier, and the number of parallel "
                                                      "jobs was set to 40 to leverage computational resources effectively. The SHAP values were computed for a sample of the feature set "
                                                        "of 10,000 observations based on class weights "
                                                      "('X_sample'). These values indicate the contribution of each feature to the model's predictions. \n"
                                                     "In panel D we plotted the percentage of cells transcribing each Htr grouped by class label, "
                                                        "additional plots related to classification accuracy were computed following the "
                                                        "instructions of the previous panel and are available as supplementary figure. "
                                                      "In panel E we plotted the correlation between transcription of different Htr genes by using the pandas 'corr()' method. \n"
                                                      "To plot the co-localization data of panel F a dictionary named 'coexp' was initialized to store the co-localization results. "
                                                      "This dictionary would eventually hold the percentage of co-localization for each pair of genes. A nested "
                                                      "loop was employed to iterate through each pair of selected genes, excluding a placeholder category labeled "
                                                      "'Any Htr'. For each target gene and gene to check, the following computations were performed: "
                                                      "Co-localization Calculation: For each gene pair, the boolean dataframe 'joined_boolean' was used to"
                                                      " check whether both genes were transcribed (True) in each sample. This was done using the '.all(axis=1)' "
                                                      "method, which returned True for rows where both genes were transcribed. The sum of these True values "
                                                      "indicated the total number of samples where both genes were co-transcribed. Normalization: This sum was "
                                                      "then normalized by dividing it by the total number of samples where the target gene was transcribed. This "
                                                      "provided the percentage of samples where the gene pair was co-transcribed relative to the transcription of the "
                                                      "target gene. Storing Results: The computed co-localization percentage for each gene pair was stored in the "
                                                      "coexp dictionary with the gene pair as the key. After computing the co-localization percentages for all gene "
                                                      "pairs, the results were converted into a pandas dataframe for further analysis and visualization. The same co-localization "
                                                      "was used in the barplots of panel G. \n"
                                                      "For panel H we aggregated Htr transcription by family. "
                                                    "These genes were grouped into four primary families: Htr1/5: Summing the transcription levels of "
                                                      "genes Htr1a, Htr1b, Htr1d, Htr1f, Htr5a, and Htr5b. Htr2: Summing the transcription levels of genes Htr2a, Htr2b, "
                                                      "and Htr2c. Htr4/6/7: Summing the transcription levels of genes Htr4, Htr6, and Htr7. Htr3: Summing the transcription "
                                                      "levels of genes Htr3a and Htr3b. These aggregated values were combined with additional columns representing "
                                                      "neuronal classifications (class, subclass, supertype, and cluster_group_name). The columns of the resulting "
                                                      "dataframe were labeled accordingly, and a new column ('Primary Htr family') was added. "
                                                      "This column identified the primary serotonin receptor family for each entry by determining the "
                                                      "family with the highest aggregated transcription. ",

            "Receptor figure":  "This figure relies on both the scRNA-seq and MERFISH datasets. "
                                                              "In panel A we plot both the prevalence and the average amount of transcription of the selected gene in the two datasets. "
                                                             "We excluded from the analysis the 'NN-IMN-GC' neighborhood because of consistently low transcription across all Htr genes. For the visualization "
                                                              "of gene transcription patterns across different 'neighborhoods', we used the seaborn 'pointplot' function to illustrate the transcription levels of a given gene across various groups. The 'violinplot' "
                                                              "function was used to create violin plots of amount of transcription per group. \n"
                                                              "In panel B we used the same co-localization data used in Figure 1 panel F (scRNA-seq dataset), This barplot is a 'sliced' version of that panel focusing on one receptor at the time. "
                                                              "To visualize the number of colocalized genes (barplot on the right), we utilized a boolean dataframe ('joined_boolean') to"
                                                              " filter for selected genes and focus on the transcription status of a particular gene. "
                                                              "We then calculated the sum of true values (indicating gene transcription) across each row "
                                                              "where the specific gene was transcribed. The distribution of these sums was normalized to obtain "
                                                              "the percentage of samples exhibiting co-transcription of the genes. \n"
                                                              "In panel C on the left we repeat the same computation of panel A but using 'class' as grouping variable. On the right, "
                                                              "we plotted the raw number of cells transcribing the selected gene across different classes. "
                                                              "We first filtered the 'joined' dataframe to include only rows where the transcription level of a specific gene exceeded a defined threshold ('threshold_expression'). "
                                                              "We then counted the occurrences of each class in this filtered dataset. The top 10 classes with the highest counts were "
                                                              "selected for visualization. Using Seaborn's barplot function, we created a bar plot to display the distribution of these classes. "
                                                              "The y-axis represented the count of occurrences, while the x-axis denoted the different classes. \n"
                                                              "In panel D we plotted the prevalence of the selected gene in brain regions at two different hierarchical levels, 'division' and 'structure'. "
                                                              "Here we take advantage of the high-confidence label integration between the scRNA-seq and MERFISH dataset {Zhang, 2023 #2887}. "
                                                              "Each cell of the MERFISH "
                                                              "dataset is assigned a cell-type label ('class', 'subclass', 'supertype' and 'cluster') from the clustering of the scRNA-seq {Yao, 2023 #2886}."
                                                              "To analyze the transcription of specific genes across different brain regions and neuronal clusters, we utilized a multi-step data processing approach. "
                                                              "First, we calculated in the scRNA-seq the percentage of cells within each cluster transcribing the target gene above a defined threshold ('threshold_expression'), grouping the data by cluster. "
                                                              "This allowed us to identify clusters with high gene transcription levels (>70%, 'threshold_enriched_clusters' in Utils.Settings.py) in the scRNA-seq. "
                                                                "Next, we focused on clusters with significant gene transcription, "
                                                              "filtering the MERFISH dataset "
                                                              "to include only cells belonging to these enriched clusters. We then computed the prevalence of cells transcribing the selected gene across different parcellation divisions and "
                                                              "structures. This was done by normalizing the number of cells transcribing the gene in each division or structure by the total number of cells in that division or "
                                                              "structure, expressed as a percentage. The results were visualized using bar plots to illustrate the top 10 parcellation divisions and structures with the highest "
                                                              "gene transcription prevalence. Additionally, "
                                                              "we included an inset pie chart to show the proportion of gene transcription attributable to the enriched clusters relative to the total gene transcription. This "
                                                              "pie chart highlighted the contribution of these enriched clusters to the overall transcription of the target gene. This pie chart showed that in some cases only a minority of cells "
                                                              "transcribing a selected gene belongs to enriched clusters, in this cases, consequetially, the majority of cells will be ignored. "
                                                              "To address this problem we include in the interactive visualizer ('Overview genes by brain structure' dashboard) "
                                                              "a data source selector that can switch the algorithm used by the dashboard from 'scRNA-seq + MERFISH'' (the one described above) "
                                                              "to 'MERFISH only'. This latter option computes the prevalence using solely the MERFISH dataset by simply calculating the proportion of cells transcibing "
                                                              "the selected genes across spatial groups. \n"
                                                              "Panel E: We calculated the percentage of cells within each cluster that expressed the target gene above a defined threshold ('threshold_expression'), "
                                                              "allowing us to identify clusters with enriched gene transcription. Next, we focused on cells within these "
                                                              "enriched clusters and calculated the prevalence of the target gene's transcription across different "
                                                              "brain sections. This was done by normalizing the number of cells transcribing the gene in each section "
                                                              "by the total number of cells in that section, expressed as a percentage. The results were plotted using a "
                                                              "line plot to illustrate the gene's prevalence across brain sections. \n"
                                                              "Panel F: To visualize the transcription of a specific gene in the top four brain sections, "
                                "we implemented a function called 'plot_4_best_sections' (in 'Figures/Figure_2.ipynb'). "
                                                              "This function aimed to identify and plot the sections with the highest gene transcription levels. "
                                                                "Data preparation: We first prepared the dataset by selecting the relevant brain sections and ensuring that unassigned "
                                                              "parcellation divisions were excluded. We merged this dataset with cluster membership information to provide context for the gene transcription data."
                                                                "Gene transcription calculation: The percentage of cells within each cluster transcribing the target gene above a defined threshold was calculated. "
                                                              "This allowed us to identify clusters with enriched gene transcription. "
                                                                "Section identification: We calculated the prevalence of the target gene's transcription in each brain section. "
                                                              "Using these prevalence values, we identified the top four sections with the highest gene transcription. "
                                                              "Peaks in the transcription data, spaced adequately apart, were determined "
                                                                 "using the 'find_peaks' function from scipy. The top four peaks were selected for visualization. "
                                                                "Plotting: For each of the top four sections, the gene transcription data was plotted. "
                                                              "The plot_slice function was used to generate the plots "
                                                                   "for each section, and the border color of each subplot was set to match the "
                                                                   "assigned color for the respective section. "
                                                                "The final figure comprised four subplots, each representing one of the top four brain sections with the highest gene transcription levels, "
                                                                                                                              "providing a clear and comparative visualization of the gene transcription "
                                                                                                                              "patterns across these key sections.",

           "Data Visualizer":   "The visualizer was built in Python using Matplotlib, Holoviews and Panel libraries. It is available as a jupyter notebook ('Figures/Interactive_vizs.ipynb') "
                                "and online (https://rdef654875678597657-5-ht-transcriptomics.hf.space). "
                                "The jupyter notebook can be used locally by following the installation instructions available in https://github.com/RobertoDF/Transcriptomics-5-HT. "
                                "The visualizer is deployed and accessible online on the Hugging Face portal. "
                                "It is organized in four different dashboards: 'Spatial MERFISH', 'Gene by class/subclass/supertype/cluster', "
                                "'Overview genes by class' and 'Overview genes by brain structure'. "
                                "The 'Spatial MERFISH' and 'Overview genes by brain structure' are associated with the MERFISH dataset, "
                                "remaining tabs are associated with the scRNA-seq dataset. Each dashboard's data source is annotated in the title.. "
                                "'Spatial MERFISH': Five interactive controls enable the selections of different datasets from {Zhang, 2023 #2887}, brain section, gene, class and subclass. "
                                "The datasets available are 2 coronal (Zhuang-ABCA-1/2) and 2 sagittal (Zhuang-ABCA-3/4). "
                                "The controls allow visualization of different slices, specific genes, and selected groups. The dashboard includes six panels: "
                                "1. Lineplot representing the proportion of cells selected across the spatial axis associated to each dataset, "
                                "2. Lineplot representing the amount of transcription across space of the selected gene, 3. Lineplot representing the percentage of "
                                "cells across space in which RNA of "
                                f"the selected gene was detected (threshold set at {threshold_expression_MERFISH }), 4. Barplot representing the percentage "
                                f"of Htr positive cells in the selected slice "
                                f"grouped by brain structure (number in each bar is the absolute number of cells), 5-6. Slice selected with "
                                f"gene transcription (left) and atlas metadata (right). \n"
                                
                                "'Gene by class/subclass/supertype/cluster': This dashboard has two interactive controls for selecting neighborhood "
                                "group and gene. For each class of neurons, three levels of visualization are provided: 1. Violinplots: Gene prevalence by subclass, "
                                "2.Violinplots: Prevalence by supertype, 3. Barplots: Prevalence by cluster. \n"
                                
                                "'Overview genes by class': This dashboard includes four interactive controls for selecting class, subclass, type of grouping, and sorting. "
                                "The plot can be grouped at different clustering depths: classes, "
                                "subclasses, supertypes and even individual clusters (the number of groups that can visualized at the same time is limited by "
                                "the maximum recursion depth of Holoviews). "
                                "The plot can be sorted by the group´s alphabetical name or gene transcription. "
                                "Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits "
                                'of the current selection. Y axis is populated by the name of the groups selected by the "Group by" selector. X axis shows each Htrs. \n'
                                
                                "'Overview genes by brain structure': This dashboard includes four interactive controls for "
                                "selecting data source, division, neurotransmitter, and sorting. "
                                "Gene prevalence is represented with a heatmap in which the colorbar is updated according to the limits "
                                f"of the current selection. Gene prevalence is limited to cluster enriched in the according gene "
                                f"(prevalence within cluster of the gene >{threshold_enriched_clusters}%). "
                                f"The y axis is populated by the brain structures belonging to the currently selected brain division. For each division we can "
                                "refine our selection by isolating neurons releasing a specific neurotransmitter. X axis shows each Htrs. \n"
                                "First, enriched clusters in the scRNA-seq dataset are identified, then the proportion of cells belonging to enriched clusters"
                                " over the total number of cells per region is analyzed. To handle cases where most cells do not belong to enriched clusters and are ignored, "
                                "a 'Data Source Selector' is used to bypass scRNA-seq data and use MERRFISH data directly."
                                "In this case we look directly at the ratio of cells transcribing each gene over the total number of cells per region. "

           }
