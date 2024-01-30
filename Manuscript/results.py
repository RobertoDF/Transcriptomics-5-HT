from Utils.Results_variables import *
from Utils.Settings import Adapt_for_Nature_style
from Utils.Utils import Naturize_text

results = {"Transcriptomic overview of 5-HT receptors landscape":
               "We analysed the single-cell scRNA-seq dataset provided by the Allen Institute {Yao, 2023 #2828} "
                "focusing on the expression of Htrs RNA across approximately 4 million brain cells. The scRNA-seq dataset comprehensively encompassed all known 14 Htr subtypes. "
               f"{round(((exp>threshold_expression).sum(axis=1).astype(bool).sum()/exp.shape[0])*100,2)}% of cells expressed RNA of at least one Htr. "
               f"Prevalence of Htrs across the entire dataset was considerably different ranging from {expression_total.min()}% of {expression_total.idxmin()} "
               f"to {expression_total.max()}% of {expression_total.idxmax()} (Figure 1A). RNA of 6 Htr was found in less than 2.5% of the cells "
               f"({', '.join(expression_total[expression_total < 2.5].index.values)}). On the other hand, RNA of Htr1f, Htr2a and Htr2c was present in "
               f"at least 1 every 5 cells. Average amount of RNA expression also varied across receptors (Supplementary Figure 1A). "
               f"Interestingly, the variation in amount of RNA shared around half of the variability with the prevalence, i.e., receptors found in more cells also tended "
               f"to be expressed more at the single cell dimension. "
               f"Beside the prevalence and amount of expression, also the distribution was considerably different. This is exemplified "
               f"by looking at the distribution of the Htr1 and Htr2 families across different cell type neighborhoods (Figure 1B). "
               f"Neighborhoods were defined both by location and neurotransmitter (Supplementary Figure 1B, Table 1). "
               f"When looking at the UMAP distributions for individual Htr, considerable differences were present also within family (Supplementary Figure 2). "
               f"We analyzed this differences more in detail by grouping cells by "
               f"neurotransmitter, neighborhoods or class. These categorizarions divided cells in a higlhy skewed manner (Supplementary Figure 1C), "
               f"for example when looking at groups by neurotransmitter release, "
               f"3 groups (Glut, Gaba and unassigned) made up for almost the totality of cells ({round((joined.groupby('neurotransmitter').size() / joined.shape[0])[['','Glut','GABA']].sum()*100,2)}%). "
               f"Expectedly, the vast majority of cells was classified as excitatory (Glut) ({total_cells_by_neurotransmitter['Glut']}%) and "
               f"around 1 every 5 cells was found to release GABA ({total_cells_by_neurotransmitter['GABA']}%). All the other neurotransmitter were found in less than 1% of the cells, "
               f"in particular, 5-HT releasing neurons (Sero) were found in only {total_cells_by_neurotransmitter['Sero']}% of the cells. "
               f"Pattern of Htrs expression across different neurotransmitter groups exhibited a relatevely high Pearson correlation coefficient "
               f"(r={round(corr_by_neurotransmitter.stack().mean(),2)}±{round(corr_by_neurotransmitter.stack().sem(),2)}). "
               f"Sero and cholinergic neurons (Chol) showed the most distinct patterns of expression with respectively r={round(corr_by_neurotransmitter.mean()['Sero'],2)}±"
               f"{round(corr_by_neurotransmitter.sem()['Sero'],2)} and "
               f"r={round(corr_by_neurotransmitter.mean()['Chol'], 2)}{round(corr_by_neurotransmitter.sem()['Chol'],2)}. "
               f"To evaluate the uniqueness of Htrs RNA expression per group we employed a Random Forest Classifier aiming at decoding the grouping variable from the Htrs expression (Figure 1C). "
               f"Overall accuracy of the model was {round(accuracy_neurotransmitter,2)}%. "
               f"Reflecting the correlation analysis, Sero and Chol were among the groups with higher true positive (TP) rate (Sero={round(cm_neurotransmitter['Sero']['Sero'], 2)}%, "
               f"Chol={round(cm_neurotransmitter['Chol']['Chol'], 2)}%). Cells not expressing any Htr were also "
               f"identified succesfully (no-Htr={round(cm_neurotransmitter[''][''], 2)}%). Noradrenaline (Nora) and glycine (GABA-Glyc) releasing neurons "
               f"were also identified almost half of the times (Nora={round(cm_neurotransmitter['Nora']['Nora'], 2)}% and "
               f"GABA-Glyc={round(cm_neurotransmitter['GABA-Glyc']['GABA-Glyc'], 2)}%). "
               f"To understand the contribution of each Htr in each prediction "
               f"we calculated the mean absolute SHAP values for each receptor and neurotransmitter. The shap values associated with the mean prevalence enable us "
               f"to understand which are the defining features of each grop. Here we can see, for example, that the identification of Sero neurons "
               f"is determined mainly by expression of Htr1a and Chol neurons by Htr4 and Htr5b. Absence of expression can also contribute to the classification, e.g., "
               f"Htr4 is rarely expressed by Nora neurons. "
               f"When looking at different neighborhoods the accuracy of the model was {round(accuracy_neighborood,2)}%. "
               f"The model could differentiate best the NN-IMN-GC, TH-EPI-Glut and Pallium-Glut groups (NN-IMN-GC={round(cm_neighborood['NN-IMN-GC']['NN-IMN-GC'], 2)}%, "
               f"TH-EPI-Glut={round(cm_neighborood['TH-EPI-Glut']['TH-EPI-Glut'], 2)}% and Pallium-Glut={round(cm_neighborood['Pallium-Glut']['Pallium-Glut'], 2)}%, Supplementary Figure 3A). "
               f"NN-IMN-GC includes all the cells not releasing any neurotransmitter, their classification is thereofre predictably strongly influenced by absence of any Htr. On the other hand, "
               f"TH-EPI-Glut cells were characterized by high expression of Htr7 and low expression of Htr2a and Htr4. "
               f"Across classes, differences in Htrs expression were even more evident (Figure 1D). {report_class.loc['recall'][report_class.loc['recall']>.4].shape[0]} groups could be identified "
               f"with a TP rate > 40%: {formatted_class_string} (Supplementary Figure 3B). 04 DG-IMN Glut were charachterized by high expression of Htr4 and absence of, the usually prevalent, Htr2c. Similarly, 05 OB-IMN GABA cells "
               f"showed virtual absence of Htr2c but expression of Htr1f; 09 CNU-LGE GABA cells show high Htr1b and low Htr7; 17 MH-LH Glut exhibited high levels of Htr5b and Htr4"
               f"; 18 TH Glut showed high levels of Htr7 and virtual absence of Htr4; 22 MB-HB Sero, miroring the results showed by Sero neurons, were charachterized by Htr1a; at last, 34 Immune cells "
               f"were identified by absence of any Htr expression. "
               f"Correlation between Htrs expression across the totality of cells ranged from {corr_by_cell.min().values[0]} ({'-'.join(corr_by_cell.idxmin())}) to {corr_by_cell.max().values[0]} "
               f"({'-'.join(corr_by_cell.idxmax())}). "
               f"Considerable correlation was also found for the {corr_by_cell.index[1]} (r={corr_by_cell.iloc[1].values[0]}) and "
               f"{corr_by_cell.index[2]} (r={corr_by_cell.iloc[2].values[0]}) pairs (Figure 1E). Interestingly, correlation patterns were not stable across neighboroods "
               f"(Supplementary Figure 4A). For example, Pallium-Glut exhibited a negative correlation between Htr4-Htr2a not evident from the entire dataset. Of note, TH-EPI-Glut showed the "
               f"highest absolute correlation across all neighboroods with r={correlation_TH_EPI.max().max()} between Htr5b-Htr4 and a unique negative correlation between Htr5b-Htr7 and Htr4-Htr7. "
               "To investigate the underlying cause of the correlations we investigated colocalization between Htrs using the same stringent threshold used by the original authors of the dataset {Yao, 2023 #2828}. "
               f"Across the entire dataset we observed that the most expressed genes, Htr1f and Htr2c, were often colocalized with other genes (Figure 1G). This was a driving factor for correlation. "
               f"After dividing the dataset across neighboroods we noticed that Pallium-Glut and TH-EPI-Glut were the most peculiar groups. they showed the lowest Htr2c colocalization, i.e., "
               f"cells expressing any other Htr colocalized with Htr2c <20%, in all other groups Htr2c colocalization was at least 40%. Moreover, Pallium-Glut was the only group showing "
               f" a uniquely low Htr7 colocalization ({coloc_matrix.loc['Htr7', 'Pallium-Glut']}%) on the other hand, TH-EPI-Glut showed the highest "
               f"({coloc_matrix.loc['Htr7', 'TH-EPI-Glut']}%).  "
               f"Only rarely a cell was found to express uniquely one Htr, {round(at_least_2_receptors.mean().values[0], 2)}±{round(at_least_2_receptors.sem().values[0], 2)}% "
               f"of cells indeed expressed at least "
               f"2 Htrs (Figure 1F). Surprisingly, {round(at_least_5_receptors.mean().values[0], 2)}±{round(at_least_5_receptors.sem().values[0], 2)}% of cells expressed at least 5 Htrs. "
               f"The extensive expression across different Htr classes and the considerable coexpression within cells point at the complexity of the 5-HT system even "
               f"at the single cell dimension. "
               f"To facilitate an understanding of the downstream  cellular effects "
               f"of 5-HT we aggregated receptors according to their main intracellular effector. We aggregated Htr1 and Htr5 due to their inhibitory effect (cAMP decrerase),"
               f" Htr4, Htr6 and Htr7 because "
               f"of the shared downstream effect of increasing cAMP. Htr2 is the only one that causes an Ca2+ increase while Htr3 is the only ionotropic receptor. For each cell we determined the "
               f"the principal pathway activated by 5-HT by looking at the amount of RNA of each Htr. Afterwards we divided the cells across the different neighborhoods (Figure 1H). "
               f"Ht3 were present only in a small minority of cortical inhibitory neurons."
               f"In the telencephalon, the absolute majority of both Pallium-Glut "
               f"and Subpallium-Gaba cells were linked to Htr1/5, around one quarter of cells instead featured Htr2 as primary effector. Subcortical cells exhibited "
               f"a more balanced division without any absolute majority. "
               f" In the following sections we will take a deeper look at "
               "Htrs grouped by effector and we will take advantage of the spatial information provided by the MERFISH dataset of {Zhang, 2023 #2887} regarding 9 Htrs. ",

           "Htr1 & Htr5": u"Receptors belonging to these two families have an inhibitory effect on the host cell, they are coupled to G\u1D62 and cause a downstream decrease of cAMP "
                          "and activation of GIRK channels {Sharp, 2020 #2888; McCorvy, 2015 #2889}. "
                          "Hr1a RNA have a stable prevalence of around 10% in the brain in the RNA-seq dataset (excluding non-neuronal cells and immature neurons), with virtual absence in the TH-EPI-Glut group (Figure 2A). "
                          "Htr1a co-localized most frequently with Htr1f, Htr2c and Htr2a (Figure 2B) and only in a minority of cases was expressed alone (<10%). "
                          "Expression across classes was highly correlated between the RNA-seq and MERFISH datasets (Figure 2A) and show an almost perfect proportional relationship. "
                          "Highest expression was found in 5-HT neurons of the mid- and hindbrain (class 22 MB-HB Sero, Figure 2C), nonetheless, cortical excitatory neurons (01 IT-ET Glut) "
                          "had the higher absolute number of cells expressing Ht1a. To pinpoint the spatial location we first identified the clusters highly enriched with Htr1a RNA with "
                          f"a threshold of {threshold_enriched_clusters}%, i.e., to be classified as enriched at least {threshold_enriched_clusters}% of cells must express the receptor. "
                          f"Only {perc_enriched_htr1a}% of Htr1a expressing cells were contained in enriched clusters, pointing at a relatively low importance in the clustering algorithm. "
                          "Looking at the spatial distribution across divisions, the highest prevalence was found in cortical region of the pallidum (PAL) and hippocampus (HPF) (Figure 2D). "
                          "At a more granular level, "
                          "the highest prevalence was observed in the dorsal raphe (DR). DR expression is reflection of the high prevalence in Sero neurons outlined above, "
                          "DR contains a substantial proportion of Sero neurons. "
                          "The hippocampal structure exhibiting the higher prevalence was medial entorhinal cortex (ENTm) while the medial septum nucleus (MS) and the diagonal band nucleus (NDB), "
                          "two structures linked to generation of theta waves {Winson, 1978 #2908}, contributed substantially to the expression in PAL. "
                          "Levels of of expression were stable across the anterior-posterior axis (Figure 2E-F). "
                          
                          "Htr1b exhibited a more diverse pattern of expression across neighboroods (Figure 3A) ranging from 10 to 30%. Highest prevalence was observed in the MB-HB-Glut-Sero-Dopa group, "
                          "glutamatergic, serotonergic and dopaminergic neurons located in midbrain and hindbrain. "
                          "Colocalization showed a similar pattern compared to Htr1a (Figure 3B) and also here only a minority of cells expressed Htr1b alone (<10%)."
                          f"Looking at expression across classes, {joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).idxmax()} "
                          f"class showed the highest prevalence ({round(joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).max(), 2)}%) closely followed by "
                          f"{joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).nlargest(2).index[1]} "
                          f"({round(joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).nlargest(2).iloc[1], 2)}%) (Figure 3C). "
                          f"High expression in {joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).idxmax()} was in sharp contrast with Htr1a that showed "
                          f"only minimal expression in this class ({round(joined.groupby('class')['Htr1a'].apply(percentage_above_threshold)['09 CNU-LGE GABA'], 2)}%). "
                          f"Also in this case 01 IT-ET Glut exhibited the highest absolute number of Htr1b expessing cells."
                          f"{perc_enriched_htr1b}% of Htr1b expressing cells belonged to highly enriched clusters and"
                          f"striatum (STR) showed by far the highest prevalence with "
                          f">60% (Figure 3D-E-F). Surprisongly, Caudoputamen (CP) showed a prevalence of >40%. Nucleus accumbens (ACB) olfactory tubercle (OT), "
                          f"lateral septal nucleus (LSc) and the parabigeminal nucleus (PBG) all exhibited a prevalence of >50%. "
                          
                          f"Htr1d was expressed at a much lower level, never exceeding 7% prevalence in any neighborhood (Supplementary Figure 5A). "
                          f"It colocalized at highest levels with "
                          f"Htr2c, Htr1f and Htr1b (Supplementary Figure 5B) and only rarely was expressed alone (<5%). "
                          f"Similarly to Htr1b, expression was highest in 09 CNU-LGE GABA and 22 MB-HB Sero (Supplementary Figure 5C). Notably, 09 CNU-LGE GABA exhibited the highest absolute number of cells surpassing  "
                          f"01 IT-ET Glut. Only a small minority of cells belonged to enriched clusters ({perc_enriched_htr1d}%). "
                          f"The paraventricular nucleus of the thalamus (PT and PVT) showed the highest prevalence of around 10% (Supplementary Figure 5D-E-F), however, the triangular nucleaus of septum (TRS) located in the pallidum (PAL) division, was "
                          f"the structure with highest prevalence. "
                         
                          f"Htr1f showed the highest expression of all 5-HT receptors in the RNA-seq dataset. Higher prevalence is found in the Pallium and Subpallium groups (Figure 4A),"
                          f" reaching almost 50%. Other groups showed a prevalence of 30-40% with TH-EPI-Glut at around 20% (Figure 4A). "
                          f"Htr1f was found to colocalize the most with Htr2a and Htr2c (Figure 4B). Such high prevalence caused, however, lower levels of colocalization (Figure 4B). "
                          f"Notably, the slope of the linear regression between values provided by RNA-seq and MERFISH was "
                          f"significantly lower (Figure 4C). The two datasets are still highly correlated, half of the variability is shared. "
                          f"Spatial distribution showed a peculiarly asymettric pattern where expression was concentrated in the most anterior regions. Highest expression was observed in OLF, "
                          f"reaching over 20% consistently (Figure 4D-E-F). Specifically, highest expression was observed in the main (MOB) and accessory (AOB) olfactory bulbs. "
                          f"Both Htr5a and Htr5b were not included in the MERFISH dataset, therefore we do not have any spatial information regarding these two receptors. "
                          f"Htr5a was expressed at 10-15% prevalence across all groups with the exception of NN-IMN-GC (Supplementary Figure #) "
                          f"and colocalized the most with Htr1f, Htr2c and Htr2a. "
                          f"Prevalence a across did not show any clear peaks. "
                          f"Htr5b was expressed at a much lower level (Supplementary Figure #). "
                          f"Interestingly two classes accounted for the majority of the expression: 17 MH-LH Glut and 22 MB-HB Sero. ",

           "Htr2": "The Htr2 family is mainly linked to  Gq/11 and causes excitation by increasing intracellular Ca2+. Htr2a, famous for being instrumental in mediating the effects of psychedelics {Nichols, 2016 #854},"
                   " is found across the brain with highest prevalence in cortical groups, Pallium-Glut and Subpallium-GABA (Figure 6A). "
                   "Colocalization was highest with Htr1f and Htr2c (Figure 6B). Considerable expression (around 40%) was found in 01 IT-ET Glut, 07 CTX-MGE GABA and 16 HY-MM Glut classes (Figure 6C). "
                   "Similarly to Htr1f, also here the MERFISH dataset hinted at a lower overall expression when compared to RNA-seq. Shared variability between the two was, nonetheless very high. "
                   "Isocortex and CTXsp showed the highest prevalence, reaching more than 8% (Figure 6D). At a more detailed level, surprisingly, two structures belonging to the mammillary complex "
                   "(dorsal premammillary nucleus, PMd and tuberomammillary nucleus,TMd) "
                   "exhibited the highest prevalence. Claustrum (CLA) also showed high prevalence. Some subclasses in IT-ET Glut exhibited a particularly high prevalence, 001 CLA-EPd-CTX Car3 Glut and "
                          "027 L6b EPd Glut both had a prevalence of more than 90% (Supplementary Figure #). "
                   "Htr2b was found only in a minority of neurons and was not included in the MERFISH dataset. Interestigly, neurons belonging to "
                          f"the Pineal Glut class showed the highest prevalence at {round(joined.groupby('class')['Htr2b'].apply(percentage_above_threshold)['25 Pineal Glut'], 2)}% (Supplementary Figure #). "
                   f"Htr2c was found at highest prevalence in the MB-HB-Glut-Sero-Dopa and Hy-EA-Glut-Gaba groups (Figure 7A). Interestingly, groups with higher prevalence showed also higher levels of expression. "
                   f"Colocalization was highest with Htr1f, Htr4 and Htr7. In similar fashion to Htr2a, also here there were discrepancies between the RNA-seq and MERFISH methods (Figure 7C). "
                   f"Expression was broadly distributed across many different classes without a clear peak. The same applied to the spatial distribution, with no clear peaks (Figure 7D-E-F). "
                   f"HIghest prevalence was found in the MB, CTXsp and HY.",
           "Htr4, Htr6 and Htr7": "These receptor are all connected to Gs {McCorvy, 2015 #2889}, leading to increasing cellular levels of cAMP and excitation. "
                                  "Htr4, similarly to htr2C, showed highest prevalence (>40%) in the MB-HB-Glut-Sero-Dopa and Hy-EA-Glut-Gaba groups (Figure 8A). It colocalized the most with Htr2c and Htr1f (Figure 8B). "
                                  "Discrepancies in amount of expression between RNA-seq and MERFISH were present also here (Figure 8C). This did not affect notably, however, the correlation between the two datasets. "
                                  "Expression across classes did not exhibited any peculiar pattern. Spatial distribution,however, was more interesting, exhibiting a high prevalence in one specific structure of the STR, "
                                  "the olfactory tubercle (OT). "
                                  "We do not have spatial information about the rarely expressed Htr6. This receptor seemed to be expressed at significant prevalence only in the 09 CNU-LGE class () Supplementry Figure #. "
                                  "On the other habd, Htr7 was expressed in many more cells. It reached 60 % in the TH-EPI Glut group, and considerable amounts (around 40%) in MB, HB and HY groups (Figure 9A). "
                                  "Colocalization was the highest with Htr2c and Htr1f (Figure 9B). Expression was broadly distributed across classes present in HY, MB and TH (Figure 9C). This was reflected in the MERFISH dataset, "
                                  "showing highest prevalence in HY and TH (Figure 9D). At a structure level, the parafascicular nucleus of TH (PF) showed the highest prevalence (>40%). "
                                  "Notable expression was present in the most anterior parts of the brain (Figure 9E-F), areas belonging to OLF. "
                                  ""

           }

if Adapt_for_Nature_style is True:
    results = Naturize_text(results)
#%%
