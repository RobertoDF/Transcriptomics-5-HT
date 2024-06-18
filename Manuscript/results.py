from Utils.Results_variables import *
from Utils.Settings import Adapt_for_Nature_style
from Utils.Utils import Naturize_text

results = {"Htrs transcription overview":
               "We analyzed the single-cell scRNA-seq dataset provided by the Allen Institute {Yao, 2023 #2886} "
                "focusing on the transcription of Htrs genes across approximately 4 million brain cells passing quality control. "
               "The scRNA-seq dataset comprehensively encompassed all known 14 Htr subtypes. "
               f"{round(((exp>threshold_expression).sum(axis=1).astype(bool).sum()/exp.shape[0])*100,2)}% "
               f"of cells transcribed RNA of at least one Htr. To evaluate transcription "
               f"we used the same stringent threshold (log(CPM)>{threshold_expression}) "
               "used by the original authors to determine neurotransmitter release {Yao, 2023 #2886}. "
               f"Prevalence of Htrs, the percentage of cells transcribing a receptor, "
               f"across the entire dataset was considerably different ranging from {expression_total.min()}% of {expression_total.idxmin()} "
               f"to {expression_total.max()}% of {expression_total.idxmax()} (Figure 1A). RNA of 6 Htr was found in less than 2.5% of the cells "
               f"({', '.join(expression_total[expression_total < 2.5].index.values)}). On the other hand, RNA of Htr1f, Htr2a and Htr2c was present in "
               f"at least 1 every 5 cells. Average amount of RNA transcription also varied across receptors (Figure S1A). "
               f"Interestingly, the variation in amount of RNA shared around half (R²={r_squared_prevalence_CPM}) of the variability with the prevalence, "
               f"i.e., genes that were more widespread across cells also exhibited higher transcription rates within individual cells. "
               f"In addition to differences in prevalence and transcription levels, "
               f"the distribution of genes across the brain also showed notable variation. "
               f"This variation is highlighted by comparing the distribution patterns of the Htr1 and Htr2 families, "
               f"as showcased through uniform manifold approximation and projection (UMAP) analysis (Figure 1B). "
               f"The UMAP visualization was color-coded according to neighborhood classification. "
               "Neighborhoods, characterized by cursory anatomical proximity and molecular signatures such "
               "as neurotransmitter-release {Yao, 2023 #2886}, offer a condensed categorization of cell types (Figure S1B, Table 1). "
               f"When looking at the UMAP distributions of individual Htr, considerable differences were also present within each family of receptors (Figure S2). "
               f"We analyzed these differences grouping cells by "
               "neurotransmitter, neighborhood or class (Figure S1B). The cells were subdivided into four nested levels of classification (as defined in {Yao, 2023 #2886}) with "
               "34 classes, 338 subclasses, 1,201 supertypes and 5,322 clusters. "
               "These categorizations divided cells in a highly skewed fashion (Figure S1C). "
               f"For example, when looking at neurotransmitter-release, "
               f"3 groups (Glut, Gaba and unassigned) made up almost the totality of cells "
               f"({round((joined.groupby('neurotransmitter').size() / joined.shape[0])[['','Glut','GABA']].sum()*100,2)}%). "
               f"Expectedly, the vast majority of cells was classified as excitatory (Glut, {total_cells_by_neurotransmitter['Glut']}%) and "
               f"around 1 every 5 cells was found to release GABA ({total_cells_by_neurotransmitter['GABA']}%). "
               f"All the other neurotransmitters were found in less than 1% of the cells, "
               f"in particular, 5-HT releasing neurons (Sero) were found in only {total_cells_by_neurotransmitter['Sero']}% of the cells. "
               f"Pattern of Htrs transcription across different neurotransmitter groups exhibited a relatively high mean Pearson correlation coefficient "
               f"(r={round(corr_by_neurotransmitter.stack().mean(),2)}±{round(corr_by_neurotransmitter.stack().sem(),2)}). "
               f"Sero and cholinergic neurons (Chol) showed the most distinct patterns of "
               f"transcription with respectively mean r={round(corr_by_neurotransmitter.mean()['Sero'],2)}±"
               f"{round(corr_by_neurotransmitter.sem()['Sero'],2)} and "
               f"{round(corr_by_neurotransmitter.mean()['Chol'], 2)}±{round(corr_by_neurotransmitter.sem()['Chol'],2)} (Figure 1C). "
               f"To better evaluate the uniqueness of Htrs RNA transcription per group, and account for differences in amplitude, "
               f"not captured by simple correlation, we employed a Random Forest "
               f"Classifier aiming at decoding the grouping variable solely from Htrs' transcription. "
               f"Overall accuracy of the model in decoding neurotransmitter was {round(accuracy_neurotransmitter,2)}%"
               f" (chance level={1/joined_boolean['neurotransmitter'].nunique()*100}%). "
               f"Reflecting the correlation analysis, the confusion matrix showed that Sero and Chol were among the groups with "
               f"higher true positive (TP) rate (Sero={round(cm_neurotransmitter['Sero']['Sero'], 2)}%, "
               f"Chol={round(cm_neurotransmitter['Chol']['Chol'], 2)}%). Cells not transcribing any neurotransmitter, not exhibiting a low r beforehand, were, "
               f"nonetheless, identified even more successfully ({round(cm_neurotransmitter[''][''], 2)}%). Moreover, "
               f"Noradrenaline (Nora) and glycine (GABA-Glyc) releasing neurons "
               f"were identified at considerable levels (Nora={round(cm_neurotransmitter['Nora']['Nora'], 2)}% and "
               f"GABA-Glyc={round(cm_neurotransmitter['GABA-Glyc']['GABA-Glyc'], 2)}%). "
               f"To understand the contribution of each Htr in each prediction "
               "we calculated the mean absolute SHAP (SHapley Additive exPlanations) values for each receptor and "
               "neurotransmitter {Lundberg, 2017 #2921; Lundberg, 2020 #2922}. "
               "The SHAP values in association with the mean prevalence enabled us "
               f"to easily understand the defining features of each group. We can appreciate, for example, that the identification of Sero neurons "
               f"is determined mainly by transcription of Htr1a and Chol neurons by Htr4 and Htr5b. Crucially, "
               f"absence of transcription can also contribute to the classification, e.g., "
               f"cells not transcribing any neurotransmitter were identified mainly by absence of any Htr, and Nora neurons detection was guided by the unique absence of Htr4. "
               f"When looking at different neighborhoods the accuracy of the model was {round(accuracy_neighborhood,2)}% (chance level={1/group_membership['cluster_group_name'].nunique()*100}%). "
               f"The model could differentiate best the NN-IMN-GC, TH-EPI-Glut and Pallium-Glut groups (NN-IMN-GC={round(cm_neighborhood['NN-IMN-GC']['NN-IMN-GC'], 2)}%, "
               f"TH-EPI-Glut={round(cm_neighborhood['TH-EPI-Glut']['TH-EPI-Glut'], 2)}% and Pallium-Glut={round(cm_neighborhood['Pallium-Glut']['Pallium-Glut'], 2)}%, "
               f"Figure S3A). "
               f"NN-IMN-GC includes all the cells not releasing any neurotransmitter, their classification was "
               f"therefore expectedly influenced by absence of any Htr. On the other hand, "
               f"TH-EPI-Glut cells were characterized by the unique combination of high transcription of Htr7 and low transcription of Htr2a and Htr4, "
               f"Pallium-Glut cells, instead, exhibited relatively low levels of Htr2c and Htr7. Notably, Htr7 and Htr1f seemed to follow opposite gradients across neighborhoods. "
               f"Across classes, differences in Htrs transcription were even more striking (Figure 1D). "
               f"{report_class.loc['recall'][report_class.loc['recall']>.4].shape[0]} groups could be identified "
               f"with a TP rate >40%: {formatted_class_string} (Figure S3B). 04 DG-IMN Glut were characterized by high transcription of Htr4 and absence of "
               f"the usually prevalent Htr2c. Similarly, 05 OB-IMN GABA cells "
               f"showed virtual absence of Htr2c as well as low Htr4 and high Htr1f transcription; "
               f"09 CNU-LGE GABA cells showed high Htr1b and low Htr7/Htr1a; 17 MH-LH Glut exhibited high levels of Htr5b and Htr4"
               f"; 18 TH Glut showed high levels of Htr7 and virtual absence of Htr4; 22 MB-HB Sero, mirroring the results showed by Sero neurons,"
               f" were characterized by high levels of Htr1a; at last, 34 Immune cells "
               f"were identified by absence of any Htr transcription. The exclusive use of Htrs transcription pattern reached an "
               f"impressive {round(accuracy_class,2)}% accuracy in decoding classes (chance level={1/joined_boolean['class'].nunique()*100}%). \n"
               f"Correlation between Htrs transcription across the totality of cells ranged from "
               f"{corr_by_cell.min().values[0]} ({'-'.join(corr_by_cell.idxmin())}) to {corr_by_cell.max().values[0]} "
               f"({'-'.join(corr_by_cell.idxmax())}). "
               f"Considerable correlation was also found for the {corr_by_cell.index[1]} (r={corr_by_cell.iloc[1].values[0]}) and "
               f"{corr_by_cell.index[2]} (r={corr_by_cell.iloc[2].values[0]}) pairs (Figure 1E). Interestingly, correlation patterns were not stable across neighborhoods "
               f"(Figure S4A). For example, Pallium-Glut exhibited a unique negative correlation between Htr4-Htr2a not visible "
               f"from the analysis of the entire dataset. Of note, TH-EPI-Glut showed the "
               f"highest absolute correlation across all neighborhoods with r={correlation_TH_EPI.max().max()} "
               f"between Htr5b-Htr4 and a unique negative correlation between Htr4-Htr7. "
               "To explore the underlying causes of the correlations we analyzed colocalization (co-transcription) between Htrs. "
               f"Across the entire dataset we observed that the most transcribed genes, Htr1f and Htr2c, "
               f"were regularly transcribed whenever the RNA of any other Htr was detected (Figure 1F). "
               f"This was a driving factor for correlation. Looking more in detail across neighborhoods, also here we noticed important differences, mainly explainable "
               f"by differential prevalence of Htrs in each neighborhood. "
               f"{round(at_least_2_receptors.mean().values[0], 2)}±{round(at_least_2_receptors.sem().values[0], 2)}% "
               f"of Htr-transcribing cells exhibited at least "
               f"2 Htrs , therefore, only in a minority of cases a cell was found to transcribe uniquely one Htr "
               f"({round(((joined_boolean['Number of Htrs'].value_counts() / joined_boolean.shape[0]) * 100)[1], 2)}% of the totality of cells, Figure 1G). "
               f'Surprisingly, {round(((joined_boolean[joined_boolean["Number of Htrs"]>0]["Number of Htrs"].value_counts()/joined_boolean[joined_boolean["Number of Htrs"]>0].shape[0])*100).loc[5:].sum(), 2)}% '
               f"of Htr-transcribing cells transcribed at least 5 Htrs. "
               f"The extensive transcription of different Htr families within the same cell points at the complexity of the 5-HT system even "
               f"at the single cell dimension. \n"
               f"To facilitate an understanding of the downstream cellular effects "
               f"of 5-HT, we aggregated receptors according to their main intracellular effector. We aggregated Htr1 and Htr5 due to their inhibitory effect (cAMP decrease);"
               f" Htr4, Htr6 and Htr7 because "
               f"of the shared downstream effect of increasing cAMP; Htr2 is the only one that causes an Ca2+ increase while Htr3 is the only ionotropic receptor. "
               f"For each cell we determined "
               f"the principal pathway activated by 5-HT by analyzing the detected RNA levels for each Htr, grouping them by intracellular effector and selecting the top-ranked. "
               f"We grouped the results by neighborhood, informed by the differential Htrs' transcription (Figure 1H). "
               f"Ht3 were present only in a small minority of subpallium inhibitory neurons. "
               f"In the telencephalon, the absolute majority of both Pallium-Glut "
               f"and Subpallium-Gaba cells were linked to Htr1/5, and around one quarter of cells featured Htr2 as primary effector. Subcortical cells exhibited "
               f"a more balanced partition without any absolute majority and a considerable presence of Htr4/6/7. "
               f"In the following sections we will take a deeper look at "
               "Htrs grouped by intracellular effector, We will take advantage of the information provided by the MERFISH dataset of {Zhang, 2023 #2887} regarding 9 Htrs "
               "to analyze in detail their spatial distribution. ",

           "Htr1 & Htr5": u"Receptors belonging to these two families have an inhibitory effect on the host cell, they are coupled to G\u1D62 and cause a downstream decrease of cAMP "
                          "and activation of GIRK channels {Sharp, 2020 #2888; McCorvy, 2015 #2889}. \n"
                          "Some Htr1a agonists are currently used as anxiolytics {Parks, 1998 #2950}  and antidepressant. Htr1b and Htr1d agonists, like triptans, "
                          "are effective in treating migraines by causing vasoconstriction of cranial blood vessels."
                          "Htr1a RNA have a stable prevalence of ≈10% across neighborhoods in the scRNA-seq dataset, with virtual absence in the TH-EPI-Glut group (Figure 2A). "
                          "Htr1a co-localized most frequently with Htr1f, Htr2c and Htr2a (Figure 2B) and only in a minority of cases was transcribed alone (<10%). "
                          "Transcription across classes was highly correlated between the scRNA-seq and MERFISH datasets (Figure 2A) "
                          "and showed a good correspondence in absolute values, this was the case for the majority of others Htrs. "
                          "Highest transcription was found in Sero neurons of the mid- and hindbrain (class 22 MB-HB Sero, Figure 2C), "
                          "nonetheless, cortical excitatory neurons (01 IT-ET Glut), "
                          "like in the majority of Htrs, "
                          "contained the highest absolute number of cells transcribing the receptor. "
                          "Subclasses located in the hippocampus (HPF, see Table 2 for a list of acronyms) "
                          "contained most of the cortical cells transcribing Htr1a (see interactive visualizer, 'Overview genes by class'). "
                          "To pinpoint the spatial location, we first identified in the scRNA-seq dataset the clusters highly enriched with Htr1a RNA with "
                          f"a threshold of {threshold_enriched_clusters}%, i.e., to be classified as enriched at least {threshold_enriched_clusters}% of cells in "
                          "a cluster must express the receptor RNA. Taking advantage of the clustering label integration between the scRNAseq and MERFISH dataset (see {Zhang, 2023 #2887}), "
                          "we could identify the spatial distribution of cells belonging to enriched clusters defined using the scRNAseq. "
                          f"Only {perc_enriched_htr1a}% of Htr1a transcribing cells were contained in enriched clusters, "
                          "pointing at a relatively low importance of this receptor in the clustering algorithm used by {Yao, 2023 #2886}. "
                          "Looking at the spatial distribution across divisions, the highest prevalence was found in the pallidum (PAL) and HPF (Figure 2D). "
                          "At a more granular level, 5 of the top 10 structures by prevalence belonged to the raphe nuclei: dorsal nucleus raphe (DR), "
                          "nucleus raphe obscurus (RO), nucleus raphe pallidus (RPA), "
                          "nucleus raphe magnus (RM) and superior central nucleus raphe (CS). The high levels of Htr1a transcription in the raphe nuclei is "
                          "reflection of the high prevalence in Sero neurons outlined beforehand, "
                          "the raphe nuclei contain the vast majority of Sero neurons of the brain. "
                          "The hippocampal structure exhibiting the higher prevalence were the medial entorhinal cortex (ENTm) and the area prostata (APr) while "
                          "the medial septum nucleus (MS) and the diagonal band nucleus (NDB), "
                          "two structures linked to generation of theta waves {Winson, 1978 #2908} and containing Chol neurons, contributed substantially to the transcription in PAL. "
                          "Notably, "
                          "all these results confirms previous reports of Htr1a expression in the raphe {Haj-Dahmane, 1991 #2924;Sprouse, 1987 #2923}, ENTm {Schmitz, 1995 #2925; de Filippo, 2021 #1086} "
                          "and MS {Kia, 1996 #2926}. "
                          "Levels of of transcription were stable across the anterior-posterior axis like in most other Htrs (Figure 2E-F). "
                          "We offer the option to bypass the scRNA-seq enriched cluster calculations and directly view the prevalence "
                          "of all cells transcribing the selected gene in the MERFISH dataset using the interactive visualizer (see 'Spatial MERFISH' and 'Overview genes by brain structure' dasboards). \n"
                          
                          "Htr1b is involved in social memory persistance in mouse {Wu, 2021 #2945}. Htr1b exhibited a more diverse pattern of transcription across neighborhoods (Figure 3A) ranging from 10 to 30%. "
                          "Highest prevalence was observed in the MB-HB-Glut-Sero-Dopa group, i.e., "
                          "glutamatergic, serotonergic and dopaminergic neurons located in midbrain and hindbrain. "
                          "Colocalization showed a similar pattern compared to Htr1a (Figure 3B), only a minority of cells transcribed Htr1b alone (<10%). "
                          f"Looking at transcription across classes, the {joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).idxmax()} "
                          f"class showed the highest prevalence ({round(joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).max(), 2)}%) closely followed by "
                          f"{joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).nlargest(2).index[1]} "
                          f"({round(joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).nlargest(2).iloc[1], 2)}%) (Figure 3C). "
                          f"High transcription in {joined.groupby('class')['Htr1b'].apply(percentage_above_threshold).idxmax()} was in sharp contrast with Htr1a that showed "
                          f"only minimal transcription in this class ({round(joined.groupby('class')['Htr1a'].apply(percentage_above_threshold)['09 CNU-LGE GABA'], 2)}%). "
                          f"Also in this case, 01 IT-ET Glut exhibited the highest absolute number of Htr1b transcribing cells, specifically, in a subclass "
                          f"of the nucleus of the lateral olfactory tract "
                          f"(NLOT, see interactive visualizer, 'Overview genes by class'). "
                          f"{perc_enriched_htr1b}% of Htr1b transcribing cells belonged to highly enriched clusters and the "
                          f"striatum (STR) showed an impressive high prevalence with "
                          ">30% (Figure 3D), in line with previous reports {Maroteaux, 1992 #2928;Pommer, 2021 #2927}. "
                          "Caudoputamen (CP), Nucleus accumbens (ACB), olfactory tubercle (OT), "
                          f"lateral septal nucleus (LSc) and the parabigeminal nucleus (PBG) all exhibited a prevalence of >20%. DR and RO of the raphe nuclei also exhibited "
                          f"considerable prevalence. Similarly to Htr1a, Htr1b seems to be specifically linked to Sero neurons, in line with this, they have been reported to "
                          "mediate self-inhibition in these neurons {Hjorth, 1991 #2932}. Distribution across the antero-posterior axes reflected the "
                          f"high prevalence in STR (Figure 3E-F). \n"
                          
                          f"Htr1d was transcribed at a much lower level, never exceeding 7% prevalence in any neighborhood (Figure S5A). "
                          f"It colocalized at highest levels with "
                          f"Htr2c and Htr1f (Figure S5B) and only rarely was transcribed alone (<5%). "
                          f"Similarly to Htr1b, transcription was highest in 09 CNU-LGE GABA and 22 MB-HB Sero (Figure S5C). "
                          f"Notably, 09 CNU-LGE GABA exhibited the highest absolute number of cells surpassing "
                          f"01 IT-ET Glut. Only a small minority of Htr1d transcribing cells belonged to enriched clusters ({perc_enriched_htr1d}%). "
                          f"The paraventricular nucleus of the thalamus (PT and PVT) showed the highest prevalence at only >4% (Figure S5D-E-F). \n"
                         
                          f"Htr1f, surprisingly, showed the highest levels of transcription of all Htrs in the scRNA-seq dataset. "
                          f"Highest prevalence was found in the Pallium and Subpallium groups (Figure 4A), "
                          f"reaching ≈50%. Other groups showed a prevalence of 30-40% with TH-EPI-Glut at ≈20% (Figure 4A). "
                          f"Htr1f was found to colocalize the most with Htr2a and Htr2c (Figure 4B). In 30% of cases Htr1f was the only Htr transcribed in a cell and colocalization "
                          f"decreased linearly with the number of co-transcribed Htrs (Figure 4B). "
                          f"Notably, the slope of the linear regression between values provided by scRNA-seq and MERFISH was "
                          f"significantly lower pointing at a difference in absolute prevalence per class (Figure 4C). "
                          f"This differences are imputed to the different technique employed (see https://community.brain-map.org/t/consistent-difference-in-expression-between-zhuang-and-zeng-merfish-datasets/2604/2). "
                          f"The two datasets are, however, still highly correlated, with 66% of shared variability. "
                          f"This was the case also for Htr2a, Htr2c and Htr4. Htr1f was broadly transcribed across almost all classes, "
                          f"including some non-neuronal cells. Pineal gland cells were a notable exception. "
                          f"In absolute numbers, cortical glutamatergic cells showed the highest transcription. "
                          f"Various subclasses located in L5, claustrum (CLA) and HPF exhibited prevalence >50% "
                          f"(see interactive visualizer, 'Overview genes by class'). "
                          f"Spatial distribution showed a peculiarly asymmetric pattern with transcription concentrated in the most anterior regions. "
                          f"Highest transcription was observed in STR, olfactory areas (OLF) and the cortical subplate (CTXsp) "
                          f"reaching >20% (Figure 4D). Specifically, the highest transcription was observed in nucleus accumbens (ACB) and olfactory tract (OT), similarly to Htr1b. "
                          f"The accessory olfactory bulb (AOB) "
                          f"was the OLF structure with the highest prevalence. CLA and the endopiriform nucleus (EPd), "
                          f"on the other hand, were the CTXsp structure exhibiting the highest prevalence. "
                          f"Interestingly, in the CTXsp, transcription in Glut and Gaba neurons was anticorrelated. High prevalence in Glut neurons "
                          f"corresponded to lower prevalence in Gaba "
                          f"and vice versa. In CLA and EPd HTR1f was transcribed mainly in Glut neurons, while in the amygdala (LA, BLA, BMA) predominantly in Gaba neurons "
                          f"(see interactive visualizer, 'Overview genes by brain structure'). "
                          f"Isocortex and HPF "
                          f"also exhibited considerable transcription both in excitatory and inhibitory neurons. The amount of RNA "
                          f"transcription per cell was not linear, with a clear peak in the frontal olfactory areas (Figure 4E-F). High transcription of Htr1f in this region "
                          "was previously observed using immunohistochemistry {Bruinvels, 1994 #2929}. The broad transcription of Htr1f observed in the scRNA-seq dataset "
                          "across the entire telencephalon is in line with earlier reports {Vila-Pueyo, 2018 #2933}. \n"
                          
                          f"Both Htr5a and Htr5b were not included in the MERFISH dataset, therefore we do not have any direct spatial visualization of their transcription. "
                          f"Htr5a was transcribed at 8-16% prevalence across all neighborhoods (Figure S6A) "
                          f"and colocalized the most with Htr1f, Htr2c and Htr2a (Figure S6B). Transcription was broadly distributed across many classes, "
                          f"although only at lower levels compared to other Htrs (Figure S6C). Only one cluster was considered enriched with Htr5a in the entire "
                          f"scRNA-seq dataset, 3453 PAG-PPN Pax5 Sox21 Gaba. This cluster was located mainly in the midbrain reticular nucleus (RR, Figure S6D-E). "
                          f"Htr5b was transcribed at a much lower level across neighborhoods (Figure S7A), with a maximum of ≈%5 in TH-EPI-Glut. "
                          f"Surprisingly, even if its overall prevalence was much lower than Htr5a, 10 clusters were found to be enriched in Htr5b. "
                          f"This receptor was transcribed at considerable levels only in the 17 MH-LH Glut class (≈50% prevalence). This was reflected by high levels of transcription "
                          "in the medial habenula (MH, Figure S7D-E), a structure involved in the response to stress and fear "
                          "{Chou, 2016 #2913;Soria-Gomez, 2015 #2910;Winson, 1978 #2908;Yamaguchi, 2013 #2909}. "
                          "Some transcription was also evident in the posterior part of the brain, "
                          "specifically in the inferior olivary complex (IO), driven by a single supertype, 253 IO Fgl2 Glut, and some structures populated by Sero neurons. \n",

           "Htr2": "The Htr2 family is mainly linked to Gq/11 and causes depolarization by increasing intracellular Ca2+. "
                   "HtraA antagonists, such as atypical antipsychotics (e.g., clozapine and risperidone), are used in treating schizophrenia and other psychiatric disorders. Htr2c"
                   " antagonists are being explored for their potential in treating obesity and metabolic disorders {He, 2022 #2942; Yao, 2021 #2943}. Htr2a, "
                   "instrumental in mediating the effects of psychedelics {Nichols, 2016 #854; De Filippo, 2024 #2904},"
                   " is found across the brain with highest prevalence in telencephalic neighborhoods, Pallium-Glut and Subpallium-GABA (Figure 5A). "
                   "Colocalization was highest with Htr1f and Htr2c (Figure 5B). Highest transcription (≈40%) was found in 01 IT-ET Glut, 07 CTX-MGE GABA "
                   "and 16 HY-MM Glut classes (Figure 5C). Interestingly somatotatin (Sst) neuron belonging to 07 CTX-MGE GABA, "
                   "while exhibiting a relatively low prevalence at the subclass level, contained various clusters "
                   "with >70% prevalence {De Filippo, 2024 #2904}. Htr2a was also prevalent across many other classes across the whole brain. "
                   "01 IT-ET Glut exhibited by far the highest absolute number "
                   "of neurons transcribing Htr2a, specifically in subclasses of L5 and CLA, resembling Htr1f (see interactive visualizer, 'Overview genes by class'). "
                   "CTXsp showed the highest prevalence, reaching >12% (Figure 5D). Isocortex and STR exhibited both ≈5% prevalence. "
                   "At a structure level, two structures belonging to the mammillary complex "
                   "(dorsal premammillary nucleus, PMd and tuberomammillary nucleus, TMd) were in the top ten by prevalence. The mammillary complex "
                   "has been linked to Alzheimer´s disease {Huang, 2023 #2915}, "
                   "and memory {Roy, 2017 #2916}. CLA and the EPd showed the highest absolute prevalence. Interestingly, CLA has been proposed to play an "
                   "important role in mediating the effects of psychedelic compounds {Doss, 2022 #2917}. "
                   "Prevalence in the STR was driven by the small bed nucleus (BA), a structure important for the "
                   "integration of limbic and environmental informations {Lebow, 2016 #2931}. "
                   "Htr2a transcription in CLA and mammillary complex is in line with a previous report in monkey {López-Giménez, 2001 #2930}. "
                   "Prevalence of Htr2a was highest in frontal regions of the brain, "
                   "decaying linearly to virtual absence in the cerebellum (Figure 5E-F). \n"
                   
                   "Htr2b was found only in a minority of neurons and was not included in the MERFISH dataset. No cluster was found to be enriched with Htr2b. "
                   "Interestingly, neurons belonging to "
                   f"the Pineal Glut class showed the highest prevalence at {round(joined.groupby('class')['Htr2b'].apply(percentage_above_threshold)['25 Pineal Glut'], 2)}% "
                   f"(Figure S8C). \n"
                   
                   f"Htr2c was found at highest prevalence in the MB-HB-Glut-Sero-Dopa and Hy-EA-Glut-Gaba neighborhoods (Figure 6A). "
                   f"Apart from Pallium-Glut, its prevalence was always >40%. "
                   f"Colocalization was highest with Htr1f, Htr4 and Htr7 (Figure 6B). "
                   f"Transcription was broadly distributed across many different classes, especially subcortically (Figure 6C). "
                   f"Many classes exhibited a >60% prevalence. "
                   "As usual, cortical excitatory neurons exhibited the highest absolute number "
                   f"of cells transcribing Htr2c. Some subclasses in OLF, amygdala and retrosplenial cortex (RSP) "
                   f"exhibited >80% prevalence (see interactive visualizer, 'Overview genes by class'). "
                   f"The majority of cells transcribing Htr2c RNA belonged to enriched clusters. "
                   f"Highest prevalence was found in STR. Similarly to Htr1b, ACB, CP and OT exhibited the highest prevalence (Figure 6D-E-F). "
                   f"Isocortex prevalence derived from the "
                   f"unique transcription in excitatory neurons of the ventral part of the RSP, "
                   f"curiously the area with lowest transcription of Htr1f, otherwise highly prevalent in all other cortical regions. Htr2a RNA was also minimally expressed in this specific area. "
                   f"High prevalence was observed also in excitatory neurons of "
                   f"the anterior olfactory nucleus (AON), piriform area (PIR and PAA) and amygdala (LA and BLA). "
                   f"Htr2c RNA was found across a variety of structures also in the MB (non in Sero neurons), pons (P), medulla (MY) and cerebellum (CB).  ",

           "Htr4, Htr6 and Htr7": "These receptors are all connected to Gs {McCorvy, 2015 #2889}, leading to increasing cellular levels of cAMP. Htr4 modulation in HPF has been found to "
                                  "bidirectionally influence memory formation in  mice {Teixeira, 2018 #924}. "
                                  "Htr4, similarly to Htr2C, showed highest prevalence (>40%) in the MB-HB-Glut-Sero-Dopa and Hy-EA-Glut-Gaba groups (Figure 7A). "
                                  "It colocalized the most with Htr2c and Htr1f (Figure 7B). "
                                  "Transcription across classes was broadly distributed, with many subcortical classes showing a prevalence >40% (Figure 7C). "
                                  "Highest prevalence was found in "
                                  "the 17 MH-LH Glut class, specifically in the Chol releasing neurons belonging to this class located in TH. "
                                  "In absolute numbers, transcription in excitatory cortical neurons was comparable to other classes but still the highest, driven specifically "
                                  "by subclasses of CA1, CA2, CA3 and subiculum (see interactive visualizer, 'Overview genes by class').  "
                                  "Spatial distribution exhibited a peculiar pattern "
                                  "with high prevalence in one specific structure of the STR: OT (Figure 7D-E-F). A subclass of interneurons present in OT (060 OT D3 Folh1 Gaba) "
                                  "showed a >98% prevalence. "
                                  "PAL and HPF also exhibited relatively high prevalence (≈10%). "
                                  "Dentate gyrus (DG) granule cells (037 DG Glut) were one of the reasons of the high prevalence in HPF. "
                                  "Excitatory cells of CA2, CA3 and indusium griseum (IG) also transcribed Htr4 RNA (see interactive visualizer, 'Overview genes by brain structure'). \n"
                                  
                                  "We do not have MERFISH information about the rarely transcribed Htr6 and no enriched cluster was "
                                  f"present in the scRNA-seq dataset. The 09 NU-LGE GABA class exhibited the highest prevalence with "
                                  f"{round(joined.groupby('class')['Htr6'].apply(percentage_above_threshold)['09 CNU-LGE GABA'], 2)}, still, the absolute majority of neurons "
                                  f"transcribing the RNA of this gene were excitatory cortical neurons (Figure S9C). \n"
                                  
                                  "Conversely, Htr7 was transcribed in >10% of the totality of cells. "
                                  "It reached ≈60% in the TH-EPI Glut group, and considerable amounts (≈40%) in MB, HB and HY groups (Figure 8A). "
                                  "Colocalization was the highest with Htr2c and Htr1f (Figure 8B). "
                                  "Transcription was broadly distributed across classes present in HY, MB and TH (Figure 8C). "
                                  "It colocalized the most with Htr2c, "
                                  "Htr1f and Htr4. Htr7 was broadly transcribed across classes, especially in subcortical structures. "
                                  "Peak prevalence was found in 10 LSX GABA, 16 MY MM Glut and 18 TH Glut with >60% "
                                  "(Figure 8C). Cortical transcription in excitatory neurons is driven primarily by subclasses in CA2 and L2 ENT (see interactive visualizer, 'Overview genes by class'). "
                                  "Htr7 enriched clusters were located mainly in HY and TH (Figure 8D). "
                                  "At a structure level, the parafascicular (PF)  and paraventricular nucleus (PVT) of TH showed the highest prevalence (>30%). ",

           "Htr3": "The Htr3 family is the only ionotropic Htr and it causes direct excitation by allowing the influx of cations. "
                   "Htr3 antagonists, such as ondansetron, are effective antiemetics used to prevent nausea and vomiting."
                   "The Htr3a subunit is required for the formation of a functional channel {Maricq, 1991 #2918} and can form functional homopentameric receptors "
                   "{Walstab, 2010 #2919}. "
                   "Heteromeric receptors containing Htr3b have an increased channel conductance and different selectivity {Davies, 1999 #2920}. "
                   "Htr3a is transcribed almost uniquely in the "
                   "Subpallium-Gaba neighborhood, with a prevalence of ≈8% (Figure 9A), specifically in the 06 CTX-CGE GABA class (Figure 9C). "
                   "It is one of the few Htr, together with Htr3b and Htr1d, "
                   "that is not transcribed the most in absolute numbers in 01 IT-ET glut. It colocalizes mainly "
                   "with Htr2c and Htr7 (Figure 9B). This Htr was mainly transcribed in OLF, CTXsp, HPF and Isocortex (Figure 9D) and "
                   "is most prevalent in the anterior part of the brain, although, puzzlingly, with slightly lower amount of RNA per cell (Figure 9E-F). "
                   "Htr3b was not included in the MERFISH dataset, and no cluster was found to be enriched with this receptor. Htr3b was the least "
                   "transcribed Htr gene in the entire RNAseq dataset. "
                   "Similarly to Htr3a, its transcription was delimited to the 06 CTX-CGE GABA class (Figure S10C)."


           }

if Adapt_for_Nature_style is True:
    results = Naturize_text(results)
#%%
