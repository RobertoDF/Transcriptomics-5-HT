from Utils.Results_variables import *
from Utils.Settings import Adapt_for_Nature_style
from Utils.Utils import Naturize_text

results = {"Transcriptomic overview of 5-HT receptors landscape":
               "We analysed the single-cell RNA-sequencing (scRNA-seq) dataset provided by Allen Institute {Yao, 2023 #2828} "
                "focusing on the expression of Htrs across 4 million cells. The scRNA-seq dataset contained informations about all the known "
                f"14 Htr. Prevalence of Htrs was considerably different ranging from {expression_total.min()} % of {expression_total.idxmin()} "
               f"to {expression_total.max()} % of {expression_total.idxmax()} (Figure 1A). RNA of 6 Htr was found in less than 2.5% of the cells "
               f"({', '.join(expression_total[expression_total < 2.5].index.values)}). On the other hand, RNA of Htr1f, Htr2a and Htr2c was present in"
               f"at least 1 every 5 cells. Beside the amount of expression, also the distribution among classes was considerably different. This is exemplified"
               f"by looking at the distribution of the Htr1 and Htr2 families across different cluster groups (Figure 1B). We could see clear areas of overlap and "
               f"separation in cortical neurons, constituted by the Pallium-Glut and Subpallium-GABA groups. Distribution within family also showed considerable differences (Supplementary Figure 1)."
               f"Htrs distribution was also markedly different across neurons "
               f"releasing different neurotransmitters (Figure 1C). Cells not found to express any transmitter made out {total_cells_by_neurotransmitter[""]} % of the total (Supplementary Figure 2). "
               f"This is the only group that did not express significant (prevalence<20%) amounts of any Htrs. All other groups expressed significant amounts of at least 2 different Htrs. "
               f"Expectedly the vast majority of cells is classified as excitatory ({total_cells_by_neurotransmitter["Glut"]} %). "
               f"Around 1 every 5 cells was found to release GABA ({total_cells_by_neurotransmitter["GABA"]} %). All the other neurotransmitter are found in less than 1% of the cells, "
               f"in particular, 5-HT releasing neurons  were found in  {total_cells_by_neurotransmitter["Sero"]} % of the cells. 5-HT neurons expressed a highest variety of Htrs. They show "
               f"the highest prevalence for all the receptors belonging to the Ht1 family (Htr1a: {expression_by_neurotransmitter["Htr1a"]["Sero"]} %, "
               f"Htr1b: {expression_by_neurotransmitter["Htr1b"]["Sero"]} %, Htr1d: {expression_by_neurotransmitter["Htr1d"]["Sero"]} % and Htr1f: {expression_by_neurotransmitter["Htr1f"]["Sero"]} %). "
               f"They also show significant amounts of Htr2a ({expression_by_neurotransmitter["Htr2a"]["Sero"]} %), Htr2c ({expression_by_neurotransmitter["Htr2c"]["Sero"]} %), "
               f"Htr4 ({expression_by_neurotransmitter["Htr4"]["Sero"]} %), Htr5a ({expression_by_neurotransmitter["Htr5a"]["Sero"]} %), Htr5b ({expression_by_neurotransmitter["Htr5b"]["Sero"]} %) "
               f"and Htr7 ({expression_by_neurotransmitter["Htr7"]["Sero"]} %). In total, 5-HT neurons showed significant expression of 8 different Htrs. GABA-Glyc neurons, constituting only "
               f"{total_cells_by_neurotransmitter["GABA-Glyc"]} % of cells, showed significant amounts of 7 different Htrs, with particualrly high prevalence of Hr1d "
               f"({expression_by_neurotransmitter["Htr1d"]["GABA-Glyc"]} %), Htr7 ({expression_by_neurotransmitter["Htr7"]["GABA-Glyc"]} %), Htr1f ({expression_by_neurotransmitter["Htr1f"]["GABA-Glyc"]} %) and "
               f"Htr4 ({expression_by_neurotransmitter["Htr4"]["GABA-Glyc"]} %). Cells expressing GABA ({total_cells_by_neurotransmitter["GABA"]} % of cells) show significant expression of "
               f"Htr2c ({expression_by_neurotransmitter["Htr2c"]["GABA"]} %), Htr1f ({expression_by_neurotransmitter["Htr1f"]["GABA"]} %), Htr7 ({expression_by_neurotransmitter["Htr7"]["GABA"]} %) and "
               f"Htr2c ({expression_by_neurotransmitter["Htr2a"]["GABA"]} %). GABA-Glyc neurons showed a similar pattern with, notably, the higher prevalence of  "
               f"Htr2c ({expression_by_neurotransmitter["Htr2a"]["GABA-Glyc"]} %) and Htr7 ({expression_by_neurotransmitter["Htr7"]["GABA-Glyc"]} %). Cholinergic neurons distinguish themself "
               f"by exhibiting the highest prevalence of Htr4 ({expression_by_neurotransmitter["Htr4"]["Chol"]} %) and Htr5b ({expression_by_neurotransmitter["Htr5b"]["Chol"]} %). "
               f"GLutamatergic neurons show significant expression of Htr1f ({expression_by_neurotransmitter["Htr1f"]["Glut"]} %), Htr2c ({expression_by_neurotransmitter["Htr2c"]["Glut"]} %) and "
               f"Htr2a ({expression_by_neurotransmitter["Htr2a"]["Glut"]} %). Dopaminergic neurons show a similar pattern with lower Htr2a and higher Htr7 ({expression_by_neurotransmitter["Htr7"]["Dopa"]} %)). "
               f"At last, Histamine neurons express significant amounts of Htr2c ({expression_by_neurotransmitter["Htr2c"]["Hist"]} %) and Htr4 ({expression_by_neurotransmitter["Htr4"]["Hist"]} %), "
               f"Noradrenergic neurons instead show high prevalence of Htr1f ({expression_by_neurotransmitter["Htr1f"]["Nora"]} %) and Htr5a ({expression_by_neurotransmitter["Htr5a"]["Hist"]} %). "
               f"Looking at expression across groups described in Figure 1B, we notice that non-neuronal cells (NN-IMN-GC) show the lowest expression, mirroring the data regarding cells without "
               f"any neurotransmitter. Interestingly the patterns of expression were less differentiated across groups "
               f"(Pearson coefficient={round(lower_triangle_groups.stack().mean(),2)}±{round(lower_triangle_groups.stack().sem(),2)}) compared to neurotransmitters "
               f"(Pearson coefficient={round(lower_triangle_neurotransmitter.stack().mean(),2)}±{round(lower_triangle_neurotransmitter.stack().sem(),2)}, Supplementary Figure 3). "
               f"The totality of cells analyzed were divided in 34 classes in the original study. We analyzed expression across this pre-identified classes (Figure 1E). Average correlation between "
               f"patterns of expression was {round(lower_triangle_class.stack().mean(), 2)}±{round(lower_triangle_class.stack().sem(), 2)}. Across classes Htr2c is the one with the highest average prevalence "
               f"({mean_expression_by_class["Htr2c"]}), followed by Htr1f ({mean_expression_by_class["Htr1f"]}), Htr7 ({mean_expression_by_class["Htr7"]}) and Htr4 ({mean_expression_by_class["Htr4"]}). "
               f"Correlation between Htrs expression across the totality of cells ranged from {lower_triangle_corr.min()} ({'-'.join(lower_triangle_corr.stack().idxmin())}) to {lower_triangle_corr.miax()} ({'-'.join(lower_triangle_corr.stack().idxmax())}). "
               f"Considerable correlation was found also for the {'-'.join(lower_triangle_corr.index[1])} (Pearson coefficient={lower_triangle_corr.sort_values(ascending=False).iloc[1]}) and "
               f"{'-'.join(lower_triangle_corr.index[2])} (Pearson coefficient={lower_triangle_corr.sort_values(ascending=False).iloc[2]})) pairs (Figure 1H). Effect of this correlation was also visible "
               f"when looking at co-expression (Figure 1H). Expectedly, Htr1f and Htr2c, the most prevalent Htrs, were found to co-localize with other receptors respectevely {mean_coloc["Htr1f"]} % "
               f"and  {mean_coloc["Htr2c"]} % of the times. Only rarely a cell was found to express only one Htr, {round(at_least_2_receptors.mean(), 2)}±{round(at_least_2_receptors.sem(), 2)} % of cells indeed expressed at least "
               f"2 Htrs (Figure 1G).Surprisongly, {round(at_least_5_receptors.mean(), 2)}±{round(at_least_5_receptors.sem(), 2)} % of cells expressed at least 5 Htrs. This is indicative of the complexity "
               f"of the 5-HT system even at a single cell level. The highest amount of co-localization (at least 2 Htrs) was present in the GABAergic neurons of midbrain, hindbrain, and cerebellum (MB-HB-Glut-Sero-Dopa, "
               f"{at_least_2_receptors_per_group["MB-HB-Glut-Sero-Dopa"]} %). The average excluding non neuronal cells was {round(at_least_2_receptors_per_group.drop("NN-IMN-GC").mean(), 2)}±{round(at_least_2_receptors_per_group.drop("NN-IMN-GC").sem(), 2)} %. "
               f"The extensive expression across different classes and the considerable coexpression within cells point at the complexity of the 5-HT sistem. In the following sections we will take a deeper look "
               f"on each receptor with a particular focus on the one more prevalent",
           "Htr1a":""




           }

if Adapt_for_Nature_style is True:
    results = Naturize_text(results)