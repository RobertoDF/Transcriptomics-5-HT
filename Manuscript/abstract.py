
from Utils.Results_variables import *

abstract = ("Serotonin (5-HT) is crucial for regulating brain functions such as mood, sleep, and cognition. "
            "This study presents a comprehensive transcriptomic analysis of 5-HT receptors (Htrs) across ≈4 million cells in the adult mouse brain "
            "using single-cell RNA sequencing (scRNA-seq) data from the Allen Institute. "
            "We observed differential transcription patterns of all 14 Htr subtypes, revealing diverse prevalence and distribution across cell classes. "
            "Remarkably, we found that "
            f"{round(((exp>threshold_expression).sum(axis=1).astype(bool).sum()/exp.shape[0])*100,2)}% "
            f"of cells transcribe RNA of at least one Htr, with frequent co-transcription of multiple Htrs, "
            "underscoring the complexity of the 5-HT system "
            "even at the single-cell dimension. "
            "Leveraging a multiplexed error-robust fluorescence in situ hybridization "
            "(MERFISH) dataset provided by Harvard University of ≈10 million cells, "
            "we analyzed the spatial distribution of each Htr, confirming previous findings and uncovering novel transcription patterns. "
            "To aid in exploring Htr transcription, we provide an online interactive visualizer (https://rdef654875678597657-5-ht-transcriptomics.hf.space).")
