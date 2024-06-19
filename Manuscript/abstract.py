
from Utils.Results_variables import *

abstract = ("Serotonin (5-HT) plays a pivotal role in regulating a wide range of brain functions, including mood, sleep, and cognition. "
            "This study presents a comprehensive transcriptomic analysis of 5-HT "
            "receptors (Htrs) covering ≈4 million cells across the whole adult mouse brain, utilizing single-cell "
            "RNA sequencing (scRNA-seq) data from the Allen Institute. We report on the differential transcription"
            " patterns of all 14 known Htr subtypes, revealing a wide diversity in their prevalence and "
            f"distribution across cell classes.  Remarkably, we found that "
            f"{round(((exp>threshold_expression).sum(axis=1).astype(bool).sum()/exp.shape[0])*100,2)}% "
            f"of cells transcribe RNA of at least one Htr, and co-transcription of multiple Htrs within single cells was frequently observed. "
            "The transcription patterns of Htrs can accurately inform a random forest classifier to identify specific classes and types "
            "of neurotransmitter-releasing cells with surprising success. "
            "Leveraging a multiplexed error-robust fluorescence in situ hybridization "
            "(MERFISH) dataset provided by Harvard University of ≈10 million cells found in a mouse brain, "
            "we analyzed "
            "the spatial distribution of each Htr confirming previous findings and uncovering novel patterns of transcription "
            f"at an unprecedented level of detail.  We show that the majority of Htr-transcribing cells "
            f'{round(((joined_boolean[joined_boolean["Number of Htrs"]>0]["Number of Htrs"].value_counts()/joined_boolean[joined_boolean["Number of Htrs"]>0].shape[0])*100).loc[5:].sum(), 2)}% '
            "contain RNA of at least one other Htr, "
            "underscoring the complexity of the 5-HT system "
            "even at the single-cell dimension. "
            "To aid the exploration of Htrs "
            "transcription in the datasets "
            "we provide an interactive visualizer available online (https://rdef654875678597657-5-ht-transcriptomics.hf.space). This tool enables in-depth analysis at various levels of granularity. ")
