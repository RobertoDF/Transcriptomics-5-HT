
from Utils.Results_variables import *

abstract = ("Serotonin (5-HT) plays a pivotal role in regulating a wide range of brain functions, including mood, sleep, and cognition. "
            "This study presents a comprehensive transcriptomic analysis of 5-HT "
            "receptors (Htrs) covering ~4 million cells across the whole adult mouse brain, utilizing single-cell "
            "RNA sequencing (scRNA-seq) data from the Allen Institute. We report on the differential expression"
            " patterns of all 14 known Htr subtypes, revealing a wide diversity in their prevalence and "
            f"distribution throughout the brain. Notably, we found that {round(((exp>threshold_expression).sum(axis=1).astype(bool).sum()/exp.shape[0])*100,2)}%"
            f" of cells transcribe RNA of at least one Htr and often "
            "Htrs were found to colocalize . "
            "The transcription patterns of Htrs can accurately inform a random forest classifier to identify specific classes and types "
            "of neurotransmitter-releasing cells with surprising success. "
            "Leveraging a a multiplexed error-robust fluorescence in situ hybridization "
            "(MERFISH) dataset provided by Harvard University of ~7% of the totality of cells found in a mouse brain, "
            "we analysed "
            "the spatial distribution of each Htr confirming previous findings and uncovering novel patterns of possible expression at an unprecedented level of detail. "
            "Our findings underscore the complexity of the 5-HT system "
            "even at the single-cell dimension and"
            " provide new insights into the receptor-mediated mechanisms that underpin diverse neural functions and behaviors. To aid the exploration of Htrs "
            "transcription in the datasets "
            "we provide a custom online visualizer. This tool enables in-depth analysis at various levels of granularity. ")
