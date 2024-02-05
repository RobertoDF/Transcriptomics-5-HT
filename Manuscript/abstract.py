
from Utils.Results_variables import *

abstract = ("Serotonin (5-HT) plays a pivotal role in regulating a wide range of brain functions, including mood, anxiety, sleep, and cognition. "
            "This study presents a comprehensive transcriptomic analysis of the 5-HT "
            "receptors (Htrs) transcription across approximately 4 million brain cells, utilizing single-cell "
            "RNA sequencing (scRNA-seq) data from the Allen Institute. We report on the differential expression"
            " patterns of all 14 known Htr subtypes, revealing a wide diversity in their prevalence and "
            f"distribution throughout the brain. Notably, we found that {round(((exp>threshold_expression).sum(axis=1).astype(bool).sum()/exp.shape[0])*100,2)}%"
            f" of cells transcribe RNA of at least one Htr subtype, "
            "with significant variability in prevalence across different depth of classification levels and distribution across brain regions. "
            "Htrs were found to colocalize in the vast majority of cases. "
            "The transcription patterns of Htrs can accurately inform a random forest classifier to identify specific classes and types "
            "of neurotransmitter-releasing cells with surprising success. Leveraging a a multiplexed error-robust fluorescence in situ hybridization (MERFISH) dataset "
            "provided by Harvard University we analysed "
            "the spatial distribution of each Htr confirming previous findings and uncovering novel patterns of possible expression at an unprecedented detailed level. "
            "Our findings underscore the complexity of the 5-HT system "
            "even at the single-cell dimension and"
            " provide new insights into the receptor-mediated mechanisms that underpin diverse neural functions and behaviors. To aid the exploration of the dataset "
            "at different level of granularity "
            "we provide a custom online visualizer.")
