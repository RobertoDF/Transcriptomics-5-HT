from Utils.Results_variables import *

discussion = ("In this work we described the main transcriptional features of each Htr across the entire brain, "
              "leveraging two datasets provided by {Yao, 2023 #2828} and {Zhang, 2023 #2887}. "
              f"We found that Htrs RNA is transcribed in around 2 every 3 cells and 6 Htrs were transcribed in >10% of cells, with Htr1f reaching a peak of {expression_total.max()}%. "
              f"Htrs patterns of transcription can be used to decode the identity of cells grouped by neurotransmitter, neighborhoods and classes at an above chanche level. "
              f"Surprisingly, it was common to detect multiple Htrs within a single cells. This points at the great complexity of the 5-HT "
              f"system even at the a unicellular level."
              "We can recapitulate our results regarding each Htr by summarizing the defining feature of each receptor: "
              "Htr1a is expressed in an important fraction of Sero neurons of the raphe and some HPF excitatory neurons; "
              "Htr1b is expressed in many inhibitory striatal neurons and Sero neurons; "
              "Htr1d, similarly to Htr1b, is expressed in the striatum, although at much lower levels; "
              "Htr1f is widely expressed in telencephalic structures,especially the Isocortex, with a peak in frontal olfactory structures; "
              "Htr2a is prevalent in glutamatergic cells of the cortical subplate (CLA and EPd) and the mammillary bodies (TMd, PMd), and hippocampal interneurons; "
              "Htr2b is rarely transcribed and is present in some neurons of the pineal gland; "
              "Htr2c is broadly transcribed, especially in the STR, excitatory neurons of the amygdala (LA, BLA and BMA) and RSPv, OLF neurons and structures in MB, P, MY and CB; "
              "Htr3a and Htr3b are uniquely observed in cortical gabaergic neurons of the 06 CTX-CGE GABA class; "
              "Htr4 is transcribed at high levels in the OT, excitatory cells of the hippocampus proper and DG, and Chol neurons of the TH (17 MH-LH Glut); "
              "Htr5a is transcribed at low levels with only one enriched cluster in the MB; "
              "Htr5b is also transcribed only in few cells, specifically in Chol neurons of the TH; "
              "Htr6 does not feature any enriched cluster, some cells in CA3 transcribed this Htr; "
              "Htr7 is widely transcribed in subcortical structures, especially in some TH nuclei (PF, PVT, IAD and PT), the mammillary complex (MM and PMd), "
              "the lateral septal nucleus (LSv) and the fasciola cinerea of the HPF. "
              "Our analysis is in no way exhaustive and it is limited in scope by the costraints of a traditional scientific article. To bypass this limit and, at the same time, "
              "provide the ability to explore the 5-HT transcription landscape at different depths, we provide a custom online visualizer. The visualizer enbles "
              "the exploration of: Htrs transcription in the MERFISH dataset; the prevalence of each Htr across neighborhoods, class, subclass, supertype and clusters; "
              "an overview of Htrs prevalence across classes and subclasses; "
              "and an overview of Htrs prevalence across all brain divisions and structures optionally filtered by neurotransmitter release. Our entire "
              "analysis pipeline can be easily modified to enable the exploration of different families of genes. "
              "Instructions are available in 'Jupyter notebooks structures' in the methods section. "
              "One constraint of our study is the indirect characterization of Htrs through the detection of RNA molecules, rather than direct assessment of their presence. "
              "However, this potential limitation is mitigated by the fact that mRNA levels are frequently a reliable indicator of receptor expression {Vilaró, 2020 #2939}. Conversely, "
              "while mapping receptors directly allows for precise localization, it fails to differentiate between pre- and postsynaptic expression, "
              "an important aspect of understanding receptor function and distribution. "
              "This lack of specificity becomes particularly problematic, for example, in the context of Sero neurons, which have "
              "extensive projections throughout the brain and exhibit diverse autoreceptors. "
              "This complexity is underscored both in our findings and in previous research, highlighting the intricate regulatory "
              "mechanisms of serotonin neurotransmission {Hjorth, 1991 #2932; Haj-Dahmane, 1991 #2924}. "
              "Our exploration of the Htrs landscape represents a substantial advancement, contributing to our understanding "
              "of the 5-HT system's role in brain function and behavior. ")


