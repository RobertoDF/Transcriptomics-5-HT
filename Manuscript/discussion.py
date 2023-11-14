from Utils.Results_variables import r_common, r_strong

discussion = "Our results show for the first time that strong ripples propagate differentially along the hippocampal " \
             "longitudinal axis. This propagation idiosyncrasy can be explained by a specific ability of the hippocampal septal pole " \
             "(medial section in our analysis) to produce longer " \
             "ripples that better entrain the hippocampal network and spread across the longitudinal axis. " \
             "It was previously observed that ripples located at the septal and temporal pole are generated independently from each other, " \
             "in addition, despite the presence of connections within the hippocampal longitudinal axis {Witter, 2007 #1155; van Strien, 2009 #1156}, " \
             "in the vast majority of cases ripples do not propagate to the opposite pole {Sosa, 2020 #1154}. "\
             "In accordance with these results, " \
             "we observed a strong effect of spatial distance on ripple strength correlation confirming a previous study {Nitzan, 2022 #1157}: " \
             "the strength correlation, predictably, was higher in CA1 pairs " \
             "closer to each other. The effect of distance was also apparent on the ripple chance of propagation, only half of the ripples generated in " \
             "the septal pole were detected additionally in the intermediate hippocampus (lateral section in our analysis). " \
             "This chance is much higher compared to the ~3.7% reported regarding propagation between " \
             "opposite poles {Sosa, 2020 #1154}, it would be interesting to understand whether the temporal pole is also able to entrain the intermediate hippocampus " \
             "in similar fashion or it is a peculiarity of the septal pole. " \
             "A limitation of our work derives from the dataset being limited to the septal and intermediate hippocampus. \n" \
             "Ripples can arise at any location along the hippocampal longitudinal axis {Patel, 2013 #1133}. " \
             "Our analysis shows that ripples are, however, " \
             "not homogeneously generated across space. We observed important differences between strong ripples and common ripples generation. " \
             "Common ripples followed a gradient with higher generation probability in the intermediate section and lowest in the septal pole. " \
             "Strong ripples, on the other hand, were mostly generated locally (i.e. a strong ripple detected in the medial section is most likely " \
             "generated in the medial section itself). Furthermore, only rarely a strong ripple generated in the intermediate hippocampus is able to propagate " \
             "towards the septal pole retaining its strong status (top 10%). Conversely strong ripples generated in the septal pole have" \
             " a significantly higher chance of " \
             "propagate longitudinally and still be in the top 10% in terms of ripple strength. Notably, this is not consequence of a simple " \
             "longitudinal gradient in ripple strength, indeed, we did not observe any" \
             " difference in ripple strength along the longitudinal axis. " \
             "Additionally, we show that ripples generated in the septal pole and in the " \
             "intermediate hippocampus have a significantly different ability to engage hippocampal networks in the 50-120 ms window " \
             "post ripple start. Ripples generated in the septal pole activate more neurons, both excitatory and inhibitory, and, moreover, " \
             "elicit an higher spiking rate per neuron. This prolonged network activation is reflected by the fact that the position on the longitudinal axis explains " \
             f"{round(r_common*100, 2)}% and {round(r_strong*100, 2)}% of the variability in ripple duration in common and strong ripples respectively. " \
             "Consistent with a duration gradient along the longitudinal axis, the temporal hippocampus has been shown to " \
             "produce shorter ripples both in awake and sleep conditions {Sosa, 2020 #1154}. \n"\
            "What is the reason that enables the septal pole to generate longer ripples? There might be for example underlying electrophysiological" \
             " differences between the septal and intermediate hippocampus. Upon closer examination of the electrophysiological features of the neurons, " \
             "we were able to discern significant differences in the shape and duration of their waveform. We can hypothesize that slower action potentials and, consequentially, " \
             "longer refractory periods hinder the ability" \
             " to sustain protracted high frequency spiking. Accordingly, we found an increased firing rate and a smaller waveform duration in putative excitatory " \
             "neurons of the septal pole. Moreover, putative inhibitory neurons in the septal pole showed reduced firing. " \
             "These differences might contribute to explain the prolonged ripples observed in the septal pole. " \
             "We can also speculate that the neuromodulatory inputs gradient, monoamine fibers have been shown to be stronger in the ventral part {Strange, 2014 #1158}, " \
             "might influence neurons responses. Serotonin {ul Haq, 2016 #1165; Wang, 2015 #1170}, noradrenaline {Ul Haq, 2012 #1166; Novitskaya, 2016 #1171} and acetylcholine {Zhang, 2021 #1167} " \
             "have all been shown to suppress ripples. In accordance with this, some ripples are coupled with a reduced activation of the locus coeruleus" \
             " and the dorsal raphe nucleus in vivo {Ramirez-Villegas, 2015 #1168}. \n" \
             "Ripples can be subdivided in different types according to the relationship between the hippocampal LFP and the ripple itself " \
             "{Ramirez-Villegas, 2015 #1168}. " \
             "Intriguingly these subtypes are associated with two different brain-wide networks, the first communicating preferentially with " \
             "the associative neocortex and a second one biased towards subcortical structures. Moreover, these different types of ripples have been " \
             "proposed to possibly fulfill different functional roles. Given the different input/output connectivity between " \
             "septal, intermediate and temporal hippocampus {Fanselow, 2010 #1150} we hypothesize that ripple generated at different points of " \
             "the hippocampal longitudinal axis might as well have functional differences, with the longer ripples generated septally possibly able to " \
             "combine the different kind of informations processed in the distinct hippocampal sections and additionally relaying the " \
             "integrated information back to the neocortex in accordance with the two-stage memory hypothesis {Diekelmann, 2010 #1172; Marr, 1971 #1173;" \
             "Buzsáki, 1989 #1174; Rasch, 2007 #1175; McClelland, 1995 #1176}. \n" \
             "Long duration ripples have been shown to be of particular importance in situations of " \
             "high-memory demand {Fernández-Ruiz, 2019 #1121}, at the same time, previous studies highlighted " \
             "the role of septal hippocampus in memory tasks and information processing {Hock, 1998 #1122;Moser, 1993 #1123;Moser, 1995 #1124;Steffenach, 2005 #1128;" \
             "Kheirbek, 2013 #1126;McGlinchey, 2018 #1125;Fanselow, 2010 #1150;Maras, 2014 #1130;Bradfield, 2020 #1131;Qin, 2020 #1132}. " \
             "Our results can contribute to explain the specific role of septal hippocampus in memory-demanding tasks with" \
             " its ability of generating particularly long " \
             "ripples that are able to strongly engage networks in the entire top half of the hippocampal formation for an extended time. \n" \
             "Functional differences between ripples are supported by our finding that the ripple origination point has an influence on the " \
             "engagement of various hippocampal subfields. Interestingly, these functional differences can be observed even before the onset of ripples. " \
             "Our analysis revealed that DG, CA1, and CA3 subfields of the hippocampus were more strongly engaged during the early phase " \
             "in response to medial ripples, the opposite was true for SUB. In the late phase, all subfields, including ProS, " \
             "showed stronger modulation in response to medial ripples. Moreover, the analysis of pre-ripple activity showed that ProS and " \
             "SUB exhibited stronger modulation in response to medial ripples, while DG and, especially, CA3 displayed weaker modulation compared to lateral ripples.  " \
             "These results indicate that CA3 may play a more critical role in initiating lateral ripples, in accordance  with the traditional view " \
             "that CA3 is the primary " \
             "generator of ripples {Buzsáki, 1986 #1196; Buzsáki, 1989 #1197; Csicsvari, 2000 #1195}. On the other hand, " \
             "medial ripples demonstrate stronger" \
             " engagement of ProS and " \
             "SUB subfields before ripple start, supporting the more recent " \
             "hypothesis that output structures in the hippocampus are also capable of generating ripples {Imbrosci, 2021 #1194}. " \
             "The subicular complex has been suggested to play a role in the transfer of information from the hippocampus to other " \
             "brain regions, such as the cortex {Aggleton, 2015 #1199; Naber, 1998 #1200}. It is possible that the weaker engagement of " \
             "the subiculum during lateral ripples " \
             "may reflect a more localized processing of information within the hippocampus, while the stronger engagement of the " \
             "subiculum during medial ripples may reflect a more global transfer of information to other brain regions. "\
             "Differences in ripple initiation and engagement patterns can provide valuable insights into the mechanisms that" \
             " underlie the dynamics of the hippocampal network.  " \
             "Gaining a better understanding of these differences may shed light on the functional significance " \
             "of ripples in the brain and their role in memory consolidation and retrieval. \n" \


