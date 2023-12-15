from docxtpl import DocxTemplate
from Manuscript.abstract import abstract
from Manuscript.results import results
from Manuscript.discussion import discussion
from Manuscript.introduction import introduction
from Manuscript.methods import methods
from Figures.legends import legends
from Figures.legends_supplementary import legends_supplementary
from Utils.Settings import manuscript_folder
# to work with citations use {Abi-Saab, 1999 #888}. A Endnote travelling library is provided in the manuscript folder.

doc = DocxTemplate(f"{manuscript_folder}/Manuscript_template.docx")

title = "Transcriptomic Mapping of the 5-HT Receptor Landscape"
authors = "Roberto De Filippo¹ and Dietmar Schmitz¹²³⁴⁵"
affiliations =  "¹ Charité Universitätsmedizin Berlin, corporate member of Freie Universität Berlin, Humboldt-Universität" \
                " zu Berlin, and Berlin Institute of Health; Neuroscience Research Center, 10117 Berlin, Germany. \n" \
                "² German Center for Neurodegenerative Diseases (DZNE) Berlin, 10117 Berlin, Germany. \n"\
                "³ Charité-Universitätsmedizin Berlin, corporate member of Freie Universität Berlin, Humboldt-Universität Berlin, and "\
                "Berlin Institute of Health, Einstein Center for Neuroscience, 10117 Berlin, Germany. \n"\
                "⁴ Charité-Universitätsmedizin Berlin, corporate member of Freie Universität Berlin, Humboldt-Universität Berlin, " \
                "and Berlin Institute of Health, NeuroCure Cluster of Excellence, 10117 Berlin, Germany. \n"\
                "⁵ Humboldt-Universität zu Berlin, Bernstein Center for Computational Neuroscience, Philippstr. 13, 10115 Berlin, Germany. "\


correspondence_to = "roberto.de-filippo@charite.de"
keywords = "5-HT receptors, serotonin"

acknowledgements = "This study was supported by the German Research Foundation Deutsche Forschungsgemeinschaft (DFG), " \
                   "project 184695641 - SFB 958, project 327654276 - SFB 1315, Germany's Excellence Strategy - Exc-2049-390688087 and by the " \
                   "European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (Grant agreement No. 810580). " \
                "We thank J.T. Tukker, N. Maier for feedback on an early version of the manuscript and the members of the Schmitz lab for " \
                   "scientific discussion. We thank Willy Schiegel and Tiziano Zito for technical help with cluster computing. We thank Federico Claudi " \
                   "for support with brainrender. " \
                "The authors declare that they have no competing interests. "

"""contributions = "Conceptualization, data curation, formal analysis, investigation, visualization: RDF. Writing - original draft: RDF. " \
                "Writing - review & editing: RDF, DS. " \
                "Funding acquisition: DS."
"""

data_availability = "All the code used to process the dataset is available at https://github.com/RobertoDF/De-Filippo-et-al-2022, pre-computed data structures "\
                    "can be downloaded at 10.6084/m9.figshare.20209913. "\
                    "All figures and text can be reproduced using code present in this repository, each number present in the text is directly "\
                    "linked to a python data structure. The original dataset is provided by the Allen Institute and available at "\
                    "https://allensdk.readthedocs.io/en/latest/visual_coding_neuropixels.html."

context = {'title': title, "authors": authors, "affiliations": affiliations, "correspondence_to": correspondence_to,
           "keywords": keywords, "abstract": abstract, "introduction": introduction,
           "discussion": discussion, "methods": methods, "results": results, "legends": legends, "supplementary_legends": legends_supplementary,
           "acknowledgements": acknowledgements,
          "data_availability": data_availability}

doc.render(context, autoescape=True)

doc.save(f"{manuscript_folder}/De Filippo et al., 2022.docx")

