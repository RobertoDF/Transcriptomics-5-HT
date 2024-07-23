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

doc = DocxTemplate(f"{manuscript_folder}/Manuscript_supplementary_template.docx")

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

correspondence_to = "roberto.de-filippo@bccn-berlin.de"


context = {'title': title, "authors": authors, "affiliations": affiliations, "correspondence_to": correspondence_to,
           "supplementary_legends": legends_supplementary,
}

doc.render(context, autoescape=True)

doc.save(f"{manuscript_folder}/Transcriptomics_5-HT.docx")


#%%
