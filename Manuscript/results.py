from Utils.Results_variables import *
from Utils.Settings import Adapt_for_Nature_style
from Utils.Utils import Naturize_text

results = {"Transcriptomic overview of 5-HT receptors landscape": ""



           }

if Adapt_for_Nature_style is True:
    results = Naturize_text(results)