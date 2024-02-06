# Transcriptomic Mapping of the 5-HT Receptor Landscape

[Online visualizer](https://huggingface.co/spaces/RDeF654875678597657/5-HT-Transcriptomics)

## Description

Our entire analysis is contained in 2 Jupyter notebooks located in the folder 'Figures': 'Figure_1.ipynb' and 'Figure_2.ipynb'. To adapt the code for the visualization of different genes, it is sufficient to change the'family_name' and 'genes_families' variables in the Utils.Settings.py file. Data is downloaded following the [instructions](https://alleninstitute.github.io/abc_atlas_access/intro.html) provided by the Allen Institute. Notebooks to download the RNA-seq and MERFISH datasets are contained in the 'Load_Data' folder. To explore the expression of different genes, it is necessary to download the associated expression matrices by changing the selected genes in the 'Download_RNAseq_data.ipynb' notebook.This can be achieved by modifying the cells underneath the headings 'Select genes RNA-seq' and 'Select genes MERFISH'.

## Installation
```
conda create --name Transciptomics-5-HT --file requirements.txt
```
## Contributing

We welcome contributions from the community! If you'd like to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes.
4. Push your branch and submit a pull request.

## License

cc-by-nc-nd-4.0
