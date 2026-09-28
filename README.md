This repository contains the code and information to reproduce the analysis presented in the manuscript "Genomically discrete and ecologically cohesive subspecies as biologically meaningful units in the global ocean microbiome" by Paoli & Priest et al.

To facilitate usability and reproducibility, we have provided detailed descriptions of the subspecies inference framework and downstream analysis on dedicated wiki pages, which accompany the scripts provided on this repository and the datasets provided on the [Zenodo repository](https://zenodo.org/records/22705862):

1. [Selecting species and identifying core genomes](https://github.com/SushiLab/ocean-subspecies-paper/wiki/01-Species-selection-and-core-genome-definition) - outlines how the set of species to run through the subspecies detection framework were selected, and their core genetic content determined
2. [Subspecies identification framework](https://github.com/SushiLab/ocean-subspecies-paper/wiki/02-Subspecies-identification-framework) - the primary workflow used to identify genomically discrete subspecies from metagenomic Single Nucleotide Variant profiles
3. [Processing subspecies framework output](https://github.com/SushiLab/ocean-subspecies-paper/wiki/03-Processing-subspecies-framework-output) - post-processing steps employed to identify the species for which subspecies structure was confidently resolved.
4. [Mechanisms of diversification](https://github.com/SushiLab/ocean-subspecies-paper/wiki/04-Mechanisms-of-diversification) - the downstream analysis to elucidate the mechanisms through which subspecies are diverging
5. [Single species analysis](https://github.com/SushiLab/ocean-subspecies-paper/wiki/05-Single-species-analysis) - brief description on how to extract information on single species from the produced datasets (which was used to generate the case studies presented in the manuscript)    

## Structure of directories

```
├── pipeline/ --> snakemake pipeline and rules used to conduct the primary analyses of the paper
├── scripts/ --> collection of scripts to process the pipeline output, perform downstream analysis and generate figures and supplementary data
├── tables/ --> scripts to generate the suppl. tables
```

