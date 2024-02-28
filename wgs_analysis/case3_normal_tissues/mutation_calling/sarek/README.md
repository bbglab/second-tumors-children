# Run sarek

To run sarek, first install and activate the conda environment ```mutcall.yml```\

Reference Genome files to run sarek can be downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/).

Notice that the Reference Genome used is without ALT contigs (this is the one used in hmf pipeline, platinum in gcloud). Explanation [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use).

- ```commands.txt```: The commands used to run sarek with mutect (tumor only) and haplotypecaller.
- ```Input_for_sarek_normal_tissues_case3.ipynb```: Notebook that explains how to prepare the input file.
