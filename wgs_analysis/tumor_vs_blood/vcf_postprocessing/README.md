# VCF files postprocessing

## Executing the python scripts to process the variant caller files

All files are executed using qmap parallelisation in an HPC environment. Please see documentation [here](https://github.com/bbglab/qmap). Qmaps are executed in the same ```qmap_files/``` folder

The notebook ```qmap_process_vcf_files_from_hmf_pipeline_and_sarek.ipynb``` contains all the info to create the processed files.

First, a folder tree is created to store the output files (```./output/```)

Then, qmap files are created for each step.

Each qmap file will execute each corresponding python script (```./python_scripts/```) to process the files from the variant callers, and will save the output in the correspoding output folder.

To execute vep, we need to install and activate vep101.yml conda environment.

Before executing the filter_and_annot.py script, we need to calculate first the CCF and the thresholds for clonality. This is calculated in thisnotebook: ```TMB_and_CCF_analysis.ipynb```. Also, this notebook creates figures representing the CCF (and are saved in ```./figures_paper/``` folder:

- **Figure 1D** (case1)
- **Figure 2C** (case2)
- **Suplementary Fig. S6A** (case3)

Folder ```./data/``` contains three files:
- ```ensembl_canonical_transcripts.tsv```: file with the canonical transcripts (1 per gene). This file is used in ```process_gridds.py``` script.
- ```genomic_positions_ensembl.txt.gz```: file with starting and ending genomic positions of all genes. This file is used in ```process_gridds.py``` script.

These two files above were downloadad [from ensembl biomart](https://www.ensembl.org/info/data/biomart/index.html) (version 101).

- ```unique_drivers.tsv```: file with the complete set of drivers from Intogen release 2023 (This file can be downloaded from [the IntOGen website](https://www.intogen.org/download/)).

## Examine the processed files

The notebook ```all_coding_alterations_per_sample.ipynb``` inspects and summarises all the somatic and germline alterations per sample. It creates the following tables:

- **Supplementary Table S1**: Table with relevant germline mutations from cancer predisposing genes. Saved in ```./table1_paper/``` as separated tables per case.
- **Supplementary Table S2**: Table with relevant somatic mutations (snv and indels, cnv and sv). Saved in ```./table2_paper/``` as separated tables per case and mutation type.

The notebook ```phylogenetic_trees_tumors.ipynb``` counts the common and unique number of mutations per sample and case. It creates phylogenetic trees per case. The trees shown in the Figures in the paper are created at biorender using the numbers calculated in this notebook. Note that case3 tree from the paper is calculated differenty (see ```case3_normal_tissues```)

- **Figure 1B**
- **Figure 2B**
- **Figure 3F**
- **Supplementary Figure S3A**
- **Supplementary Figure S4A**
- **Supplementary Figure S10A**
