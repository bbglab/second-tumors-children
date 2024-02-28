# Mutational signature analysis with mSigAct

This folder contains the mutational signatrue analysis with [mSigAct](https://github.com/steverozen/mSigAct).

- ```Analysis_msigact_signatures.ipynb```: notebook with all the code to do prepare the data and to generate all the plots.
- ```./input/```: this folder is created in the notebook, and it will store the tables with the clonals and subclonals mutations from all the samples.
- ```./msigact/```: this folder contains the information to run mSigAct and generate ther results:
    - ```clonal_commands.txt```: commands used to generate the results for the clonal mutations.
    - ```subclonal_commands.txt```: commands used to generate the results for the subclonal mutations.
- ```./data/```: this folder contains necessary data files tu run mSigAct scripts:
    - ```COSMIC_v3.3.1_SBS_GRCh38.txt```: file containing all the probabilities from COSMIC signatures. Downloaded from [COSMIC database](https://cancer.sanger.ac.uk/signatures/downloads/).
- ```msigact.yml```: conda environment to be used when running mSigAct scripts.
- ```./scripts/```: folder containing the necessary scripts to run mSigAct:
    - ```create_matrix_input.py```: script to create the matrix to be used with mSigAct
    - ```create_matrix_sigs.py```: script to create the matrix with the signatures to be tested
    - ```run_pipeline.R```: script to run mSigAct. Note that the paths to the files have to absolute (/root/to/path/...)
- ```./figures_paper/```: folder to store all the plots for the following figures:
    - **Figure 1D**
    - **Figure 2C**
    - **Supplementary Figure S3B**
    - **Supplementary Figure S4B**
    - **Supplementary Figure S6A**
