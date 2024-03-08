# Mutational profiles of duplex sequencing samples
- ```Mutational_profiles.ipynb```: notebook to generate the mutational profiles.
- ```./data/```: json file with the triplet counts in the mappable regions in the genome.
- ```mu_sample_sheet.csv```: csv file used as input in dnanexus Twinstrand Mutagenesis pipeline. First, cram files have to be converted into fastq files.
Then, follow the instructions at the Twinstrand Mutagenesis App at dnanexus. You will need to log in in dnanexus ([more info here](https://twinstrandbio.com/wp-content/uploads/TwinStrand-Mutagenesis-Brochure-V2.pdf))
- ```./output_dnanexus/```: variant calling files from dnanexus Twinstrand Mutagenesis pipeline (consensus.variant-calls.genome.mut and consensus.variant-calls.mut).
- ```./scripts/```: scripts to run the profiles
- ```./figures_paper/```: mutational profile plots.
  - **Supplementary Figure S11**
