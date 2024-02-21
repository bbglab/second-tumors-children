# VCF files postprocessing

## Executing the python scripts to process the variant caller files

All files are executed using qmap parallelisation in an HPC environment. Please see documentation here: ```https://github.com/bbglab/qmap```

The notebook ```qmap_process_vcf_files_from_hmf_pipeline_and_sarek.ipynb``` contains all the info to create the processed files.

First, a folder tree is created to store the output files (```./output/```)

Then, qmap files are created for each step.

Each qmap file will execute each corresponding python script (```./python_scripts/```) to process the files from the variant callers, and will save the output in the correspoding output folder.

To execute vep, we need to install and activate vep101.yml conda environment.

Before executing the filter_and_annot.py script, we need to calculate first the CCF and the thresholds for clonality. This is calculated in thisnotebook: ```TMB_and_CCF_analysis.ipynb```

## Examine the processed files

The notebook ```all_coding_alterations_per_sample.ipynb``` inspects and summarises all the somatic and germline alterations per sample.

The notebook ```phylogenetic_trees_tumors.ipynb``` counts the common and unique number of mutations per sample and case.