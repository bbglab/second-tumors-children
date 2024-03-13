# VCF post-processing of the clones derived from the MRT from case 3

This analysis is perfomed as follows:

1. ```qmap_process_vcf_files_from_hmf_pipeline_and_sarek.ipynb```: Notebook to create the post-processed mafs frmo the VCF files from strelka.
2. ```strelka_process.qmap```: qmap file to process vc files from strelka.
3. ```./output/```: folder to save the processed files from strelka.
4. ```TMB_and_CCF_analysis_case3_t2_clones.ipynb```: Notebook to get the tumor mutational burden and clonal cell fraction (ccf) of the MRT clones.
5. ```ccf_thresholds_case3_t2_clones.json```: the calculated ccf thresholds for the clones.