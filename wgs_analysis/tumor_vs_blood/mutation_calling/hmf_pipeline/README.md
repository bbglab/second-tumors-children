# Run Hartwig Medical Foundation pipeline

This pipline is run in google cloud with [Platinum](https://github.com/hartwigmedical/platinum), which contains the [HMF pipeline 5](https://github.com/hartwigmedical/pipeline5).
1. First, revert the downloaded CRAM files into FASTQ files.
2. Second, upload the FASTQ files in gcloud.
3. Then, create an input file with this notebook: ```prepare_input_yaml_file_for_gcloud_platinum_hmf_pipeline.ipynb```
4. And follow the instructions in Platinum repository.
