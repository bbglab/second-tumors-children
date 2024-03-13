# Case 3, analysis with all normal tissues

This analysis is distributed in these 3 parts:

1. ```./mutation_calling/```: mutation calling using sarek, haplotypecaller and mutect tumor only mode
2. ```./vcf_postprocessing```: all the postprocessng steps after the mutation calling, to identify the somatic mutations per tumor and tissues.
3.. ```dPCR```: analysis of the digital PCR assays results.