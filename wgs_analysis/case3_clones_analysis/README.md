# WGS analysis of 2 clones derived from the MRT (tumor2) cell line from case 3

This analysis is organised as follows:

1. ```./mutation_calling/```: mutation calling instructions using sarek to run strelka.
2. ```./vcf_postprocessing/```: postprocessing after mutation calling, to get the clonal cell fraction and the clonality thresholds of the clones.
3. ```./signature_analysis/```: mSigAct analysis to evaluate the presence of SBS31 in the clones derived from the MRT tumor.