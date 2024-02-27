# WGS, Tumor vs. Blood analyses

There are 4 folders containing all the analyses performed with the WGS Tumor vs. Blood. The order of the analyses is the following:

1. ```mutation_calling```: mutational calling with hmf pipeline and sarek pipeline, starting with the CRAM files from the EGA repository.

2. ```vcf_postprocessing```: postprocessing of all the files from the mutational calling (snv and indels, cnv and sv).

3. ```mutational_profiles```: mutational profile plots of clonal and subclonal snvs from all samples.

4. ```signature_analysis```: mutational signature analysis with mSigAct.