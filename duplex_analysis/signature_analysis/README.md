# Mutational signatures analysis on duplex sequencing samples

Analysis to identify which samples have significant amount of SBS31 mutations, and in which range (confidence intervals).

To identify the list of mutational signatures for the fitting for each sample, we used a methodology called "combinatorial refitting".

Then, we evaluate with msigact the presence of SBS31.

To be sure that the method is specific enough, we evaluate the presence of SBS31 in a set of synthetic samples containing the proportional amounts of the non-SBS31 signatures. The sample will have a significant (specific) number of SBS31 mutations when no synthetic sample reaches the number of SBS31 detected mutations in the sample.

To calculate a range of sensitivity of the method, we evaluate the number of SBS31 mutations in a set of synthetic samples with increasing number of SBS31 injected mutations and proportional number of mutations for the non-SBS31 mutational signatures. The range in which the method detects the same number of SBS31 mutations as in the sample, is the confidence interval. We do the same analysis for the aging mutations (SBS1, SBS5 and SBS40).

- ```01_selection_of_signatures_with_combinatorial_refitting.ipynb```: notebook to select the list of mutational signatures to fit in each sample, using combinatorial refitting.
- ```02_signature_analysis_with_mSigAct.ipynb```: notebook to analyse the presence of SBS31 signature with mSigAct, once we have already selected the signatures per sample.
- ```03_qmap_msigact.ipynb```: notebook to create the qmap files necessary to run mSigAct per sample (with different set of signatures per sample) and the specificity and sensitibity tests.
- ```04_test_sbs31_specificity_and_sensibility.ipynb```: notebook to perform specificity test (significant presence of SBS31), and sensitibity test (confidence interval).
- ```./data/```: folder containing the triplet counts for the sequencing panel used in duplex sequencing (Twinstrand Mutagenesis).
- ```./msigact/```: folder with the analysis of mSigAct in each sample. The commands are in the qmap files.
- ```./synthetic_samples_specificity/```: folder with the created synthetic samples based on the proportion of signatures for each sample without SBS31, and the corresponding mSigAct analysis.
- ```./synthetic_samples_sensitivity/```: folder with the created synthetic samples with increasing injections of SBS31 mutations (and proportionally the rest of mutational signatures), and the corresponding mSigAct analysis.
- ```./qmap_files/```: qmap files to run mSigAct per sample and in the synthetic samples for the specificity test and the sensitivity test.
- ```./confidence_intervals/```: json files with the calculated confidence intervals for SBS31 and aging signatures.
- ```./figures_paper/```: folder containing the figures for the paper. Note that the specificity and sensitivity tests are not exactly the same figures due to the randomisation of the syntehtic samples everytime we execute the code.
  - **Figure 4**: Number of mutations (per mappable genome) and mutation rate for each tissue sequenced with duplex sequecing (Case2, Case3 and unrelated donors) and WGS (Case1 tAML sample, clonal mutations).
  - **Supplementary Figure S6B**: Number of mutations (per mappable genome) of duplex sequencing MRT sample from Case 3.
  - **Supplementary Figure S12**: Relative exposure of the different signatures, testing the fitting with either SBS31 or SBS35. Note that the asterisk added in the paper figure denote significance after specificity test.
  - **Supplementary Figure S13**: Specificity test for the presence of SBS31 signature in the samples. Note that these plots may not be exactly as the ones in the figure of the paper due to the randomness in the generation of the synthetic samples. 
  - **Supplementary Figure S14**: Sensitivity test for the confidecne interval of SBS31 exposure in the samples. Note that these plots may not be exactly as the ones in the figure of the paper due to the randomness in the generation of the synthetic samples.
  