# second-tumors-children
This repository contains the necessary code to reproduce the analysis of "Origin of second tumors in four children and the mutational footprint of chemotherapy in normal tissues"
Sánchez-Guixé M. et al 2024

Original data (CRAM files) can be found at: EGAS50000000167

## Repository organisation:  
```
- wgs_analysis
	|- tumor_vs_blood  
	|	|- mutation_calling
	|	|	|- pipeline5  
	|	|	|- sarek 
	|	|- vcf_postprocessing  
	|	|- ccf_analysis 
	|	|- mutational_profiles  
	|	|- signature_analysis
	|- case3_normal_tissues
	|	|- mutation_calling
	|	|	|- sarek 
	|	|- vcf_postprocessing
	|- case2_clones_analysis
	|	|- signature_analysis
- duplex_analysis
	|- mutational_profiles
	|- signature_analysis
```
1. Whole Genome Sequencing analysis:  
1.1 Tumor vs. Blood:  
    - From CRAM files to variant calling VCF files
        - Pipeline5 Hartwig Medical Foundation pipeline (SAGE, Purple, Haplotypecaller)  
**Supplementary Figure S3C**, **Supplementary Figure S4C**, **Supplementary Figure S9A**, **Supplementary Figure S10B**
        - Sarek pipeline (Mutect2, Strelka2, ASCAT)  
    - From VCF files to postprocessed annotated files  
**Figure 1B**, **Figure 2B**, **Figure 3F**,**Supplementary Figure S3A**, **Supplementary Figure S4A**,**Supplementary Figure S10A**
    - Clonal Cell Fraction plots, definition of clonality threshold
**Figure 1D**, **Figure 2C**  
    - Mutational profiles of clonal and subclonal mutations
**Figure 1C**, **Supplementary Figure S3D**, **Supplementary Figure S4D**,**Supplementary Figure S9B**,**Supplementary Figure S10C**  
    - Signature analysis of clonal and subclonal mutations
**Figure 1D**, **Figure 2C**,**Supplementary Figure S3B**, **Supplementary Figure S4B**,**Supplementary Figure S6A**   

    1.2. Case 3; several normal tissues analysis  
    - From CRAM files to variant calling VCF files  
        - Sarek pipeline (Mutect2, Haplotypecaller)  
    - From VCF files to postprocessed annotated files (Figure 3B, Supplementary Figure 1A-C)  
    - Shared somatic mutations among tissues (Supplementary Figure 5D)  
    
    1.3. Case 3; expanded clones analysis   
    - Mutational signatures analysis

2. **Duplex Sequencing analysis**:  
2.1 Mutational profiles of private mutations  
2.2 Signature analysis of private mutations  


