# second-tumors-children
This repository contains the necessary code to reproduce the analysis of "Origins of second tumors in children and mutational footprint of chemotherapy in normal tissues"
Sánchez-Guixé M. et al 2024

Original data (CRAM files) can be found at: EGAS50000000167

## Repository organisation:  
```
- wgs_analysis
	|- tumor_vs_blood  
	|	|- mutation_calling
	|	|	|- hmf_pipeline  
	|	|	|- sarek 
	|	|- vcf_postprocessing  
	|	|- ccf_analysis 
	|	|- mutational_profiles  
	|	|- signature_analysis
	|- case3_normal_tissues
	|	|- mutation_calling
	|	|	|- sarek 
	|	|- vcf_postprocessig
	|	|- dPCR
	|- case3_clones_analysis
	|	|- signature_analysis
- duplex_analysis
	|- mutational_profiles
	|- signature_analysis
```
## Folders descriptions and allocation of Figures:  
1. Whole Genome Sequencing analysis ```wgs_analysis```:  
1.1 Tumor vs. Blood ```tumor_vs_blood```:  
    - From CRAM files to variant calling VCF files ```mutaton_calling```:  
        - Hartwig Medical Foundation pipeline (SAGE, Purple, Haplotypecaller) ```hmf_pipeline```:  
**Supplementary Figure S3C**, **Supplementary Figure S4C**, **Supplementary Figure S9A**, **Supplementary Figure S10B**  
        - Sarek pipeline (Mutect2, Strelka2, ASCAT) ```sarek```:  
    - From VCF files to postprocessed annotated files ```vcf_postprocessing```:  
**Figure 1B**, **Figure 2B**, **Figure 3F**,**Supplementary Figure S3A**, **Supplementary Figure S4A**,**Supplementary Figure S10A**  
    - Clonal Cell Fraction plots, definition of clonality threshold ```ccf_analysis```:  
**Figure 1D**, **Figure 2C**  
    - Mutational profiles of clonal and subclonal mutations ```mutational_profiles```:  
**Figure 1C**, **Supplementary Figure S3D**, **Supplementary Figure S4D**,**Supplementary Figure S9B**,**Supplementary Figure S10C**  
    - Signature analysis of clonal and subclonal mutations ```signature_analysis```:  
**Figure 1D**, **Figure 2C**, **Supplementary Figure S3B**, **Supplementary Figure S4B**,**Supplementary Figure S6A**   

    1.2. Case 3; several normal tissues analysis ```case3_normal_tissues```:  
    - From CRAM files to variant calling VCF files ```mutation_calling```:  
        - Sarek pipeline (Mutect2, Haplotypecaller) ```sarek```:  
    - From VCF files to postprocessed annotated files ```vcf_postprocessing```:  
**Figure 3B**, **Supplementary Figure S5B,C,D**
    - Digital PCR results ```dPCR```:  
**Supplementary Figure S8A**

    1.3. Case 3; expanded clones analysis ```case3_clones_analysis```:   
    - Mutational signatures analysis ```signature_analysis```:  
**Figure 3C**, **Supplementary Figure S6C**  

2. Duplex Sequencing analysis ```duplex_sequencing```:
    - Mutational profiles of private mutations ```mutational_profiles```:  
**Supplementary Figure S11**  
    - Signature analysis of private mutations ```signature_analysis```:  
**Supplementary Figure S6B**, **Supplementary Figure S12**, **Supplementary Figure S13**, **Supplementary Figure S14**  


