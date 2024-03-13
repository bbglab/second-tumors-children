# second-tumors-children
This repository contains the necessary code to reproduce the analysis of "Origins of Second Malignancies in Children and Mutational Footprint of Chemotherapy in Normal Tissues"
Sánchez-Guixé M. et al 2024 (doi: 10.1158/2159-8290.CD-23-1186).

## Original data

Original data (CRAM files) can be found at: EGAS50000000167 (Datasets: EGAD50000000237, EGAD50000000238, EGAD50000000239, EGAD50000000240).  

## Conda environments

```process_vc.yml```: This is the main conda environment used in this project. This is used in all jupyter notebooks as kernel.  
Occasionaly, other conda environments are used in some analysis. This is annotated accordingly.  

## Repository organisation:  
```
- wgs_analysis
	|- tumor_vs_blood  
	|	|- mutation_calling
	|	|	|- hmf_pipeline  
	|	|	|- sarek 
	|	|- vcf_postprocessing   
	|	|- mutational_profiles  
	|	|- signature_analysis
	|- case3_normal_tissues
	|	|- mutation_calling
	|	|	|- sarek 
	|	|- vcf_postprocessig
	|	|- dPCR
	|- case3_clones_analysis
	|	|- mutation_calling
        |       |       |- sarek
        |       |- vcf_postprocessing
        |       |- signature_analysis
- duplex_analysis
	|- mutational_profiles
	|- signature_analysis
```
## Folders descriptions and allocation of Figures:  
1. ```wgs_analysis```, Whole Genome Sequencing analysis:  
1.1 ```tumor_vs_blood```, Tumor vs. Blood:  
    - ```mutaton_calling```, From CRAM files to variant calling VCF files:  
        - ```hmf_pipeline```, Hartwig Medical Foundation pipeline (SAGE, Purple, Haplotypecaller):
        	- 	**Supplementary Figure S3C**
        	- 	**Supplementary Figure S4C**
        	- 	**Supplementary Figure S9A**
      		-  	**Supplementary Figure S10B**  
        - ```sarek```, Sarek pipeline (Mutect2, Strelka2, ASCAT):  
    - ```vcf_postprocessing```, From VCF files to postprocessed annotated files, Clonal Cell Fraction and trees:
    	-	**Figure 1B**
    	-	**Figure 2B**
    	-	**Figure 3F**
    	-	**Supplementary Figure S3A**
    	-	**Supplementary Figure S4A**
    	-	**Supplementary Figure S10A**
    	-	**Figure 1D**
    	-	**Figure 2C**
        - 	**Supplementary Figure S6A**	
     	-	**Supplementary Table S1**
      	-	**Supplementary Table S2**
    - ```mutational_profiles```, Mutational profiles of clonal and subclonal mutations:
    	- 	**Figure 1C**
     	- 	**Supplementary Figure S3D**
      	- 	**Supplementary Figure S4D**
      	- 	**Supplementary Figure S9B**
      	- 	**Supplementary Figure S10C**  
    - ```signature_analysis```, Signature analysis of clonal and subclonal mutations:
    	- 	**Figure 1D**
     	- 	**Figure 2C**
      	- 	**Supplementary Figure S3B**
      	- 	**Supplementary Figure S4B**
      	- 	**Supplementary Figure S6A**   

    1.2. ```case3_normal_tissues```, Case 3; several normal tissues analysis:  
    - ```mutation_calling```, From CRAM files to variant calling VCF files:  
        - ```sarek```, Sarek pipeline (Mutect2, Haplotypecaller):  
    - ```vcf_postprocessing```, From VCF files to postprocessed annotated files:
    	- 	**Figure 3B**
     	- 	**Supplementary Figure S5A,B,C,D**
    - ```dPCR```, Digital PCR results:
    	- 	**Supplementary Figure S8A**

    1.3. ```case3_clones_analysis```, Case 3; expanded clones analysis:
    - ```mutation_calling```, From CRAM files to variant calling VCF files:
    	- ```sarek```, Sarek pipeline (Strelka):
    - ```vcf_postprocessing```, From VCF files to postprocessed annotated files, calculate clonality threshold:
    - ```signature_analysis```, Mutational signatures analysis:
    	- 	**Figure 3C**
     	- 	**Supplementary Figure S6C**  

3. ```duplex_sequencing```, Duplex Sequencing analysis:
    - ```mutational_profiles```, Mutational profiles of private mutations:
    	- 	**Supplementary Figure S11**  
    - ```signature_analysis```, Signature analysis of private mutations:
    	- 	**Figure 4**
    	- 	**Supplementary Figure S6B**
     	- 	**Supplementary Figure S12**
      	- 	**Supplementary Figure S13**
      	- 	**Supplementary Figure S14**  


