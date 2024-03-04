# WGS analysis of case 3 tumor samples with all normal tissues and parents' blood samples

Here is described the vcf postprocessing analysis on case 3 using several tissues as normal matching controls, including WGS from parents' blood.

- ```Case3_with_parents_all_tissues.ipynb```: notebook with all the step-by-step process with all te filterings to obtain the final list of mutations. 
- ```Case3_final_list_of_muts_heatmap_plot```: notebook to create the plots: Clonality thresold, depth threshold and heatmap.
- ```qmap_samtools_depth_vep_and_artifacts.ipynb```: notebook to create all the qmaps necessary to obtain the depth per mutation position, to run vep and to run the betabinomial test for the artifacts.
- ```calculate_sequencing_artifacts_with_betabinomial.ipynb```: notebook to annotate the artifacts, using a cohort of blood samples from 26 unrelated patients.
- ```./depth/```: folder to save the files and plots with the depth information
- ```./scripts/```: folder with the necessary scripts to run the analysis.
- ```./dictionaries/```: folder to save the necesary dictionaries for the analysis
- ```./output/```: folder to save intermediate and final mutation files from the analysis
- ```./qmap_files/```: folder to save all the qmap files (created in ```qmap_samtools_depth_vep_and_artifacts.ipynb```)
- ```./vep/```: folder with the input, output and processed files from vep.
- ```./rescued_muts/```: folder to annotate mutations frmo bams.
- ```./betabinom/```: folder to calculate the betabinomial test for the artifacts filter.
- ```./figures_paper/```: folder to save the figures for the paper:
  - **Figure 3B**
  - **Suypplementary Figure S5A**: Number of mutations in the step-by-step filterings. Note that the final plot is builded in biorender, this plots just shows the numbers.
  - **Supplementary Figure S5B**: Depth of chromosome 1 in blood sample,  showing the thresholds of low and high depth
  - **Supplementary Figure S5C**: Clonality thresholds of NB and MRT.
