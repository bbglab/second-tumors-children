# WGS analysis of case 3 tumor samples with all normal tissues and parents' blood samples

Here is described the vcf postprocessing analysis on case 3 using several tissues as normal matching controls, including WGS from parents' blood.

- ```Case3_with_parents_all_tissues.ipynb```: All the step-by-step process is annotated in this notebook
- ```qmap_samtools_depth.ipynb```: notebook to create all the qmaps necessary to obtain the depth per mutation position
- ```./depth/```: folder to save the files and plots with the depth information
- ```./scripts/```: folder with the necessary scripts to run the depth files.
- ```./dictionaries/```: folder to save the necesary dictionaries for the analysis
- ```./output/```: folder to save intermediate and final mutation files from the analysis
