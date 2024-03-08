#!/bin/bash
#SBATCH --no-requeue
set -e

source "/workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/qmap_files/06_msigact_script_per_sample_synthetic_sensitibity_injections_20240308/execution.env"
if [ -f "/workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/qmap_files/06_msigact_script_per_sample_synthetic_sensitibity_injections_20240308/8.env" ]; then
	source "/workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/qmap_files/06_msigact_script_per_sample_synthetic_sensitibity_injections_20240308/8.env"
fi

. "/home/$USER/miniconda3/etc/profile.d/conda.sh"
conda activate msigact
Rscript ../../../wgs_analysis/tumor_vs_blood/signature_analysis/scripts/run_pipeline.R /workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/synthetic_samples_sensitibity/AZ6371/count_matrix.tsv /workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/synthetic_samples_sensitibity/AZ6371/sigs_SBS31.tsv /workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/synthetic_samples_sensitibity/AZ6371/SBS31/ 4

