#!/bin/bash
#SBATCH --no-requeue
set -e

source "/workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/qmap_files/06_msigact_script_per_sample_synthetic_sensitibity_injections_20240306/execution.env"
if [ -f "/workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/qmap_files/06_msigact_script_per_sample_synthetic_sensitibity_injections_20240306/09.env" ]; then
	source "/workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/qmap_files/06_msigact_script_per_sample_synthetic_sensitibity_injections_20240306/09.env"
fi

. "/home/$USER/miniconda3/etc/profile.d/conda.sh"
conda activate msigact
Rscript ../../../wgs_analysis/tumor_vs_blood/signature_analysis/scripts/run_pipeline.R /workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/synthetic_samples_sensitibity/AZ4614/count_matrix.tsv /workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/synthetic_samples_sensitibity/AZ4614/sigs_SBS31.tsv /workspace/projects/sjd_pediatric_tumors/second-tumors-children/duplex_analysis/signature_analysis/synthetic_samples_sensitibity/AZ4614/SBS31/ 4

