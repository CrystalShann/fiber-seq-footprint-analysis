#!/bin/bash
#SBATCH --job-name=cluster_fp_extract
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/cluster_fp_extract_%j.err
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/cluster_fp_extract_%j.out

set -euo pipefail

export LANG=C
export LC_ALL=C

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/Jupyter-notebook

python /project/spott/cshan/fiber-seq/code/PolII_footprints_code/cluster_reads_code/cluster_reads_by_FP_features_extract.py
