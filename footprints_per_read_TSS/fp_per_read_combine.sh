#!/bin/bash
#SBATCH --job-name=fp_per_read_combine
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_combine_%j.err
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_combine_%j.out

export LANG=C
export LC_ALL=C

module load R/4.4.1

Rscript /project/spott/cshan/fiber-seq/code/footprints_per_read_TSS/footprints_TSS_profiles_per_read.r \
  sample=AL10_bc2178_19130 \
  mode=combine
