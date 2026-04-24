#!/bin/bash
#SBATCH --job-name=fp_per_read_tss
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_tss_%A_%a.err
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_tss_%A_%a.out

export LANG=C
export LC_ALL=C

module load R/4.4.1

sizes=(10-30 30-45 45-60 60-80 140-160)

Rscript /project/spott/cshan/fiber-seq/code/footprints_per_read_TSS/footprints_TSS_profiles_per_read.r \
  sample=AL10_bc2178_19130 \
  mode=size_chunk \
  fp_size=${sizes[$SLURM_ARRAY_TASK_ID]} \
  chrom=chr1
