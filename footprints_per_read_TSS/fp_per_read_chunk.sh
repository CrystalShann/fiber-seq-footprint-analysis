#!/bin/bash
#SBATCH --job-name=fp_per_read_chunk
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_chunk_%j.err
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_chunk_%j.out

export LANG=C
export LC_ALL=C

module load R/4.4.1

# FP_SIZE must be passed via --export=FP_SIZE=<size> when submitting
# e.g. sbatch --export=FP_SIZE=10-30 fp_per_read_chunk.sh

Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
  sample=AL10_bc2178_19130 \
  mode=size_chunk \
  fp_size=${FP_SIZE}
