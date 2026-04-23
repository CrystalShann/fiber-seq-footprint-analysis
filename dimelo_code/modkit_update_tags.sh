#!/bin/bash
#SBATCH --job-name=modkit_convert
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=6
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/modkit_convert_%A_%a.out

module load samtools
module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/dimelo

bam_file="/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples/AL10_bc2178_19130.5mC.6mA.aligned.phased.bam"
converted_bam="/project/spott/cshan/fiber-seq/results/labels/converted/AL10_bc2178_19130.5mC.6mA.aligned.phased.mk.bam"

modkit update-tags \
    --mode implicit \
    --threads 8 \
    "$bam_file" \
    "$converted_bam"

samtools index "$converted_bam"
