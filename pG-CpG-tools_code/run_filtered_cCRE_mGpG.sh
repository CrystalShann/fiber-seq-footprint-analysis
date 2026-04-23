#!/bin/bash
#SBATCH --job-name=filtered_mCpG
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/filtered_mCpG_%j.err
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/filtered_mCpG_%j.out


export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/code/pG-CpG-tools_code/run_filtered_accessible_cCRE_CpG_methylation.R