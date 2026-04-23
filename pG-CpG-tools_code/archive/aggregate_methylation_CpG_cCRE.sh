#!/bin/bash
#SBATCH --job-name=CpG_met
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/CpG_met_%j.err


export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/code/pG-CpG-tools_code/aggregate_methylation_CpG_cCRE.R