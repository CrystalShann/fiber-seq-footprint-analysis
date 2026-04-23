#!/bin/bash
#SBATCH --job-name=CAGE_all_chr
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/CAGE_all_chr_%j.err
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/CAGE_all_chr_%j.out


# Fix locale - use C which is always available on all nodes
export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/CAGE_TSS_gene_aggregation_all_chr.r