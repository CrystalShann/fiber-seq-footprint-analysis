#!/bin/bash
#SBATCH --job-name=convert_hic_formats
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/convert_hic_formats_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/convert_hic_formats_%j.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/hicexplorer

# hic to cool
hicConvertFormat \
  -m /project/spott/cshan/annotations/hic/ENCFF216QQM.hic \
  --inputFormat hic \
  --outputFormat cool \
  -o /project/spott/cshan/annotations/hic/ENCFF216QQM.cool

