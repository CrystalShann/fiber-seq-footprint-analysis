#!/bin/bash
#SBATCH --job-name=msp_block
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/msp_block_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/msp_block_%j.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/Jupyter-notebook

python /project/spott/cshan/fiber-seq/code/TADs_Fiber_MSP_code/get_msp_block.py
