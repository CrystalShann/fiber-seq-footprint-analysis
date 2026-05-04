#!/bin/bash
#SBATCH --job-name=convert_hic_formats
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/convert_hic_formats_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/convert_hic_formats_%j.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/hictk
# source activate /project/spott/cshan/envs/hicexplorer 
## hicexplort does not work with v9 .hic file


TMPDIR=/scratch/midway3/$USER/hictk_tmp_$SLURM_JOB_ID
mkdir -p "$TMPDIR"
export TMPDIR

echo "Using TMPDIR: $TMPDIR"

# -------------------------
# Input / output
# -------------------------
HIC_FILE="/project/spott/cshan/annotations/hic/ENCFF216QQM.hic"
OUT_FILE="/project/spott/cshan/annotations/hic/ENCFF216QQM.mcool"

# -------------------------
# Run conversion
# -------------------------
echo "Starting conversion..."

hictk convert \
  "$HIC_FILE" \
  "$OUT_FILE" \
  --tmpdir "$TMPDIR"

echo "Conversion finished."

# -------------------------
# Clean up temp 
# -------------------------
rm -rf "$TMPDIR"

