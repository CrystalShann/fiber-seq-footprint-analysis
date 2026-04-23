#!/bin/bash
#SBATCH --job-name=plot_pgt_overlay
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/plot_pgt_overlay_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/plot_pgt_overlay_%j.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/Jupyter-notebook

TRACKS_INI=${TRACKS_INI:-/project/spott/cshan/fiber-seq/results/TADs_Fiber_MSP/pygenometracks_overlay/chr1.overlay.tracks.ini}
REGION=${REGION:-chr1:100000-150000}
OUT_PNG=${OUT_PNG:-/project/spott/cshan/fiber-seq/results/TADs_Fiber_MSP/pygenometracks_overlay/chr1_2707500_2757500.overlay.png}
OUT_PDF=${OUT_PDF:-/project/spott/cshan/fiber-seq/results/TADs_Fiber_MSP/pygenometracks_overlay/chr1_2707500_2757500.overlay.pdf}

echo "Tracks: $TRACKS_INI"
echo "Region: $REGION"
echo "Output PNG: $OUT_PNG"
echo "Output PDF: $OUT_PDF"

pyGenomeTracks \
    --tracks "$TRACKS_INI" \
    --region "$REGION" \
    --outFileName "$OUT_PNG" \
    --dpi 220

pyGenomeTracks \
    --tracks "$TRACKS_INI" \
    --region "$REGION" \
    --outFileName "$OUT_PDF"

echo "Done. Rendered:"
echo "  $OUT_PNG"
echo "  $OUT_PDF"
