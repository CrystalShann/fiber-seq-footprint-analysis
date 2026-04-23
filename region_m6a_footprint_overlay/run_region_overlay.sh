#!/bin/bash
#SBATCH --job-name=region_m6a
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/region_m6a_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/region_m6a_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=4

SAMPLE=""
REGION=""
CHROM=""
START=""
END=""
OUTPUT_ROOT="/project/spott/cshan/fiber-seq/results/region_m6a_footprint_overlay"
METADATA="/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
MAX_READS="250"
MOD_PROB_THRESHOLD="0.9"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --region) REGION="$2"; shift 2 ;;
    --chrom) CHROM="$2"; shift 2 ;;
    --start) START="$2"; shift 2 ;;
    --end) END="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --max-reads) MAX_READS="$2"; shift 2 ;;
    --mod-prob-threshold) MOD_PROB_THRESHOLD="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done


if [[ -z "$SAMPLE" ]]; then
  echo "Provide --sample." >&2
  exit 1
fi

if [[ -n "$REGION" ]]; then
  CHROM="${REGION%%:*}"
  COORDS="${REGION#*:}"
  START="${COORDS%-*}"
  END="${COORDS#*-}"
else
  if [[ -z "$CHROM" || -z "$START" || -z "$END" ]]; then
    echo "Provide --region or --chrom/--start/--end." >&2
    exit 1
  fi
  REGION="${CHROM}:${START}-${END}"
fi

if [[ "$CHROM" != chr* ]]; then
  CHROM="chr${CHROM}"
  REGION="${CHROM}:${START}-${END}"
fi

REGION_SLUG="${CHROM}_${START}_${END}"
REGION_DIR="${OUTPUT_ROOT}/${SAMPLE}/${CHROM}/${REGION_SLUG}"
R_SCRIPT="/project/spott/cshan/fiber-seq/code/region_m6a_footprint_overlay/extract_region_m6a_summary.R"
PY_SCRIPT="/project/spott/cshan/fiber-seq/code/region_m6a_footprint_overlay/build_region_overlay_tracks.py"
SUMMARY_BW="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.m6A_fraction.bw"
SUMMARY_BEDGRAPH="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.m6A_fraction.bedgraph"
SUMMARY_5MC_BW="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.5mC_fraction.bw"
SUMMARY_5MC_BEDGRAPH="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.5mC_fraction.bedgraph"
TRACKS_INI="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.tracks.ini"
OUT_PNG="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.overlay.png"
OUT_PDF="${REGION_DIR}/${SAMPLE}.${REGION_SLUG}.overlay.pdf"

mkdir -p "${REGION_DIR}"

module load R
Rscript "${R_SCRIPT}" \
  --sample "${SAMPLE}" \
  --region "${REGION}" \
  --metadata "${METADATA}" \
  --output_root "${OUTPUT_ROOT}" \
  --mod_prob_threshold "${MOD_PROB_THRESHOLD}"

M6A_SUMMARY_TRACK="${SUMMARY_BEDGRAPH}"
if [[ -f "${SUMMARY_BW}" ]]; then
  M6A_SUMMARY_TRACK="${SUMMARY_BW}"
fi

MC5_SUMMARY_TRACK="${SUMMARY_5MC_BEDGRAPH}"
if [[ -f "${SUMMARY_5MC_BW}" ]]; then
  MC5_SUMMARY_TRACK="${SUMMARY_5MC_BW}"
fi

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/Jupyter-notebook

python "${PY_SCRIPT}" \
  --sample "${SAMPLE}" \
  --region "${REGION}" \
  --metadata "${METADATA}" \
  --output-root "${OUTPUT_ROOT}" \
  --m6a-summary-track "${M6A_SUMMARY_TRACK}" \
  --5mc-summary-track "${MC5_SUMMARY_TRACK}" \
  --max-reads "${MAX_READS}" >/dev/null

pyGenomeTracks \
  --tracks "${TRACKS_INI}" \
  --region "${REGION}" \
  --outFileName "${OUT_PNG}" \
  --dpi 220

pyGenomeTracks \
  --tracks "${TRACKS_INI}" \
  --region "${REGION}" \
  --outFileName "${OUT_PDF}"

echo "Tracks ini: ${TRACKS_INI}"
echo "PNG: ${OUT_PNG}"
echo "PDF: ${OUT_PDF}"
