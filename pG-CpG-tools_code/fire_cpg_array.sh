#!/bin/bash
#SBATCH --job-name=fire_cpg
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fire_cpg_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fire_cpg_%A_%a.err

set -euo pipefail

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/pb-CpG-tools

FIRE_dir="/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
ref_genome="/project/spott/reference/human/GRCh38/hg38.fa"
cpg_output_base="/project/spott/cshan/fiber-seq/results/fire_CpG"
log_dir="/project/spott/cshan/fiber-seq/results/logs"
threads="${SLURM_CPUS_PER_TASK:-2}"

mkdir -p "${cpg_output_base}" "${log_dir}"

build_pending_sample_list() {
  local sample_dir sample pending_file run_tag
  run_tag="$(date +%Y%m%d_%H%M%S)_$$"
  pending_file="${log_dir}/fire_cpg_pending_samples.${run_tag}.txt"
  : > "${pending_file}"

  while IFS= read -r sample_dir; do
    sample="$(basename "${sample_dir}")"

    if [[ -f "${cpg_output_base}/${sample}/${sample}_CPG.combined.bw" ]]; then
      echo "Already complete: ${sample}" >&2
      continue
    fi

    if ! find "${sample_dir}" -maxdepth 1 -type f \
      \( -name "*-fire-v0.1-filtered.cram" -o -name "*-fire-v0.1.1-filtered.cram" \) \
      | grep -q .; then
      echo "Skipping ${sample}: no matching CRAM found" >&2
      continue
    fi

    printf '%s\n' "${sample_dir}" >> "${pending_file}"
  done < <(find "${FIRE_dir}" -mindepth 1 -maxdepth 1 -type d | sort)

  printf '%s\n' "${pending_file}"
}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  pending_samples_file="$(build_pending_sample_list)"

  mapfile -t pending_sample_dirs < "${pending_samples_file}"
  pending_count="${#pending_sample_dirs[@]}"

  echo "Found ${pending_count} sample(s) still needing CpG output."

  if [[ "${pending_count}" -eq 0 ]]; then
    echo "Nothing to do."
    exit 0
  fi

  array_end=$((pending_count - 1))

  echo "Submitting array job for pending samples: 0-${array_end}"
  sbatch \
    --array="0-${array_end}" \
    --export=ALL,FIRE_CPG_PENDING_SAMPLES_FILE="${pending_samples_file}" \
    "$0"
  exit 0
fi

if [[ -z "${FIRE_CPG_PENDING_SAMPLES_FILE:-}" ]]; then
  echo "ERROR: FIRE_CPG_PENDING_SAMPLES_FILE is not set for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

if [[ ! -f "${FIRE_CPG_PENDING_SAMPLES_FILE}" ]]; then
  echo "ERROR: Pending sample list not found: ${FIRE_CPG_PENDING_SAMPLES_FILE}"
  exit 1
fi

mapfile -t sample_dirs < "${FIRE_CPG_PENDING_SAMPLES_FILE}"

if [[ "${#sample_dirs[@]}" -eq 0 ]]; then
  echo "ERROR: Pending sample list is empty."
  exit 1
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${#sample_dirs[@]}" ]]; then
  echo "No sample for task ${SLURM_ARRAY_TASK_ID}; total pending samples=${#sample_dirs[@]}. Exiting."
  exit 0
fi

sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]}"
sample="$(basename "${sample_dir}")"

echo "========================================="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${sample}"
echo "Sample directory: ${sample_dir}"
echo "========================================="

cram_file="$(find "${sample_dir}" -maxdepth 1 -type f \
  \( -name "*-fire-v0.1-filtered.cram" -o -name "*-fire-v0.1.1-filtered.cram" \) \
  | sort \
  | head -n 1)"

if [[ -z "${cram_file}" ]]; then
  echo "ERROR: No CRAM file found in ${sample_dir}"
  exit 1
fi

echo "Found CRAM: ${cram_file}"

cpg_output_dir="${cpg_output_base}/${sample}"
mkdir -p "${cpg_output_dir}"

if [[ -f "${cpg_output_dir}/${sample}_CPG.combined.bw" ]]; then
  echo "CpG combined bigWig already exists; skipping ${sample}"
  exit 0
fi

echo "Running aligned_bam_to_cpg_scores..."
aligned_bam_to_cpg_scores \
  --bam "${cram_file}" \
  --output-prefix "${cpg_output_dir}/${sample}_CPG" \
  --ref "${ref_genome}" \
  --threads "${threads}"

echo "CpG scores complete for ${sample}"
