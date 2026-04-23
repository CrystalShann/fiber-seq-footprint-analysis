#!/bin/bash
#SBATCH --job-name=ctcf_fp
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
#SBATCH --account=pi-spott
#SBATCH --partition=spott
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/ctcf_fp_%A_%a.out

#===============================================================================
# CTCF Footprinting Pipeline — SLURM Array Job
#
# Usage:
#   1. First, count  BAM files:
#      ls /project/spott/1_Shared_projects/LCL_Fiber_seq/final_merged_samples/*.bam | grep -v ".bai" | wc -l
#
#   2. Submit with array range matching number of BAMs (0-indexed):
#      sbatch --array=0-N ctcf_footprinting_slurm.sh
#      (where N = number_of_BAMs - 1)
#
#   Task 0 runs the setup steps (extract motifs, filter, create YAML).
#   All tasks run ft footprint on their assigned BAM file
#===============================================================================



#===============================================================================
# PATHS
#===============================================================================
ANNOTATIONS="/project/spott/cshan/annotations"
RESULTS="/project/spott/cshan/fiber-seq/results"
BAM_DIR="/project/spott/1_Shared_projects/LCL_Fiber_seq/final_merged_samples"
VIERSTRA="/project/spott/kevinluo/Fiber_seq/data/motifs/Vierstra_motif_clustering/hg38.archetype_motifs.v1.0.bed.gz"
FIBERHMM_PEAKS="/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_peaks/refined_peaks/combined_all_chrs_10-30bp_macs3_q0.01_minratio0.5_refined_peaks.bed.gz"
CTCF_CHIPSEQ="${ANNOTATIONS}/CTCF_ENCFF960ZGP.bed"

THREADS=${SLURM_CPUS_PER_TASK:-16}

#===============================================================================
# Build BAM file list
#===============================================================================

# map bam files in the directory to create a bam list and exclude .bai files
mapfile -t BAM_LIST < <(ls ${BAM_DIR}/*.bam | grep -v ".bai")
NUM_BAMS=${#BAM_LIST[@]}

echo "Job array task: ${SLURM_ARRAY_TASK_ID}"
echo "Total BAM files: ${NUM_BAMS}"

if (( SLURM_ARRAY_TASK_ID >= NUM_BAMS )); then
    echo "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} >= number of BAMs (${NUM_BAMS}). Exiting."
    exit 1
fi

# index the bam array to get the corresponding id from slrum task
BAM="${BAM_LIST[${SLURM_ARRAY_TASK_ID}]}"

# extract sample name
SAMPLE=$(basename "$BAM" .bam)
echo "Sample: ${SAMPLE}"
echo "BAM: ${BAM}"

#===============================================================================
# Create directories
#===============================================================================
mkdir -p ${RESULTS}/{motifs,yamls,logs,per_sample/${SAMPLE}/footprints}

#===============================================================================
# SETUP 
#===============================================================================
SETUP_DONE="${RESULTS}/motifs/.setup_complete"

# if the current task = 0, run set up chunk
if (( SLURM_ARRAY_TASK_ID == 0 )); then
    echo ""
    echo "========================================="
    echo "TASK 0: Running setup"
    echo "========================================="

    # Activate bedtools environment
    module load python/miniforge-25.3.0
    source activate /project/spott/cshan/envs/bedtools

    # Extract CTCF motifs from Vierstra 
    echo "Extracting CTCF motifs..."
    
    # extract CTCF froms the files
    zcat "$VIERSTRA" | awk -F'\t' '$4 == "CTCF"' \ 
        > ${RESULTS}/motifs/CTCF_vierstra_14bp.bed
    echo "  Total: $(wc -l < ${RESULTS}/motifs/CTCF_vierstra_14bp.bed)"
    echo "  Widths:"
    awk '{print $3-$2}' ${RESULTS}/motifs/CTCF_vierstra_14bp.bed | sort | uniq -c

    # Create single-module YAML (everything until EOF is passed into cat)
    cat > ${RESULTS}/yamls/CTCF_14bp.yaml << 'EOF'
modules:
  - [0, 14]
EOF
    echo "  YAML:"
    cat ${RESULTS}/yamls/CTCF_14bp.yaml

    # Filter: FiberHMM peaks
    echo "Filtering to FiberHMM 10-30bp peaks..."
    bedtools intersect \
        -a ${RESULTS}/motifs/CTCF_vierstra_14bp.bed \
        -b "$FIBERHMM_PEAKS" \
        -u \
        > ${RESULTS}/motifs/CTCF_14bp_in_fiberhmm_peaks.bed
    echo "  In FiberHMM peaks: $(wc -l < ${RESULTS}/motifs/CTCF_14bp_in_fiberhmm_peaks.bed)"

    # Filter: CTCF ChIP-seq
    echo "Filtering to CTCF ChIP-seq peaks..."
    bedtools intersect \
        -a ${RESULTS}/motifs/CTCF_14bp_in_fiberhmm_peaks.bed \
        -b "$CTCF_CHIPSEQ" \
        -u \
        > ${RESULTS}/motifs/CTCF_footprinting_sites.bed
    echo "  Final sites: $(wc -l < ${RESULTS}/motifs/CTCF_footprinting_sites.bed)"
    echo "  Widths:"
    awk '{print $3-$2}' ${RESULTS}/motifs/CTCF_footprinting_sites.bed | sort | uniq -c

    conda deactivate

    # Signal setup complete
    touch "$SETUP_DONE"
    echo "Setup complete."

else
    echo "Waiting for task 0 setup..."
    while [[ ! -f "$SETUP_DONE" ]]; do
        sleep 10
    done
    echo "Setup detected."
fi

#===============================================================================
# RUN ft footprint
#===============================================================================
echo ""
echo "========================================="
echo "Running ft footprint: ${SAMPLE}"
echo "========================================="

# Activate fibertools environment
module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/fire-env

SITES="${RESULTS}/motifs/CTCF_footprinting_sites.bed"
YAML="${RESULTS}/yamls/CTCF_14bp.yaml"
OUTDIR="${RESULTS}/per_sample/${SAMPLE}/footprints"
FP_OUT="${OUTDIR}/CTCF_footprint_results.tsv"

if [[ -f "$FP_OUT" ]]; then
    echo "SKIP: output already exists"
else
    echo "Sites: $(wc -l < "$SITES")"

    ft footprint \
        "$BAM" \
        --bed "$SITES" \
        --yaml "$YAML" \
        -t "$THREADS" \
        > "$FP_OUT" \
        2> "${RESULTS}/logs/ft_footprint_${SAMPLE}.log"

    echo "Output: $(tail -n +2 "$FP_OUT" | wc -l) sites"
fi

#===============================================================================
# Filter to >= 20 spanning fibers
#===============================================================================
echo ""
echo "Filtering to >= 20 fibers..."

FILTERED="${OUTDIR}/CTCF_footprint_ge20fibers.tsv"
head -1 "$FP_OUT" > "$FILTERED"
awk -F'\t' 'NR>1 && $5 >= 20' "$FP_OUT" >> "$FILTERED"

echo "  All sites:   $(tail -n +2 "$FP_OUT" | wc -l)"
echo "  >=20 fibers: $(tail -n +2 "$FILTERED" | wc -l)"

#===============================================================================
# Parse footprint codes
#===============================================================================
# 1 module --> 2 bits:
#   bit 0 (value 1) = spanning MSP (accessible)
#   bit 1 (value 2) = module footprinted (bound)
#
#   code 0 -> inaccessible (nucleosome)
#   code 1 -> accessible, unbound
#   code 3 -> accessible, bound
#===============================================================================
echo ""
echo "Parsing footprint codes..."

PARSED="${OUTDIR}/CTCF_parsed_states.tsv"

python3 - "$FILTERED" "$PARSED" "$SAMPLE" << 'PYEOF'
import sys

infile, outfile, sample = sys.argv[1], sys.argv[2], sys.argv[3]

total_fibers = 0
bound = 0
accessible_unbound = 0
inaccessible = 0

with open(infile) as fin, open(outfile, 'w') as fout:
    header = fin.readline().strip().split('\t')
    fp_idx = header.index('footprint_codes')

    fout.write("chrom\tstart\tend\tstrand\tn_fibers\t"
               "n_bound\tn_accessible_unbound\tn_inaccessible\t"
               "pct_bound\tpct_accessible_unbound\tpct_inaccessible\n")

    for line in fin:
        f = line.strip().split('\t')
        chrom, start, end, strand = f[0], f[1], f[2], f[3]

        codes_str = f[fp_idx]
        if codes_str in ('', '.'):
            continue

        codes = [int(x) for x in codes_str.split(',')]
        sb, su, si = 0, 0, 0

        for code in codes:
            has_msp = (code & 1) > 0
            has_fp  = (code & 2) > 0

            if not has_msp:
                si += 1
                inaccessible += 1
            elif has_fp:
                sb += 1
                bound += 1
            else:
                su += 1
                accessible_unbound += 1

        n = len(codes)
        total_fibers += n
        fout.write(f"{chrom}\t{start}\t{end}\t{strand}\t{n}\t"
                   f"{sb}\t{su}\t{si}\t"
                   f"{sb/n*100:.1f}\t{su/n*100:.1f}\t{si/n*100:.1f}\n")

print(f"\n=== {sample} SUMMARY ===", file=sys.stderr)
print(f"  Total fibers: {total_fibers}", file=sys.stderr)
print(f"    Bound:              {bound} ({bound/max(total_fibers,1)*100:.1f}%)", file=sys.stderr)
print(f"    Accessible unbound: {accessible_unbound} ({accessible_unbound/max(total_fibers,1)*100:.1f}%)", file=sys.stderr)
print(f"    Inaccessible:       {inaccessible} ({inaccessible/max(total_fibers,1)*100:.1f}%)", file=sys.stderr)
PYEOF

echo ""
echo "========================================="
echo "Task ${SLURM_ARRAY_TASK_ID} (${SAMPLE}) COMPLETE"
echo "========================================="
echo "  ${FP_OUT}"
echo "  ${FILTERED}"
echo "  ${PARSED}"
