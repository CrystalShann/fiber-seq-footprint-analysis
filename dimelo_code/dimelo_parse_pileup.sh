#!/bin/bash
#SBATCH --job-name=dimelo_parse_pileup
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=6
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/dimelo_parse_pileup_%j.out

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/dimelo

python - <<'PY'
import os
from dimelo import parse_bam

# define parameters for parsing
ref_genome_file = '/project/spott/reference/human/GRCh38/genome/hg38.fa'
bam_file = '/project/spott/cshan/fiber-seq/results/modkit_converted/AL10_bc2178_19130.5mC.6mA.aligned.phased.mk.bam'
output_dir = '/project/spott/cshan/fiber-seq/results/modkit_converted'

# parsing can optionally specify mod codes.
# default is Y/a for adenine and m/Z for cytosine methylation
motifs = ['A,0', 'CG,0']

parse_prefix = 'parsed'
met_thresh = 190
ncore=1

print(f"Running dimelo with cores={ncore}")
print("Running dimelo pileup...")

pileup_file, pileup_regions = parse_bam.pileup(
    input_file=bam_file,
    output_name=parse_prefix,
    ref_genome=ref_genome_file,
    output_directory=output_dir,
    # regions=regions_bed,
    motifs=motifs,
    thresh=met_thresh,
    # window_size=window_size,
    cores=ncore,
    quiet=True,
    cleanup=cleanup_pileup,
    override_checks=True,
    log=True,
)

print("pileup_file:", pileup_file)
print("Expected output directory:", f"{output_dir}/{parse_prefix}")
PY
