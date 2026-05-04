#!/bin/bash
set -euo pipefail

# Extract FIRE calls for one region. This is a local copy of the helper used by
# coaccess_fire_cres_high_low.Rmd, with the project FIRE output version set to v0.1.

fire_name=$1
bam_file=$2
fire_dir=$3
region=$4
plot_name=$5
out_dir=$6

fire_version="v0.1"

module load python/anaconda-2022.05

if [[ $region =~ ^(chr[0-9XYM]+):([0-9]+)-([0-9]+)$ ]]; then
    chrom="${BASH_REMATCH[1]}"
    start_1based="${BASH_REMATCH[2]}"
    end_1based="${BASH_REMATCH[3]}"
else
    echo "Error: Invalid region format. Expected format is 'chrN:start-end'."
    exit 1
fi

echo "Processing sample: $fire_name"
echo "BAM file: $bam_file"
echo "FIRE directory: $fire_dir"
echo "Region: $region"
echo "Output directory: $out_dir"
echo "Chromosome: $chrom"
echo "Plot name: $plot_name"

mkdir -p "$out_dir/$plot_name/parsed"

region_bed="$out_dir/$plot_name/region.bed"
start_0based=$((start_1based - 1))
printf "%s\t%s\t%s\n" "$chrom" "$start_0based" "$end_1based" > "$region_bed"

echo "Subsetting BAM file to region $region..."
module load samtools
region_bam="$out_dir/$plot_name/region.bam"
samtools view -b "$bam_file" "$region" > "$region_bam"

echo "Computing FIRE over input BAM file..."
source activate /scratch/midway3/kaixuan/conda/envs/fibertools
region_fire_bed="$out_dir/$plot_name/parsed/fire.bed"
ft add-nucleosomes "$region_bam" | ft fire --min-msp 10 --min-ave-msp-size 10 --skip-no-m6a - | ft fire --extract - > "$region_fire_bed"

echo "Filtering FIRE elements to the region of interest..."
source activate /scratch/midway3/kaixuan/conda/envs/fiberseq
subset_fire_bed="$out_dir/$plot_name/parsed/subset_FIRE.bed"
fire_elements_bed="$fire_dir/additional-outputs-$fire_version/fire-peaks/$fire_name-$fire_version-fire-elements.bed.gz"
bedtools intersect -a "$fire_elements_bed" -b "$region_bed" -wa > "$subset_fire_bed"

echo "Done. Results are in $out_dir/$plot_name/parsed"
