import argparse
import csv
import gzip
import os
import re
from glob import glob
from collections import defaultdict


METADATA_DEFAULT = "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
OUTPUT_ROOT_DEFAULT = "/project/spott/cshan/fiber-seq/results/TADs_m6a_footprint_overlay"
TAD_BED = "/project/spott/cshan/fiber-seq/results/TADs_Fiber_MSP/pygenometracks_overlay/inputs/tads.bed"
TAD_BOUNDARY_BED = "/project/spott/cshan/annotations/hic/ENCFF156ECM.bedpe.gz"
CAGE_BED = "/project/spott/cshan/fiber-seq/results/TADs_m6a_footprint_overlay/inputs/cage_lcl_extended.bed.gz"
HIC_MATRIX = "/project/spott/cshan/annotations/hic/ENCFF216QQM.cool"
CTCF_BW = "/project/spott/cshan/annotations/CTCF_ENCFF636ODI.bigWig"
DNASE_BW = "/project/spott/cshan/annotations/DNase.ENCFF743ULW.bigWig"
SMALL_REGION_TAD_THRESHOLD = 500_000
POLII_FP_BW_ROOT = "/project/spott/cshan/fiber-seq/results/PolII/FIRE_combined_footprints/joint_trained_tracks"

FP_BINS = [
    ("fp_10_30bp", 10, 30, "#08306b", "Footprint 10-30bp"),
    ("fp_40_60bp", 40, 60, "#c6dbef", "Footprint 40-60bp"),
    ("fp_60_80bp", 60, 80, "#e377c2", "Footprint 60-80bp"),
    ("fp_140_160bp", 140, 160, "#ff7f0e", "Footprint 140-160bp"),
]

POLII_COLORS = [
    "#08306b",
    "#2171b5",
    "#6baed6",
    "#9ecae1",
    "#fdd0a2",
    "#f16913",
    "#a50f15",
]


def discover_polii_bigwigs(chrom, root=POLII_FP_BW_ROOT):
    """Find PolII footprint coverage bigWigs for a chromosome and parse bp ranges from filenames."""
    range_re = re.compile(r"combined_(chr[^_]+)_(\d+)-(\d+)bp_fp_cov\.(?:bw|bigwig)$", re.IGNORECASE)
    candidates = []
    search_globs = [
        os.path.join(root, chrom, "*.bw"),
        os.path.join(root, chrom, "*.bigWig"),
        os.path.join(root, f"*{chrom}*.bw"),
        os.path.join(root, f"*{chrom}*.bigWig"),
    ]
    for pattern in search_globs:
        candidates.extend(glob(pattern))

    tracks = []
    for file_path in sorted(set(candidates)):
        m = range_re.search(os.path.basename(file_path))
        if not m:
            continue
        bw_chrom, start_bp, end_bp = m.group(1), int(m.group(2)), int(m.group(3))
        if bw_chrom != chrom:
            continue
        tracks.append(
            {
                "path": file_path,
                "start_bp": start_bp,
                "end_bp": end_bp,
                "title": f"PolII FP coverage {start_bp}-{end_bp} bp",
            }
        )

    tracks.sort(key=lambda t: (t["start_bp"], t["end_bp"], t["path"]))
    for i, track in enumerate(tracks):
        track["color"] = POLII_COLORS[i % len(POLII_COLORS)]
    return tracks


def normalize_chrom(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def parse_region(region):
    chrom, coords = region.split(":", 1)
    start, end = coords.split("-", 1)
    return normalize_chrom(chrom), int(start), int(end)


def region_slug(chrom, start, end):
    return f"{chrom}_{start}_{end}"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build per-region per-read m6A and FiberHMM tracks plus a pyGenomeTracks config."
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--region", default=None, help="chr:start-end")
    parser.add_argument("--chrom", default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--metadata", default=METADATA_DEFAULT)
    parser.add_argument("--output-root", default=OUTPUT_ROOT_DEFAULT)
    parser.add_argument("--m6a-summary-track", default=None, help="Optional region-level m6A bedGraph/bigWig from the R step.")
    parser.add_argument("--5mc-summary-track", dest="mc5_summary_track", default=None, help="Optional region-level 5mC bedGraph/bigWig from the R step.")
    parser.add_argument("--max-reads", type=int, default=250)
    parser.add_argument("--hic-depth", type=int, default=5000, help="Hi-C triangle depth in bp (default: 5000).")
    return parser.parse_args()


def load_sample_row(metadata_path, sample):
    with open(metadata_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["sample_name"] == sample:
                return row
    raise SystemExit(f"Sample {sample} was not found in {metadata_path}")


def parse_int_list(value):
    return [int(x) for x in value.strip(",").split(",") if x]


def overlap(start, end, region_start, region_end):
    return max(0, min(end, region_end) - max(start, region_start))


def clip_interval(start, end, region_start, region_end):
    clipped_start = max(start, region_start)
    clipped_end = min(end, region_end)
    if clipped_end <= clipped_start:
        return None
    return clipped_start, clipped_end


def collect_m6a_records(path, chrom, region_start, region_end):
    records = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[0] != chrom:
                continue
            row_start = int(fields[1])
            row_end = int(fields[2])
            if overlap(row_start, row_end, region_start, region_end) <= 0:
                continue
            block_sizes = parse_int_list(fields[10])
            block_starts = parse_int_list(fields[11])
            all_blocks = list(zip(block_starts, block_sizes))
            # First and last blocks are read-span anchors, not real m6A calls
            real_blocks = all_blocks[1:-1] if len(all_blocks) >= 2 else []
            records[fields[3]] = {
                "chrom": chrom,
                "name": fields[3],
                "strand": fields[5],
                "start": row_start,
                "end": row_end,
                "item_rgb": fields[8],
                "blocks": real_blocks,
            }
    return records


def collect_fp_records(path, chrom, region_start, region_end):
    records = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[0] != chrom:
                continue
            row_start = int(fields[1])
            row_end = int(fields[2])
            if overlap(row_start, row_end, region_start, region_end) <= 0:
                continue
            read_id = fields[3]
            block_starts = parse_int_list(fields[8])
            block_sizes = parse_int_list(fields[9])
            records[read_id].append(
                {
                    "chrom": chrom,
                    "name": read_id,
                    "start": row_start,
                    "end": row_end,
                    "item_rgb": fields[6],
                    "blocks": list(zip(block_starts, block_sizes)),
                }
            )
    return records


def make_read_spans(m6a_records, fp_records, region_start, region_end):
    read_spans = {}
    for read_id, rec in m6a_records.items():
        clipped = clip_interval(rec["start"], rec["end"], region_start, region_end)
        if clipped is None:
            continue
        read_spans[read_id] = {
            "start": clipped[0],
            "end": clipped[1],
            "overlap_bp": clipped[1] - clipped[0],
            "feature_count": len(rec["blocks"]),
        }

    for read_id, recs in fp_records.items():
        starts = [rec["start"] for rec in recs]
        ends = [rec["end"] for rec in recs]
        clipped = clip_interval(min(starts), max(ends), region_start, region_end)
        if clipped is None:
            continue
        current = read_spans.setdefault(
            read_id,
            {
                "start": clipped[0],
                "end": clipped[1],
                "overlap_bp": clipped[1] - clipped[0],
                "feature_count": 0,
            },
        )
        current["start"] = min(current["start"], clipped[0])
        current["end"] = max(current["end"], clipped[1])
        current["overlap_bp"] = current["end"] - current["start"]
        current["feature_count"] += sum(len(rec["blocks"]) for rec in recs)

    return read_spans


def select_read_order(read_spans, max_reads):
    ordered = sorted(
        read_spans.items(),
        key=lambda item: (-item[1]["overlap_bp"], -item[1]["feature_count"], item[1]["start"], item[0]),
    )
    if max_reads and len(ordered) > max_reads:
        ordered = ordered[:max_reads]
    return [read_id for read_id, _ in ordered]


def aggregate_clipped_blocks(block_records, chrom_start, region_start, region_end, size_range=None, chrom_end=None):
    """
    Return blocks as (relative_start, size) pairs relative to chrom_start.

    chrom_start is the chromStart (col 2) of the BED12 row being written —
    i.e. span_start after clipping.  BED12 requires that the first block
    starts at 0 and the last block reaches chromEnd, so we insert 1-bp
    anchor blocks at those positions when no real block covers them.
    """
    clipped_blocks = []
    for record in block_records:
        for rel_start, size in record["blocks"]:
            if size_range is not None and not (size_range[0] <= size <= size_range[1]):
                continue
            block_start = record["start"] + rel_start
            block_end = block_start + size
            clipped = clip_interval(block_start, block_end, region_start, region_end)
            if clipped is None:
                continue
            relative = clipped[0] - chrom_start
            clipped_size = clipped[1] - clipped[0]
            clipped_blocks.append((relative, clipped_size))

    clipped_blocks.sort()

    # BED12 spec: first block relative start must be 0.
    if not clipped_blocks or clipped_blocks[0][0] != 0:
        clipped_blocks.insert(0, (0, 1))

    # BED12 spec: last block must reach chromEnd.
    if chrom_end is not None:
        required_end = chrom_end - chrom_start
        last_start, last_size = clipped_blocks[-1]
        if last_start + last_size < required_end:
            clipped_blocks.append((required_end - 1, 1))

    return clipped_blocks


def write_bed12(path, rows):
    with open(path, "w") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def build_tracks(sample, chrom, region_start, region_end, read_order, read_spans, m6a_records, fp_records, region_dir):
    baseline_rows = []
    m6a_rows = []
    fp_rows = {key: [] for key, *_ in FP_BINS}

    for read_id in read_order:
        span = read_spans[read_id]
        span_start = span["start"]
        span_end = span["end"]
        span_size = span_end - span_start
        if span_size <= 0:
            continue

        baseline_rows.append(
            [
                chrom,
                span_start,
                span_end,
                read_id,
                0,
                ".",
                span_start,
                span_end,
                "220,220,220",
                1,
                f"{span_size},",
                "0,",
            ]
        )

        if read_id in m6a_records:
            blocks = aggregate_clipped_blocks(
                [m6a_records[read_id]], span_start, region_start, region_end, chrom_end=span_end
            )
            if blocks:
                m6a_rows.append(
                    [
                        chrom,
                        span_start,
                        span_end,
                        read_id,
                        0,
                        m6a_records[read_id]["strand"],
                        span_start,
                        span_end,
                        "123,50,148",
                        len(blocks),
                        ",".join(str(size) for _, size in blocks) + ",",
                        ",".join(str(start) for start, _ in blocks) + ",",
                    ]
                )

        if read_id in fp_records:
            for key, lower, upper, _, _ in FP_BINS:
                blocks = aggregate_clipped_blocks(
                    fp_records[read_id], span_start, region_start, region_end, size_range=(lower, upper), chrom_end=span_end
                )
                # Only write a row if there are real footprint blocks beyond the anchor
                if len(blocks) > 1 or (len(blocks) == 1 and blocks[0][0] == 0 and blocks[0][1] > 1):
                    fp_rows[key].append(
                        [
                            chrom,
                            span_start,
                            span_end,
                            read_id,
                            0,
                            ".",
                            span_start,
                            span_end,
                            "0,0,0",
                            len(blocks),
                            ",".join(str(size) for _, size in blocks) + ",",
                            ",".join(str(start) for start, _ in blocks) + ",",
                        ]
                    )

    paths = {
        "baseline": os.path.join(region_dir, f"{sample}.{region_slug(chrom, region_start, region_end)}.reads_baseline.bed12"),
        "m6a": os.path.join(region_dir, f"{sample}.{region_slug(chrom, region_start, region_end)}.m6a_reads.bed12"),
    }
    write_bed12(paths["baseline"], baseline_rows)
    write_bed12(paths["m6a"], m6a_rows)

    for key in fp_rows:
        paths[key] = os.path.join(region_dir, f"{sample}.{region_slug(chrom, region_start, region_end)}.{key}.bed12")
        write_bed12(paths[key], fp_rows[key])

    return paths


def auto_cpg_bed(sample):
    path = f"/project/spott/cshan/fiber-seq/results/fire_CpG/{sample}/{sample}_CPG.combined.bed.gz"
    return path if os.path.exists(path) else None


def summary_track_type(path):
    lower = path.lower()
    if lower.endswith(".bw") or lower.endswith(".bigwig"):
        return "bigwig"
    if lower.endswith(".bedgraph") or lower.endswith(".bdg"):
        return "bedgraph"
    raise SystemExit(f"Unsupported m6A summary track type for {path}")


def tad_track_type(path):
    lower = path.lower()
    if lower.endswith(".bedpe") or lower.endswith(".bedpe.gz"):
        return "links"
    return "bed"


def write_tracks_ini(
    path,
    sample,
    chrom,
    track_paths,
    m6a_summary_track=None,
    mc5_summary_track=None,
    cpg_bed=None,
    hic_depth=200000,
    tad_bed=TAD_BED,
    tad_title="TADs (Arrowhead, GM12878)",
    show_hic=True,
):
    polii_tracks = discover_polii_bigwigs(chrom)
    with open(path, "w") as handle:
        handle.write("[x-axis]\nfontsize = 10\nwhere = top\n\n")

        if show_hic:
            handle.write("[hic]\n")
            handle.write(f"file = {HIC_MATRIX}\n")
            handle.write("title = ENCFF216QQM Hi-C\n")
            handle.write(f"depth = {hic_depth}\nheight = 6\ncolormap = Reds\nfile_type = hic_matrix\n\n")
            handle.write("[spacer]\nheight = 0.3\n\n")

        handle.write("[tads]\n")
        handle.write(f"file = {tad_bed}\n")
        handle.write(f"title = {tad_title}\n")
        handle.write(
            "height = 2\ncolor = none\nborder_color = #e41a1c\nline_width = 1.5\nlabels = false\ndisplay = interleaved"
            f"\nfile_type = {tad_track_type(tad_bed)}\n\n"
        )

        for title, file_path, color in [
            ("CTCF", CTCF_BW, "black"),
            ("DNase", DNASE_BW, "#1f77b4"),
        ]:
            handle.write("[spacer]\nheight = 0.3\n\n")
            handle.write(f"[{title.lower()}]\n")
            handle.write(f"file = {file_path}\n")
            handle.write(f"title = {title}\n")
            handle.write(f"height = 2\ncolor = {color}\nmin_value = 0\nshow_data_range = true\nfile_type = bigwig\n\n")

        handle.write("[spacer]\nheight = 0.3\n\n")
        handle.write("[cage]\n")
        handle.write(f"file = {CAGE_BED}\n")
        handle.write("title = FANTOM5 CAGE\n")
        handle.write("height = 1.5\ncolor = #2ca02c\nborder_color = #1b6e1b\nlabels = false\ndisplay = collapsed\nfile_type = bed\n\n")

        if cpg_bed is not None:
            handle.write("[spacer]\nheight = 0.3\n\n")
            handle.write("[cpg]\n")
            handle.write(f"file = {cpg_bed}\n")
            handle.write(f"title = CpG methylation pb-CpG ({sample})\n")
            handle.write("height = 1.5\ncolor = #8856a7\nmin_value = 0\nmax_value = 100\nshow_data_range = true\nfile_type = bedgraph\n\n")

        for i, track in enumerate(polii_tracks, start=1):
            handle.write("[spacer]\nheight = 0.3\n\n")
            handle.write(f"[polii_fp_cov_{i}]\n")
            handle.write(f"file = {track['path']}\n")
            handle.write(f"title = {track['title']}\n")
            handle.write(
                f"height = 1.2\ncolor = {track['color']}\nmin_value = 0\nshow_data_range = true\nfile_type = bigwig\n\n"
            )

        if (
            m6a_summary_track is not None
            and os.path.exists(m6a_summary_track)
            and os.path.getsize(m6a_summary_track) > 0
        ):
            handle.write("[spacer]\nheight = 0.3\n\n")
            handle.write("[m6a_summary]\n")
            handle.write(f"file = {m6a_summary_track}\n")
            handle.write("title = m6A fraction (footprintR)\n")
            handle.write("height = 1.6\ncolor = #7B3294\nmin_value = 0\nmax_value = 1\nshow_data_range = true\n")
            handle.write(f"file_type = {summary_track_type(m6a_summary_track)}\n\n")

        if (
            mc5_summary_track is not None
            and os.path.exists(mc5_summary_track)
            and os.path.getsize(mc5_summary_track) > 0
        ):
            handle.write("[spacer]\nheight = 0.3\n\n")
            handle.write("[5mc_summary]\n")
            handle.write(f"file = {mc5_summary_track}\n")
            handle.write("title = 5mC fraction (footprintR)\n")
            handle.write("height = 1.6\ncolor = #2c7bb6\nmin_value = 0\nmax_value = 1\nshow_data_range = true\n")
            handle.write(f"file_type = {summary_track_type(mc5_summary_track)}\n\n")

        handle.write("[spacer]\nheight = 0.3\n\n")
        handle.write("[reads_baseline]\n")
        handle.write(f"file = {track_paths['baseline']}\n")
        handle.write("title = Per-read overlay\n")
        handle.write("height = 14\ncolor = #d9d9d9\nborder_color = #d9d9d9\nlabels = false\ndisplay = stacked\nfile_type = bed\n\n")

        handle.write("[m6a_reads]\n")
        handle.write(f"file = {track_paths['m6a']}\n")
        handle.write("title = \n")
        handle.write("height = 14\ncolor = #7B3294\nborder_color = #7B3294\nlabels = false\ndisplay = stacked\noverlay_previous = share-y\nfile_type = bed\n\n")

        for key, _, _, color, title in FP_BINS:
            handle.write(f"[{key}]\n")
            handle.write(f"file = {track_paths[key]}\n")
            handle.write("title = \n")
            handle.write(f"height = 14\ncolor = {color}\nborder_color = {color}\nlabels = false\ndisplay = stacked\noverlay_previous = share-y\nfile_type = bed\n\n")


def main():
    args = parse_args()
    if args.region:
        chrom, region_start, region_end = parse_region(args.region)
    else:
        if args.chrom is None or args.start is None or args.end is None:
            raise SystemExit("Provide --region or --chrom/--start/--end")
        chrom = normalize_chrom(args.chrom)
        region_start = args.start
        region_end = args.end

    if region_start < 1 or region_end < region_start:
        raise SystemExit("Invalid genomic coordinates")

    region_span = region_end - region_start

    sample_row = load_sample_row(args.metadata, args.sample)
    fire_dir = sample_row["fire_dir"]
    region_name = region_slug(chrom, region_start, region_end)
    region_dir = os.path.join(args.output_root, args.sample, chrom, region_name)
    os.makedirs(region_dir, exist_ok=True)

    m6a_path = os.path.join(
        fire_dir, "extracted_results", "m6a_by_chr", f"{args.sample}.ft_extracted_m6a.{chrom}.bed.gz"
    )
    fp_path = os.path.join(
        "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged",
        args.sample,
        "joint_trained_results",
        f"{args.sample}_{chrom}_fp.bed.gz",
    )

    if not os.path.exists(m6a_path):
        raise SystemExit(f"Missing sample m6A BED12: {m6a_path}")
    if not os.path.exists(fp_path):
        raise SystemExit(f"Missing sample FiberHMM footprint BED: {fp_path}")

    m6a_records = collect_m6a_records(m6a_path, chrom, region_start, region_end)
    fp_records = collect_fp_records(fp_path, chrom, region_start, region_end)
    read_spans = make_read_spans(m6a_records, fp_records, region_start, region_end)
    read_order = select_read_order(read_spans, args.max_reads)
    if not read_order:
        raise SystemExit("No reads overlapped the requested region.")

    track_paths = build_tracks(
        args.sample, chrom, region_start, region_end, read_order, read_spans, m6a_records, fp_records, region_dir
    )

    read_table_path = os.path.join(region_dir, f"{args.sample}.{region_name}.selected_reads.tsv")
    with open(read_table_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["read_id", "clipped_start", "clipped_end", "overlap_bp", "feature_count"])
        for read_id in read_order:
            row = read_spans[read_id]
            writer.writerow([read_id, row["start"], row["end"], row["overlap_bp"], row["feature_count"]])

    ini_path = os.path.join(region_dir, f"{args.sample}.{region_name}.tracks.ini")
    use_tad_boundaries = region_span < SMALL_REGION_TAD_THRESHOLD
    write_tracks_ini(
        ini_path,
        sample=args.sample,
        chrom=chrom,
        track_paths=track_paths,
        m6a_summary_track=args.m6a_summary_track,
        mc5_summary_track=args.mc5_summary_track,
        cpg_bed=auto_cpg_bed(args.sample),
        hic_depth=args.hic_depth,
        tad_bed=TAD_BOUNDARY_BED if use_tad_boundaries else TAD_BED,
        tad_title="TAD boundaries (ENCFF156ECM)" if use_tad_boundaries else "TADs (Arrowhead, GM12878)",
        show_hic=not use_tad_boundaries,
    )

    print(ini_path)


if __name__ == "__main__":
    main()
