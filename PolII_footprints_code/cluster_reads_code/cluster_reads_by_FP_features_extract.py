"""Per-chromosome (read, TSS) feature extraction.

Mirrors cells fp01code + fp02code of cluster_reads_by_FP_features.ipynb.
Outputs one parquet per chromosome under
  {OUT_ROOT}/per_chrom_features/{chrom}_features.parquet
Skips chromosomes whose parquet already exists.
"""

import gc
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

SAMPLE = "AL10_bc2178_19130"
HALF_WINDOW = 500   # change to 250 for a narrower window

CAGE_TSS_BED = Path(
    "/project/spott/cshan/annotations/fantom5/"
    "fantom5.hg38.LCL.consensus.CAGE_peaks.bed.gz"
)
FP_BED_ROOT = Path(
    "/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds"
)
M6A_BED_DIR = Path(
    f"/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/{SAMPLE}/"
    "extracted_results/m6a_by_chr"
)

FP_CLASSES = ["PPP", "PIC", "nucleosome_140-160", "FIRE_nucleosome", "unknown"]
NUCLEOSOME_CLASSES = {"FIRE_nucleosome", "nucleosome_140-160"}

OUT_ROOT = Path(
    f"/project/spott/cshan/fiber-seq/results/PolII/cluster_reads_by_FP_features/"
    f"{SAMPLE}_w{HALF_WINDOW}"
)
PER_CHROM_DIR = OUT_ROOT / "per_chrom_features"
OUT_ROOT.mkdir(parents=True, exist_ok=True)
PER_CHROM_DIR.mkdir(parents=True, exist_ok=True)

# Load FANTOM5 CAGE peaks (BED9). Restrict to chromosomes that have footprint BEDs.
cage = pd.read_csv(
    CAGE_TSS_BED, sep="\t", header=None,
    usecols=[0, 1, 2, 3, 4, 5],
    names=["chrom", "start", "end", "peak_id", "score", "strand"],
    dtype={"chrom": str, "start": np.int64, "end": np.int64,
           "peak_id": str, "score": str, "strand": str},
)
fp_chroms = {p.name.replace("_footprints.bed.gz", "")
             for p in (FP_BED_ROOT / FP_CLASSES[0]).glob("*_footprints.bed.gz")}
cage = cage[cage["chrom"].isin(fp_chroms)].copy()
cage["tss_coordinate"] = cage["start"].astype(np.int64)
cage = cage[["chrom", "tss_coordinate", "strand", "peak_id"]].reset_index(drop=True)
print(f"CAGE peaks on FP-supported chromosomes: {len(cage):,} "
      f"across {cage['chrom'].nunique()} chroms", flush=True)

DIST_GENERIC = ["nearest_mid_dist", "nearest_edge_dist"]
DIST_NUC = ["plus1_dist", "plus1_span", "minus1_dist", "minus1_span"]
COUNT_FEATURES = ["count", "upstream_count", "downstream_count",
                  "span_bp", "m6a_mod_in_class"]


def _expected_feature_columns():
    cols = []
    for cls in FP_CLASSES:
        for f in COUNT_FEATURES + DIST_GENERIC:
            cols.append(f"{cls}_{f}")
        if cls in NUCLEOSOME_CLASSES:
            for f in DIST_NUC:
                cols.append(f"{cls}_{f}")
    cols.append("n_m6a_mod_window")
    return cols


FEATURE_COLS = _expected_feature_columns()
INT_COLS = {c for c in FEATURE_COLS
            if c.endswith("_count") or c.endswith("_span_bp")
            or c.endswith("_m6a_mod_in_class") or c == "n_m6a_mod_window"}
print(f"{len(FEATURE_COLS)} feature columns", flush=True)


def _fp_per_tss(chrom, tss, strand, half_window):
    """dict[read_core] -> dict[class] -> list[(rel_mid, rel_start, rel_end, span)]."""
    out = defaultdict(lambda: {cls: [] for cls in FP_CLASSES})
    for cls in FP_CLASSES:
        bed_path = FP_BED_ROOT / cls / f"{chrom}_footprints.bed.gz"
        if not bed_path.exists():
            continue
        try:
            tbx = pysam.TabixFile(str(bed_path))
        except Exception:
            continue
        try:
            if chrom not in tbx.contigs:
                continue
            try:
                fetch = tbx.fetch(chrom, max(0, tss - half_window), tss + half_window)
            except ValueError:
                continue
            for rec in fetch:
                fields = rec.split("\t")
                start = int(fields[1]); end = int(fields[2])
                rc = fields[3]
                mid = float(fields[5])
                if strand == "+":
                    rel_mid = mid - tss
                    rel_start = start - tss
                    rel_end = end - tss
                else:
                    rel_mid = tss - mid
                    rel_start = tss - end
                    rel_end = tss - start
                if abs(rel_mid) > half_window:
                    continue
                out[rc][cls].append((rel_mid, rel_start, rel_end, end - start))
        finally:
            tbx.close()
    return out


def _m6a_modified_positions(chrom, tss, strand, half_window, m6a_tbx, read_cores):
    """dict[read_core] -> list[strand-oriented relative modified-m6A positions]."""
    out = defaultdict(list)
    if m6a_tbx is None or not read_cores:
        return out
    try:
        if chrom not in m6a_tbx.contigs:
            return out
        fetch = m6a_tbx.fetch(chrom, max(0, tss - half_window), tss + half_window)
    except (ValueError, AttributeError):
        return out
    lo_abs = tss - half_window
    hi_abs = tss + half_window
    for rec in fetch:
        fields = rec.split("\t")
        if len(fields) < 12:
            continue
        rc = fields[3]
        if rc not in read_cores:
            continue
        chrom_start = int(fields[1])
        block_sizes = fields[10].rstrip(",").split(",")
        block_starts = fields[11].rstrip(",").split(",")
        rels = out[rc]
        for bs, bst in zip(block_sizes, block_starts):
            try:
                sz = int(bs); st = int(bst)
            except ValueError:
                continue
            if sz == 0:
                continue
            pos = chrom_start + st
            if pos < lo_abs or pos >= hi_abs:
                continue
            rels.append((pos - tss) if strand == "+" else (tss - pos))
    return out


def _features_for_tss(chrom, tss, strand, half_window, m6a_tbx):
    fp_data = _fp_per_tss(chrom, tss, strand, half_window)
    if not fp_data:
        return {}
    read_cores = set(fp_data.keys())
    m6a_data = _m6a_modified_positions(chrom, tss, strand, half_window,
                                       m6a_tbx, read_cores)
    feats = {}
    for rc, cls_dict in fp_data.items():
        d = {}
        m6a_rel = m6a_data.get(rc, [])
        d["n_m6a_mod_window"] = len(m6a_rel)
        for cls in FP_CLASSES:
            recs = cls_dict.get(cls, [])
            if not recs:
                d[f"{cls}_count"] = 0
                d[f"{cls}_upstream_count"] = 0
                d[f"{cls}_downstream_count"] = 0
                d[f"{cls}_span_bp"] = 0
                d[f"{cls}_m6a_mod_in_class"] = 0
                continue
            ups = sum(1 for r in recs if r[0] < 0)
            dns = sum(1 for r in recs if r[0] > 0)
            span_total = sum(r[3] for r in recs)
            nearest_mid = min(abs(r[0]) for r in recs)
            edges = []
            for r in recs:
                rs, re_ = r[1], r[2]
                lo, hi = (rs, re_) if rs <= re_ else (re_, rs)
                edges.append(0.0 if lo <= 0 <= hi else min(abs(lo), abs(hi)))
            nearest_edge = min(edges)
            m6a_in = 0
            for p in m6a_rel:
                for r in recs:
                    rs, re_ = r[1], r[2]
                    lo, hi = (rs, re_) if rs <= re_ else (re_, rs)
                    if lo <= p < hi:
                        m6a_in += 1
                        break
            d[f"{cls}_count"] = len(recs)
            d[f"{cls}_upstream_count"] = ups
            d[f"{cls}_downstream_count"] = dns
            d[f"{cls}_span_bp"] = int(span_total)
            d[f"{cls}_m6a_mod_in_class"] = m6a_in
            d[f"{cls}_nearest_mid_dist"] = nearest_mid
            d[f"{cls}_nearest_edge_dist"] = nearest_edge
            if cls in NUCLEOSOME_CLASSES:
                downstream = [r for r in recs if r[0] > 0]
                upstream = [r for r in recs if r[0] < 0]
                if downstream:
                    p = min(downstream, key=lambda r: r[0])
                    d[f"{cls}_plus1_dist"] = p[0]
                    d[f"{cls}_plus1_span"] = p[3]
                if upstream:
                    m = max(upstream, key=lambda r: r[0])
                    d[f"{cls}_minus1_dist"] = m[0]
                    d[f"{cls}_minus1_span"] = m[3]
        feats[rc] = d
    return feats


def _process_chrom(chrom, cage_chrom, half_window):
    m6a_path = M6A_BED_DIR / f"{SAMPLE}.ft_extracted_m6a.{chrom}.bed.gz"
    m6a_tbx = pysam.TabixFile(str(m6a_path)) if m6a_path.exists() else None
    if m6a_tbx is None:
        print(f"    [warn] no m6A BED for {chrom}; m6A features set to 0", flush=True)
    rows = []
    try:
        for tss, strand, peak_id in cage_chrom[
            ["tss_coordinate", "strand", "peak_id"]
        ].itertuples(index=False, name=None):
            feats = _features_for_tss(chrom, int(tss), strand, half_window, m6a_tbx)
            for rc, d in feats.items():
                d["read_core"] = rc
                d["peak_id"] = peak_id
                d["chromosome"] = chrom
                d["tss_coordinate"] = int(tss)
                d["tss_strand"] = strand
                rows.append(d)
    finally:
        if m6a_tbx is not None:
            m6a_tbx.close()
    if not rows:
        return pd.DataFrame(columns=["read_core", "peak_id", "chromosome",
                                     "tss_coordinate", "tss_strand"] + FEATURE_COLS)
    df = pd.DataFrame(rows)
    for c in FEATURE_COLS:
        if c not in df.columns:
            df[c] = np.nan
        if c in INT_COLS:
            df[c] = df[c].fillna(0).astype(np.int32)
    return df[["read_core", "peak_id", "chromosome",
               "tss_coordinate", "tss_strand"] + FEATURE_COLS]


chrom_list = sorted(cage["chrom"].unique())
print(f"Processing {len(chrom_list)} chromosomes", flush=True)
for chrom in chrom_list:
    out_path = PER_CHROM_DIR / f"{chrom}_features.parquet"
    if out_path.exists():
        print(f"  {chrom}: cached -> {out_path.name}", flush=True)
        continue
    cage_chrom = cage[cage["chrom"] == chrom]
    print(f"  {chrom}: {len(cage_chrom):,} CAGE peaks", flush=True)
    feat_df = _process_chrom(chrom, cage_chrom, HALF_WINDOW)
    feat_df.to_parquet(out_path, index=False)
    print(f"    -> {len(feat_df):,} (read, TSS) feature rows", flush=True)
    del feat_df, cage_chrom
    gc.collect()
print("Done with feature extraction.", flush=True)
