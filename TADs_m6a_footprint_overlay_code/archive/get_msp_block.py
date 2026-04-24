import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# read example accessible region from MSP file (BED12)
msp_file = "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/AL10_bc2178_19130/extracted_results/msp_by_chr/AL10_bc2178_19130.ft_extracted_msp.chr1.bed.gz"
msp_cols = [
    "chr", "start", "end", "fiber", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
]

msp_df = pd.read_csv(msp_file, sep='\t', header=None, names=msp_cols)


# expand BED12 blocks into one BED row per block for pyGenomeTracks
def _to_int_list(v):
    vals = [x for x in str(v).strip(',').split(',') if x != '']
    return [int(x) for x in vals]

records = []
for row in msp_df.itertuples(index=False):
    rel_starts = _to_int_list(row.blockStarts)
    sizes = _to_int_list(row.blockSizes)

    if len(rel_starts) != len(sizes):
        continue

    abs_starts = [row.start + s for s in rel_starts]
    abs_ends = [s + size for s, size in zip(abs_starts, sizes)]

    for i, (bs, be) in enumerate(zip(abs_starts, abs_ends), start=1):
        records.append({
            "chr": row.chr,
            "start": bs,
            "end": be,
            "name": f"{row.fiber}_block{i}",
            "score": row.score,
            "strand": row.strand
        })

msp_blocks_df = pd.DataFrame.from_records(records)
msp_blocks_df = msp_blocks_df.sort_values(["chr", "start", "end"]).reset_index(drop=True)

# write outputs for pyGenomeTracks
out_dir = "/project/spott/cshan/fiber-seq/results/TADs_Fiber_MSPs/pyGenomeTracks_inputs"
os.makedirs(out_dir, exist_ok=True)

# 1) expanded block-level BED6 (one row per block)
blocks_bed = os.path.join(out_dir, "msp_blocks_expanded.bed")
msp_blocks_df[["chr", "start", "end", "name", "score", "strand"]].to_csv(
    blocks_bed, sep='\t', header=False, index=False
)

# 2) cleaned original BED12
bed12_file = os.path.join(out_dir, "msp_fibers.bed12")
msp_df[msp_cols].to_csv(bed12_file, sep='\t', header=False, index=False)

print(f"Wrote block-level BED: {blocks_bed}")
print(f"Wrote BED12: {bed12_file}")
msp_blocks_df.head()
