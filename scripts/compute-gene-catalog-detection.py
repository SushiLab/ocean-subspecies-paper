#!/usr/bin/env python

import gzip
import sys
from datetime import datetime
from pathlib import Path

import polars as pl

def log_msg(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)

def fmt_n(n: int) -> str:
    return f"{n:,}"

#####
# Parameters -- edit these paths as needed
#####
profile_files = [
    "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/quantification/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-genes-cds-w_eukarya/test-merge/merged_profiles_prok-metag-insert.tsv.gz",
    "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/quantification/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-genes-cds-w_eukarya/test-merge/merged_profiles_lsf-metag-insert.tsv.gz",
]
sample_map_file = "samples.map"   # two columns, no header: old_name<TAB>new_name
out_dir = Path("GENE-CATALOG-FULL-DETECTION.parquet")   # a directory of part files (one
                                                          # per profile file), read back as
                                                          # one logical dataset via
                                                          # pl.scan_parquet(out_dir / "*.parquet")
#out_tsv = "GENE-CATALOG-FULL-DETECTION.tsv.gz"           # single-file tsv.gz conversion of
                                                           # the above, written at the end --
                                                           # matches what
                                                           # compute-horizontal-coverage.py
                                                           # currently reads

#####
# 1. Load the old -> new sample rename map, keeping only the METAG entries -- this is the
#    definitive list of the 1,671 valid samples for this analysis. Any raw column NOT
#    listed here gets excluded entirely below, regardless of what it looks like (even an
#    already-real, correctly METAG-suffixed name): sample_map is a whitelist of valid
#    samples, not just a source of rename rules.
#####
sample_map_df = pl.read_csv(sample_map_file, separator="\t", has_header=False, new_columns=["old", "new"])
n_total = sample_map_df.height
sample_map_df = sample_map_df.filter(pl.col("new").str.ends_with("_METAG"))
n_metat_dropped = n_total - sample_map_df.height
sample_map = dict(zip(sample_map_df["old"].to_list(), sample_map_df["new"].to_list()))
log_msg(f"sample_map: {fmt_n(sample_map_df.height)} METAG old->new mappings loaded from "
        f"{sample_map_file} ({fmt_n(n_metat_dropped)} METAT entries filtered out)")

#####
# 2. Process each profile file with a fully lazy scan -> rename -> unpivot -> filter ->
#    sink_parquet pipeline -- polars' streaming engine processes this without ever
#    materializing the full wide (or unpivoted) table in memory, so there's no manual
#    chunking to write or tune, unlike the R version this replaces.
#####

def find_header_skip(path: str, marker: str = "#reference", max_lines: int = 20) -> int:
    """Scan forward for the line whose first tab-separated field is `marker`, returning
    how many lines precede it. Mirrors the R version's skip="#reference" behavior: some
    profile files have a leading comment line, some don't."""
    with gzip.open(path, "rt") as f:
        for i in range(max_lines):
            line = f.readline()
            if not line:
                break
            if line.split("\t", 1)[0] == marker:
                return i
    raise RuntimeError(f"Could not find a '{marker}' header line in the first {max_lines} lines of {path}")

def peek_header(path: str, skip_rows: int) -> list[str]:
    with gzip.open(path, "rt") as f:
        for _ in range(skip_rows):
            f.readline()
        header_line = f.readline().rstrip("\n")
    return header_line.split("\t")

out_dir.mkdir(exist_ok=True)
samples_seen: set[str] = set()
raw_names_seen: set[str] = set()
gene_ids_seen: set[str] = set()
samples_detected: set[str] = set()
n_pairs_total = 0

for i, path in enumerate(profile_files):
    log_msg(f"Reading profile table: {path}")
    n_skip = find_header_skip(path)
    all_names = peek_header(path, n_skip)
    ref_col = all_names[0]
    has_length_col = len(all_names) > 1 and all_names[1] == "length"
    raw_sample_names_all = all_names[2:] if has_length_col else all_names[1:]
    raw_names_seen |= set(raw_sample_names_all)

    raw_sample_names = [n for n in raw_sample_names_all if n in sample_map]
    n_excluded = len(raw_sample_names_all) - len(raw_sample_names)
    if n_excluded > 0:
        log_msg(f"  excluding {n_excluded} sample(s) not present in sample_map's "
                f"{fmt_n(len(sample_map))} valid entries")

    new_names_preview = [sample_map.get(n, n) for n in raw_sample_names]
    dup_mask = [n in samples_seen for n in new_names_preview]
    raw_names_to_keep = [n for n, d in zip(raw_sample_names, dup_mask) if not d]
    if any(dup_mask):
        dropped_new_names = [n for n, d in zip(new_names_preview, dup_mask) if d]
        log_msg(f"  dropping {len(dropped_new_names)} duplicate sample(s) already present in "
                f"an earlier file: {', '.join(dropped_new_names[:10])}")

    new_names = [sample_map[n] for n in raw_names_to_keep]   # every n here is guaranteed to
                                                                # be a sample_map key now
    rename_map = {ref_col: "Gene_catalog_rep", **dict(zip(raw_names_to_keep, new_names))}
    samples_seen |= set(new_names)

    lf = (
        pl.scan_csv(path, separator="\t", skip_rows=n_skip)
          .select([ref_col] + raw_names_to_keep)   # drop "length" and duplicate columns
                                                     # before they're ever parsed -- polars
                                                     # pushes this projection down into the
                                                     # CSV scan itself
          .rename(rename_map)
          .unpivot(index="Gene_catalog_rep", variable_name="Sample", value_name="Abundance")
          .filter(pl.col("Abundance") > 0)
          .select("Gene_catalog_rep", "Sample")
          .unique()
    )

    part_path = out_dir / f"part-{i}.parquet"
    lf.sink_parquet(part_path)

    part = pl.read_parquet(part_path)   # small (sparse, detected-only) -- fine to read back
                                          # just for the summary counts below
    n_pairs_total += part.height
    gene_ids_seen |= set(part["Gene_catalog_rep"].unique().to_list())
    samples_detected |= set(part["Sample"].unique().to_list())

    log_msg(f"  {path}: {len(new_names)} samples kept (all via sample_map), "
            f"{fmt_n(part.height)} detected (Gene_catalog_rep, Sample) pairs")

#####
# 3. Confirm every sample_map entry was actually found in at least one profile file
#####
map_not_found = set(sample_map_df["old"].to_list()) - raw_names_seen
if map_not_found:
    log_msg(f"{fmt_n(len(map_not_found))}/{fmt_n(sample_map_df.height)} sample_map entries "
            f"were never found as a column in any profile file (e.g. "
            f"{', '.join(list(map_not_found)[:10])})")
else:
    log_msg(f"All {fmt_n(sample_map_df.height)} sample_map entries were found in at least one profile file")

#####
# 4. Summary (already written incrementally above, one Parquet part file per profile file)
#####
log_msg(f"Combined: {fmt_n(n_pairs_total)} detected (Gene_catalog_rep, Sample) pairs, "
        f"{fmt_n(len(gene_ids_seen))} distinct gene catalog representatives, "
        f"{len(samples_seen)} samples processed ({len(samples_detected)} with at least one detection)")
log_msg(f"Wrote {out_dir}/part-*.parquet ({len(profile_files)} files)")

#####
# 5. Convert the Parquet dataset to a single tsv.gz file, for compatibility with scripts
#    (like compute-horizontal-coverage.py) that read this format. sink_csv streams the
#    conversion (its own docs describe it as suited to "results that are larger than RAM"),
#    so this doesn't re-materialize the full dataset the way a plain read+write would.
#####
if out_tsv:
    log_msg(f"Converting {out_dir}/*.parquet -> {out_tsv}")
    pl.scan_parquet(f"{out_dir}/*.parquet").sink_csv(out_tsv, separator="\t", compression="gzip")
    log_msg(f"Wrote {out_tsv}")
