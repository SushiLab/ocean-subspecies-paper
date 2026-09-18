#!/usr/bin/env python

import gzip
import sys
from datetime import datetime

import polars as pl

def log_msg(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)

def fmt_n(n: int) -> str:
    return f"{n:,}"

#####
# Parameters -- edit these paths as needed
#####
species_motu_file = "species_over10g_x_motus.tsv"   # local export of the "Species (genomes
                                              # >= 10) x mOTUs" Google Sheet
motus_base_file   = "/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/share/paolil/OCEAN_POPGEN/gom-basecoverage.motus.gz"
genome_summary_url = "https://www.microbiomics.io/ocean/suppl_data/genomes-summary.csv.gz"
sample_map_file    = "samples.map"   # same METAG whitelist as preprocess-profiles.py
out = "vertical-coverage-per-sample-per-species.tsv"

#####
# 1. Load the METAG old->new sample map -- same file as preprocess-profiles.py.
#####
sample_map_df = pl.read_csv(sample_map_file, separator="\t", has_header=False, new_columns=["old", "new"])
sample_map_df = sample_map_df.filter(pl.col("new").str.ends_with("_METAG"))
sample_map = dict(zip(sample_map_df["old"].to_list(), sample_map_df["new"].to_list()))
log_msg(f"sample_map: {fmt_n(len(sample_map))} METAG old->new mappings loaded from {sample_map_file}")

# These two are listed in sample_map (as OLD-column values) but were never actually found
# as columns in either profile file in the horizontal-coverage pipeline (confirmed there:
# "2/1,671 sample_map entries were never found..."). Excluding them here too keeps both
# pipelines computing over the exact same population. NOTE: this must check the raw/old
# name, not the resolved new one -- an earlier version of this script subtracted these OLD
# values from a set of NEW values, which is a silent no-op (they were never in that set to
# begin with).
samples_confirmed_absent = {"TARA_N000000807-S100_METAG", "TARA_N000002223_METAG"}
# valid_samples is for genome_summary's derived Sample values below, which are already in
# "new"/real form (unlike motus_base's raw barcode-style columns, handled via
# resolve_sample() instead) -- so this needs the NEW-style equivalents of the confirmed-
# absent OLD entries, looked up from sample_map itself rather than hardcoded, since only
# one of the two ("TARA_N000000807-S100_METAG" -> "TARA_SAMEA2620959_METAG") is known here.
valid_samples = set(sample_map.values()) - {
    sample_map[old] for old in samples_confirmed_absent if old in sample_map
}

# motus_base carries at least one raw column name that doesn't match sample_map's "old"
# convention by a simple suffix difference -- it needs an actual character fix, the same
# one the original R script applied by hand:
#   mutate(samples = ifelse(samples == "TARA_N000000807SUB_METAG",
#                            "TARA_N000000807-S100_METAG", samples))
# Confirmed against the real sample_map: its "old" column has the hyphenated form
# ("TARA_N000000807-S100_METAG"), while motus_base's header has the unhyphenated one.
KNOWN_NAME_FIXUPS = {
    "TARA_N000000807SUB_METAG": "TARA_N000000807-S100_METAG",
}

def resolve_sample(raw_name: str) -> str | None:
    """Translate a raw motus_base column name to its real sample name via sample_map,
    trying the name as given first, then with a trailing _METAG/_METAT stripped (the
    original R script derived its own join key this way: barcode = gsub("_METAG", "",
    sample)). Returns None if neither resolves, or if the (post-fixup) raw name is one of
    the confirmed-absent entries -- excluded, not kept as-is, matching the strict-whitelist
    approach used in preprocess-profiles.py: sample_map is the definitive list of valid
    samples, not just a source of optional rename rules."""
    raw_name = KNOWN_NAME_FIXUPS.get(raw_name, raw_name)
    if raw_name in samples_confirmed_absent:
        return None
    if raw_name in sample_map:
        return sample_map[raw_name]
    for suffix in ("_METAG", "_METAT"):
        if raw_name.endswith(suffix):
            stripped = raw_name[: -len(suffix)]
            if stripped in sample_map:
                return sample_map[stripped]
    return None

#####
# 2. Load the species -> motu whitelist, expanding semicolon-separated multi-motu matches
#    (e.g. "majority vote after agreg n=2") into one row per (species, single motu)
#####
log_msg(f"Reading species-motu whitelist: {species_motu_file}")
species_motu_raw = pl.read_csv(species_motu_file, separator="\t", columns=["species", "motu", "type"])
# NOTE: only these 3 columns are read -- the others (n, intersect, motu_n, species_n,
# species_r, motu_r) become semicolon-separated per-motu lists on multi-motu rows (e.g.
# motu_n = "183;58" for a 2-motu match), not plain integers, so letting polars infer their
# dtype from the whole file fails. Since none of them are needed here, they're never parsed.
species_motu = (
    species_motu_raw
    .with_columns(pl.col("motu").str.split(";"))
    .explode("motu", empty_as_null=True)
    .rename({"species": "Species", "motu": "motu_match", "type": "motu_match_type"})
)
log_msg(f"species_motu: {fmt_n(species_motu.height)} (Species, motu_match) rows across "
        f"{species_motu['Species'].n_unique()} species after expanding multi-motu matches")

#####
# 3. Load motus_base: wide (#consensus_taxonomy + one column per sample), restrict columns
#    to valid METAG samples and rows to the motus actually referenced above, then unpivot
#    to long. This is a DENSE matrix (every motu has a value -- possibly 0 -- for every
#    sample), unlike the sparse gene-detection data in the other pipeline, so this join
#    already produces every sample's vertical_coverage without a separate grid expansion.
#####
def find_header_skip(path: str, marker: str, max_lines: int = 5) -> int:
    with gzip.open(path, "rt") as f:
        for i in range(max_lines):
            line = f.readline()
            if not line:
                break
            if line.split("\t", 1)[0] == marker:
                return i
    raise RuntimeError(f"Could not find a '{marker}' header line in the first {max_lines} lines of {path}")

n_skip = find_header_skip(motus_base_file, marker="#consensus_taxonomy")
header_cols = pl.read_csv(motus_base_file, separator="\t", skip_rows=n_skip, n_rows=0).columns
all_sample_cols = [c for c in header_cols if c != "#consensus_taxonomy"]

resolved = [(c, resolve_sample(c)) for c in all_sample_cols]
n_direct = sum(1 for c, r in resolved if r is not None and r == sample_map.get(c))
n_via_strip = sum(1 for c, r in resolved if r is not None and c not in sample_map)
raw_to_keep = [c for c, r in resolved if r is not None and r not in samples_confirmed_absent]
new_names   = [r for c, r in resolved if r is not None and r not in samples_confirmed_absent]
log_msg(f"motus_base header: {len(all_sample_cols)} sample columns, {len(raw_to_keep)} resolved "
        f"to valid samples ({fmt_n(n_direct)} direct match, {fmt_n(n_via_strip)} matched after "
        f"stripping a _METAG/_METAT suffix, {len(all_sample_cols) - len(raw_to_keep)} unresolved/excluded)")
if raw_to_keep:
    example = next((c for c, r in resolved if r is not None), None)
    log_msg(f"  example resolved column: {example!r} -> {resolve_sample(example)!r}")

rename_map = dict(zip(raw_to_keep, new_names))
referenced_motus = set(species_motu["motu_match"].to_list())

motus_base = (
    pl.scan_csv(motus_base_file, separator="\t", skip_rows=n_skip, infer_schema_length=100000)
      .select(["#consensus_taxonomy"] + raw_to_keep)
      .rename(rename_map)
      .with_columns(pl.col("#consensus_taxonomy").str.extract(r"\[([^\]]+)\]$", 1).alias("motu"))
      .filter(pl.col("motu").is_in(list(referenced_motus)))
      .drop("#consensus_taxonomy")
      .unpivot(index="motu", variable_name="Sample", value_name="vertical_coverage")
      .collect()
)
log_msg(f"motus_base: {fmt_n(motus_base.height)} (motu, Sample) rows for "
        f"{fmt_n(motus_base['motu'].n_unique())}/{fmt_n(len(referenced_motus))} referenced motus")

motus_not_found = referenced_motus - set(motus_base["motu"].unique().to_list())
if motus_not_found:
    log_msg(f"{fmt_n(len(motus_not_found))}/{fmt_n(len(referenced_motus))} motu_match values "
            f"were never found in motus_base (e.g. {', '.join(list(motus_not_found)[:10])}) -- "
            f"rows for these will be absent from the output entirely (no sample data exists "
            f"for them), not present-with-zero.")

#####
# 4. Load genome_summary and identify (Species, Sample) pairs where a MAG was
#    reconstructed -- a Species-level fact (via dRep Species Cluster, the same ID space as
#    our Species column) independent of which motu_match a given output row concerns.
#####
log_msg(f"Reading genome summary: {genome_summary_url}")
genome_summary = pl.read_csv(genome_summary_url)
mag_recon = (
    genome_summary
    .filter(pl.col("Genome").str.contains(r"METAG_[A-Z]{8}$"))
    .with_columns([
        pl.col("Genome").str.replace(r"METAG.*", "METAG").alias("Sample"),
        pl.col("dRep Species Cluster").alias("Species"),
    ])
    .filter(pl.col("Sample").is_in(list(valid_samples)))
    .select("Species", "Sample")
    .unique()
    .with_columns(pl.lit(True).alias("is_mag_reconstructed"))
)
log_msg(f"mag_recon: {fmt_n(mag_recon.height)} (Species, Sample) pairs with a reconstructed MAG")

#####
# 5. Collapse from (Species, motu_match, Sample) down to (Species, Sample): mOTUs are only
#    a proxy for estimating per-species coverage, not a reporting grain in their own right.
#    Taking the MAX across a species' matched motus per sample -- rather than summing --
#    avoids double-counting reads that split across near-duplicate reference clusters (the
#    "majority vote after agreg n=2+" case), and matches the union/OR logic the source
#    script converges on for its own final sample counts (unique samples, not summed
#    per-motu counts): a sample is relevant if ANY of the species' candidate motus shows
#    coverage there.
#####
per_motu = species_motu.join(motus_base, left_on="motu_match", right_on="motu", how="inner")

vertical_coverage = (
    per_motu
    .group_by("Species", "Sample")
    .agg(pl.col("vertical_coverage").max().alias("vertical_coverage"))
)

result = (
    vertical_coverage
    .join(mag_recon, on=["Species", "Sample"], how="left")
    .with_columns(pl.col("is_mag_reconstructed").fill_null(False))
    .select("Species", "Sample", "vertical_coverage", "is_mag_reconstructed")
)

log_msg(f"Output: {result['Species'].n_unique()} species, {result['Sample'].n_unique()} "
        f"samples, {fmt_n(result.height)} rows")

#####
# 6. Write outputs
#####
log_msg(f"Writing output: {out}")
result.write_csv(out, separator="\t")
log_msg("Done.")
