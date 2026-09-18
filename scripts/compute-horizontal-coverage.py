#!/usr/bin/env python

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
pangenome_file  = "../../data/processed/summaries/PANGENOME-GENE-CATEGORY-FULL.tsv.gz"
detection_dir   = "../../data/processed/summaries/GENE-CATALOG-FULL-DETECTION.parquet"   # directory of part-*.parquet files
                                                            # from preprocess-profiles.py --
                                                            # covers the FULL gene catalog
                                                            # (core+shell+cloud), not just
                                                            # core. Read via a lazy scan with
                                                            # filter pushdown below: Parquet's
                                                            # per-row-group statistics let
                                                            # this skip whole row groups that
                                                            # can't contain a matching
                                                            # Gene_catalog_rep without ever
                                                            # decompressing them.

out = "../../data/processed/summaries/horizontal-coverage-per-sample-per-species.tsv"

#####
# 1. Process Pan table
#####

log_msg(f"Reading pangenome table: {pangenome_file}")
pan = pl.read_csv(pangenome_file, separator="\t")
log_msg(f"pan: {fmt_n(pan.height)} rows loaded")

# PANGENOME-GENE-CATEGORY-FULL.tsv.gz ships as (Species, gene, representative, pangenome)
# rather than (Species, Gene_catalog_rep, Pangenome_group). It's `representative` that maps
# to Gene_catalog_rep here, not `gene`: the profile files are quantified against the
# non-redundant representative catalog (standard practice -- map reads against one
# representative sequence per cluster, not every redundant instance), so `representative`
# is the ID space that actually populates the profile files' rows, confirmed against a real
# run: 623,594 distinct matches in BOTH independent profile files, landing almost exactly on
# the 641,503 total core representatives rather than the 14.8M individual gene instances.
# Using `gene` here would silently cap every frac_core_detected near
# representative_count/instance_count (~4%) regardless of true detection.
pan = pan.rename({"representative": "Gene_catalog_rep", "pangenome": "Pangenome_group"})

# The Species column in pangenome_file is already restricted to the species of interest
# (filtered upstream), so it's used directly here rather than cross-checking against a
# separate species list file.
species_list = pan["Species"].unique().to_list()
log_msg(f"species_list: {fmt_n(len(species_list))} distinct species found in pangenome_file")

core_genes = (
    pan.filter(pl.col("Pangenome_group") == "core")
       .select("Species", "Gene_catalog_rep")
       .unique()
)

# --- representatives shared as core across more than one Species ---
# This is expected biology (two species can legitimately share a core gene), not a data
# error -- each shared representative is deliberately kept once per Species below, so a
# single detected copy counts independently toward every species that claims it as core.
shared_rep_counts = (
    core_genes.group_by("Gene_catalog_rep").agg(pl.len().alias("n")).filter(pl.col("n") > 1)
)
if shared_rep_counts.height > 0:
    examples = shared_rep_counts["Gene_catalog_rep"].head(5).to_list()
    log_msg(f"{fmt_n(shared_rep_counts.height)} core representative(s) are shared as core "
            f"across more than one Species (e.g. {', '.join(examples)}) -- each is counted "
            f"independently for every species that claims it.")

n_core_per_species = (
    core_genes.group_by("Species").agg(pl.col("Gene_catalog_rep").n_unique().alias("n_core"))
)
core_gene_ids = core_genes["Gene_catalog_rep"].unique().to_list()
log_msg(f"core_genes: {fmt_n(core_genes.height)} rows, {fmt_n(len(core_gene_ids))} distinct "
        f"core representatives across {n_core_per_species.height}/{len(species_list)} species")

# --- species present in pan but with zero rows labeled "core" ---
species_no_core = sorted(set(species_list) - set(n_core_per_species["Species"].to_list()))
if species_no_core:
    log_msg(f"{len(species_no_core)} species are in pangenome_file but have no "
            f"Pangenome_group == 'core' genes: {', '.join(species_no_core)}")

assert n_core_per_species.height > 0, \
    "No species have any core genes -- check Pangenome_group values in pangenome_file"

#####
# 2. Lazily scan the precomputed detection dataset (all part-*.parquet files as one logical
#    table), filtering to core genes via predicate pushdown -- Parquet stores per-row-group
#    statistics, so this can skip whole row groups that can't contain a matching
#    Gene_catalog_rep without decompressing them, and the scan never materializes the full
#    ~2+ billion row dataset at any point.
#####
log_msg(f"Scanning detection table: {detection_dir}/*.parquet")

detected_lazy = pl.scan_parquet(f"{detection_dir}/*.parquet")
all_samples_seen = detected_lazy.select("Sample").unique().collect()["Sample"].to_list()

core_genes_lazy = pl.LazyFrame(core_genes)
prof_core_only = (
    detected_lazy
    .filter(pl.col("Gene_catalog_rep").is_in(core_gene_ids))
    .join(core_genes_lazy, on="Gene_catalog_rep", how="inner")   # fans out naturally for
                                                                    # shared representatives,
                                                                    # no special flag needed
                                                                    # (unlike data.table)
)

detected_core = (
    prof_core_only
    .group_by("Species", "Sample")
    .agg(pl.col("Gene_catalog_rep").n_unique().alias("n_detected"))
    .collect()
)

core_genes_detected_seen = set(
    prof_core_only.select("Gene_catalog_rep").unique().collect()["Gene_catalog_rep"].to_list()
)
missing_genes = set(core_gene_ids) - core_genes_detected_seen
if missing_genes:
    species_affected = (
        core_genes.filter(pl.col("Gene_catalog_rep").is_in(list(missing_genes)))["Species"].unique().to_list()
    )
    log_msg(f"{fmt_n(len(missing_genes))}/{fmt_n(len(core_gene_ids))} core Gene_catalog_reps "
            f"({100 * len(missing_genes) / len(core_gene_ids):.1f}%) were never detected in "
            f"any sample, affecting {len(species_affected)} species -- could be a real "
            f"absence, an ID mismatch, or a gene that was always zero.")
# NOTE: detection_dir only ever contains DETECTED pairs (Abundance > 0 already applied
# upstream), so this can no longer distinguish "never appeared as a row in either profile
# file" from "appeared but was zero in every sample" -- both look identical here.

log_msg(f"Detected core genes: {fmt_n(detected_core.height)} Species x Sample groups, "
        f"{fmt_n(len(core_genes_detected_seen))}/{fmt_n(len(core_gene_ids))} core "
        f"Gene_catalog_reps detected")

#####
# 3. Expand to the full Species x Sample grid so species/samples with zero detected core
#    genes appear as frac_core_detected = 0 rather than being silently dropped by the join
#####
species_col = n_core_per_species["Species"].to_list()
full_grid = pl.DataFrame({
    "Species": [s for s in species_col for _ in all_samples_seen],
    "Sample":  all_samples_seen * len(species_col),
})

detected_core = (
    full_grid
    .join(detected_core, on=["Species", "Sample"], how="left")
    .with_columns(pl.col("n_detected").fill_null(0))
    .join(n_core_per_species, on="Species", how="inner")
    .with_columns((pl.col("n_detected") / pl.col("n_core")).alias("frac_core_detected"))
)

species_missing_from_output = set(n_core_per_species["Species"].to_list()) - set(detected_core["Species"].to_list())
assert not species_missing_from_output, "Some species with core genes are still missing from the final output"

log_msg(f"Output: {detected_core['Species'].n_unique()} species x "
        f"{detected_core['Sample'].n_unique()} samples = {detected_core.height} rows.")

#####
# 4. Write outputs
#####
log_msg(f"Writing output: {out}")
detected_core.write_csv(out, separator="\t")
log_msg("Done.")
