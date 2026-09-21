#!/usr/bin/env python3
"""
subset_species_tables.py

Given a Species label, extract everything relevant to that species from the larger
multi-speices tables and write standalone per-species subset files.

Usage:
    python3 subset_species_tables.py \\
        --species 1064_1 \\
        --species_sample species_samples.tsv \\
        --subspecies_core_genes subspecies_core_gene_list.tsv \\
        --subspecies_core_kegg subspecies_core_kegg_list.tsv \\
        --subspecies_core_pfam subspecies_core_pfam_list.tsv \\
        --species_gene_composition species_gene_composition.tsv \\
        --pnps_per_gene_sample species_pnps.tsv.gz \\
        --gene_catalog_profile gene_catalog_profile.tsv.gz \\
        --snv_distances POPULATIONS-DIST-SUMMARY.tsv.gz \\
        --output_dir subset_output/ \\
        [--output_prefix 1064_1]

Column names are hardcoded to match the source tables exactly (tab-separated,
header row required, .gz inputs handled transparently):
    species_sample            : Species  Subspecies  Sample
    subspecies_core_genes     : Species  Subspecies  Gene_catalog_rep
                                 (also the source of the gene set used to
                                 filter the gene catalog profile below)
    subspecies_core_kegg      : Species  Subspecies  KEGG_ko
    subspecies_core_pfam      : Species  Subspecies  PFAM_accession
    species_gene_composition  : Species  Gene_catalog_rep  KEGG_ko  KEGG_description
                                 PFAM_accession  PFAM_name  PFAM_description
    pnps_per_gene_sample       : Species  ref_genome  anvio_id  gene  gene_catalog_representative
                                 sample  pNpS_gene_reference
                                 (filtered on Species only; the gene ID column
                                 here is 'gene_catalog_representative', not
                                 'Gene_catalog_rep', but that distinction
                                 doesn't matter for this filter)
    gene_catalog_profile      : Gene_catalog_rep  Sample   (large - streamed)
    snv_distances              : species  refg  sample_1  sample_2  dist  cluster_1  cluster_2  comp
                                 (note: lowercase 'species' here, unlike every other table)
"""

import sys
import gzip
import argparse
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--species", required=True, help="Species label to subset, e.g. 1064_1")
    ap.add_argument("--species_sample", required=True, help="Species / Subspecies / Sample table")
    ap.add_argument("--subspecies_core_genes", required=True, help="Species / Subspecies / Gene_catalog_rep table")
    ap.add_argument("--subspecies_core_kegg", required=True, help="Species / Subspecies / KEGG_ko table")
    ap.add_argument("--subspecies_core_pfam", required=True, help="Species / Subspecies / PFAM_accession table")
    ap.add_argument("--species_gene_composition", required=True, help="Species / Gene_catalog_rep / KEGG_ko / KEGG_description / PFAM_accession / PFAM_name / PFAM_description table")
    ap.add_argument("--pnps_per_gene_sample", required=True, help="Species / ref_genome / anvio_id / gene / gene_catalog_representative / sample / pNpS_gene_reference table")
    ap.add_argument("--gene_catalog_profile", required=True, help="Gene_catalog_rep / Sample table (large - streamed)")
    ap.add_argument("--snv_distances", required=True, help="Sample-pair SNV distance table (lowercase 'species' column)")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--output_prefix", default=None, help="Defaults to the species label")
    return ap.parse_args()


def opener(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def get_col_indices(header_fields, required_names, source_path):
    """Locate each required column name in the header, failing loudly and
    clearly if the file doesn't have the expected columns."""
    idx = {}
    for name in required_names:
        if name not in header_fields:
            raise ValueError(
                f"Expected column '{name}' not found in header of {source_path}: "
                f"{header_fields}"
            )
        idx[name] = header_fields.index(name)
    return idx


def filter_small_table(path, species_col, species, out_path, extra_cols_to_collect=None):
    """Load a small table fully, keep rows matching `species` in
    `species_col`, write the filtered table, and optionally return sets of
    values from other columns (e.g. Gene_catalog_rep, Sample) collected in
    the same pass, for use in the large-file streaming step below."""
    extra_cols_to_collect = extra_cols_to_collect or []
    collected = {c: set() for c in extra_cols_to_collect}
    n_total = 0
    n_kept = 0
    with opener(path) as fh, open(out_path, "w") as fout:
        header_line = fh.readline().rstrip("\n")
        header = header_line.split("\t")
        idx = get_col_indices(header, [species_col] + extra_cols_to_collect, path)
        sp_i = idx[species_col]
        fout.write(header_line + "\n")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            fields = line.split("\t")
            if fields[sp_i] != species:
                continue
            n_kept += 1
            fout.write(line + "\n")
            for c in extra_cols_to_collect:
                collected[c].add(fields[idx[c]])
    return n_total, n_kept, collected


def filter_gene_catalog_profile(path, gene_col, sample_col, valid_genes, valid_samples, out_path):
    """Stream the large gene catalog profile, keeping rows whose gene is in
    the species' core gene set AND whose sample is in the species' sample
    set - never loaded into memory as a whole."""
    n_total = 0
    n_kept = 0
    with opener(path) as fh, open(out_path, "w") as fout:
        header_line = fh.readline().rstrip("\n")
        header = header_line.split("\t")
        idx = get_col_indices(header, [gene_col, sample_col], path)
        gene_i = idx[gene_col]
        sample_i = idx[sample_col]
        fout.write(header_line + "\n")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            fields = line.split("\t")
            if fields[gene_i] in valid_genes and fields[sample_i] in valid_samples:
                n_kept += 1
                fout.write(line + "\n")
    return n_total, n_kept


def main():
    args = parse_args()
    species = args.species
    prefix = args.output_prefix or species
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def out(name):
        return output_dir / f"{prefix}.{name}"

    def log(msg):
        print(msg, file=sys.stderr)

    # 1. Species -> Sample -> Subspecies
    log("[1/8] Filtering species/sample table...")
    n_total, n_kept, collected = filter_small_table(
        args.species_sample, "Species", species, out("species_to_sample.tsv"),
        extra_cols_to_collect=["Sample"]
    )
    valid_samples = collected["Sample"]
    log(f"      -> {n_kept:,}/{n_total:,} rows kept, {len(valid_samples):,} unique samples")

    # 2. Subspecies core genes - also the source of valid_genes for step 5,
    # since the species-level core gene list is no longer used at all
    log("[2/8] Filtering subspecies core gene list...")
    n_total, n_kept, collected = filter_small_table(
        args.subspecies_core_genes, "Species", species, out("subspecies_core_genes.tsv"),
        extra_cols_to_collect=["Gene_catalog_rep"]
    )
    valid_genes = collected["Gene_catalog_rep"]
    log(f"      -> {n_kept:,}/{n_total:,} rows kept, {len(valid_genes):,} unique core genes across all subspecies")

    if not valid_samples:
        log(f"      WARNING: no samples found for species '{species}' in {args.species_sample} - "
            f"the gene catalog profile subset below will be empty.")
    if not valid_genes:
        log(f"      WARNING: no core genes found for species '{species}' in {args.subspecies_core_genes} - "
            f"the gene catalog profile subset below will be empty.")

    # 3. Subspecies core KEGG
    log("[3/8] Filtering subspecies core KEGG list...")
    n_total, n_kept, _ = filter_small_table(
        args.subspecies_core_kegg, "Species", species, out("subspecies_core_kegg.tsv")
    )
    log(f"      -> {n_kept:,}/{n_total:,} rows kept")

    # 4. Subspecies core Pfam
    log("[4/8] Filtering subspecies core Pfam list...")
    n_total, n_kept, _ = filter_small_table(
        args.subspecies_core_pfam, "Species", species, out("subspecies_core_pfam.tsv")
    )
    log(f"      -> {n_kept:,}/{n_total:,} rows kept")

    # 5. Species gene composition (KEGG + Pfam annotations per gene)
    log("[5/8] Filtering species gene composition table...")
    n_total, n_kept, _ = filter_small_table(
        args.species_gene_composition, "Species", species, out("species_gene_composition.tsv")
    )
    log(f"      -> {n_kept:,}/{n_total:,} rows kept")

    # 6. pNpS per gene per sample
    log("[6/8] Filtering pNpS-per-gene-per-sample table...")
    n_total, n_kept, _ = filter_small_table(
        args.pnps_per_gene_sample, "Species", species, out("pnps_per_gene_sample.tsv")
    )
    log(f"      -> {n_kept:,}/{n_total:,} rows kept")

    # 7. Gene catalog profile (large - streamed, never loaded fully)
    log("[7/8] Streaming gene catalog profile (large file)...")
    n_total, n_kept = filter_gene_catalog_profile(
        args.gene_catalog_profile, "Gene_catalog_rep", "Sample",
        valid_genes, valid_samples, out("gene_catalog_profile.tsv")
    )
    log(f"      -> {n_kept:,}/{n_total:,} rows kept")

    # 8. SNV distances (note lowercase 'species' column)
    log("[8/8] Filtering SNV distance table...")
    n_total, n_kept, _ = filter_small_table(
        args.snv_distances, "species", species, out("snv_distances.tsv")
    )
    log(f"      -> {n_kept:,}/{n_total:,} rows kept")

    log(f"Done. All subset tables written to: {output_dir} (prefix: '{prefix}.')")


if __name__ == "__main__":
    main()