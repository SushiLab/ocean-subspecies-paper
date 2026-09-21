#!/usr/bin/env python3
"""
filter_gene_catalog.py

Four-step pipeline:
  1. Read the Species/Subspecies/Sample table -> unique Species list,
     and a Sample -> (Species, Subspecies) lookup.
  2. Read the Species/gene/representative/pangenome table -> for each
     species of interest (from step 1), the set of valid Gene_catalog_rep
     IDs that actually belong to that species' pangenome.
  3. Stream the (large) Gene_catalog_rep/Sample table, keeping only rows
     where BOTH: the Sample belongs to a species of interest, AND the
     Gene_catalog_rep is a valid representative for that *same* species
     (from step 2) - not just any species of interest. Writes a combined
     table: Gene_catalog_rep, Sample, Species, Subspecies. Also collects
     the set of matched Gene_catalog_rep IDs for step 4.
  4. Stream the FASTA file, writing out only the records whose header ID
     is in the set from step 3.

Handles multi-line FASTA and .gz inputs transparently. Only the small
Sample->Species/Subspecies lookup, the per-species valid-representative
sets, and the set of matched gene IDs are held in memory -- the large
gene/sample TSV and the FASTA file are both streamed, and the filtered
catalog table is written directly to disk as it's built.

Usage:
    python3 filter_gene_catalog.py \\
        --species_to_sample species_sample.tsv \\
        --species_gene_rep_map species_gene_representative.tsv \\
        --gene_rep_to_sample gene_sample.tsv \\
        --fasta_in genes.fasta \\
        --fasta_out output.fasta \\
        --filtered_catalog_out filtered_gene_catalog.tsv \\
        [--species-list-out species.txt]

Column names are hardcoded to match the source tables exactly (all files
must have a header row with these names, tab-separated):
    species_to_sample    : Species    Subspecies    Sample
    species_gene_rep_map : Species    gene    representative    pangenome
    gene_rep_to_sample   : Gene_catalog_rep    Sample
"""

import sys
import gzip
import argparse
import time

# Hardcoded column names -- must match the header row of each input file exactly.
SPECIES_COL = "Species"
SUBSPECIES_COL = "Subspecies"
SPECIES_TABLE_SAMPLE_COL = "Sample"

REP_MAP_SPECIES_COL = "Species"
REP_MAP_REP_COL = "representative"

GENE_COL = "Gene_catalog_rep"
GENE_TABLE_SAMPLE_COL = "Sample"

# Subspecies labels ending in "_-1" (e.g. "Subspecies_-1") mark clustering
# noise / unassigned samples, not a genuine subspecies population, and must
# be excluded before any prevalence/grouping work downstream.
NOISE_SUBSPECIES_SUFFIX = "_-1"


def is_noise_subspecies(subspecies):
    return subspecies.endswith(NOISE_SUBSPECIES_SUFFIX)


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--species_to_sample", required=True,
                     help="Species / Subspecies / Sample table")
    ap.add_argument("--species_gene_rep_map", required=True,
                     help="Species / gene / representative / pangenome table")
    ap.add_argument("--gene_rep_to_sample", required=True,
                     help="Gene_catalog_rep / Sample table (large)")
    ap.add_argument("--fasta_in", required=True,
                     help="Input FASTA file with gene sequences")
    ap.add_argument("--fasta_out", required=True,
                     help="Output FASTA file (matches only)")
    ap.add_argument("--filtered_catalog_out", required=True,
                     help="Output TSV: Gene_catalog_rep, Sample, Species, Subspecies "
                          "for the matched rows")
    ap.add_argument("--species-list-out", default=None,
                     help="Optional: write the unique species list here")
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


def load_sample_to_species(path):
    """Step 1: unique species + Sample -> list of (Species, Subspecies).
    A Sample is NOT unique in this table - a metagenomic sample commonly
    hosts more than one species, so the same Sample can legitimately
    appear on multiple rows here, once per species detected in it. We
    must keep all of them, not just the last one seen.

    Rows whose Subspecies is clustering noise (ends in "_-1") are
    excluded entirely - they are not a genuine subspecies population and
    must not be treated as one, nor counted toward a species being "of
    interest" if that's the only kind of row it has."""
    sample_to_species_list = {}
    species_set = set()
    n_total = 0
    n_noise = 0
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = get_col_indices(
            header, [SPECIES_COL, SUBSPECIES_COL, SPECIES_TABLE_SAMPLE_COL], path
        )
        sp_i = idx[SPECIES_COL]
        sub_i = idx[SUBSPECIES_COL]
        sample_i = idx[SPECIES_TABLE_SAMPLE_COL]
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            fields = line.split("\t")
            species = fields[sp_i]
            subspecies = fields[sub_i]
            sample = fields[sample_i]
            if is_noise_subspecies(subspecies):
                n_noise += 1
                continue
            species_set.add(species)
            sample_to_species_list.setdefault(sample, []).append((species, subspecies))
    return sample_to_species_list, species_set, n_total, n_noise


def load_species_to_valid_reps(path, species_set):
    """Step 2: for each species of interest, the set of Gene_catalog_rep
    IDs that are genuinely part of that species' pangenome. This is what
    lets step 3 confirm a detected gene rep actually belongs to the same
    species as the sample it was detected in, rather than just trusting
    that any hit in a sample-of-interest must be relevant."""
    species_to_reps = {}
    n_total = 0
    n_kept = 0
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = get_col_indices(header, [REP_MAP_SPECIES_COL, REP_MAP_REP_COL], path)
        sp_i = idx[REP_MAP_SPECIES_COL]
        rep_i = idx[REP_MAP_REP_COL]
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            fields = line.split("\t")
            species = fields[sp_i]
            if species not in species_set:
                continue
            rep = fields[rep_i]
            species_to_reps.setdefault(species, set()).add(rep)
            n_kept += 1
    return species_to_reps, n_total, n_kept


def filter_and_write_catalog(path, sample_to_species_list, species_to_reps, out_path):
    """Step 3: stream the large gene/sample TSV. For each row, check EVERY
    (species, subspecies) that this sample is associated with (there can be
    more than one), and keep the row once per species for which this gene
    is a valid representative. A single input row can therefore produce
    zero, one, or multiple output rows."""
    keep_ids = set()
    species_output = set()
    n_total = 0
    n_sample_matched = 0
    n_matched = 0
    with opener(path) as fh, open(out_path, "w") as fout:
        header = fh.readline().rstrip("\n").split("\t")
        idx = get_col_indices(header, [GENE_COL, GENE_TABLE_SAMPLE_COL], path)
        gene_i = idx[GENE_COL]
        sample_i = idx[GENE_TABLE_SAMPLE_COL]

        fout.write("Gene_catalog_rep\tSample\tSpecies\tSubspecies\n")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            fields = line.split("\t")
            gene = fields[gene_i]
            sample = fields[sample_i]
            species_subspecies_list = sample_to_species_list.get(sample)
            if not species_subspecies_list:
                continue
            n_sample_matched += 1
            for species, subspecies in species_subspecies_list:
                valid_reps = species_to_reps.get(species)
                if valid_reps and gene in valid_reps:
                    fout.write(f"{gene}\t{sample}\t{species}\t{subspecies}\n")
                    keep_ids.add(gene)
                    species_output.add(species)
                    n_matched += 1
    return keep_ids, species_output, n_total, n_sample_matched, n_matched


def extract_fasta(fasta_path, keep_ids, out_path):
    """Step 4: stream FASTA, write only matching records."""
    n_seen = 0
    n_written = 0
    keep = False
    buf = []
    with opener(fasta_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if buf:
                    fout.writelines(buf)
                    buf.clear()
                n_seen += 1
                header_id = line[1:].split()[0].rstrip("\n")
                keep = header_id in keep_ids
                if keep:
                    n_written += 1
                    buf.append(line)
            elif keep:
                buf.append(line)
        if buf:
            fout.writelines(buf)
    return n_seen, n_written


def main():
    args = parse_args()

    # Step 1
    t0 = time.time()
    print(f"[1/4] Loading species/sample table from {args.species_to_sample} ...",
          file=sys.stderr)
    sample_to_species_list, species_set, n_species_rows, n_noise = load_sample_to_species(args.species_to_sample)
    print(f"      -> {len(species_set):,} unique species, "
          f"{len(sample_to_species_list):,} unique samples in {time.time() - t0:.1f}s",
          file=sys.stderr)
    print(f"      -> {n_species_rows:,} rows scanned, {n_noise:,} excluded as noise "
          f"(Subspecies ending in '{NOISE_SUBSPECIES_SUFFIX}')", file=sys.stderr)
    if args.species_list_out:
        with open(args.species_list_out, "w") as f:
            for sp in sorted(species_set):
                f.write(sp + "\n")
        print(f"      -> species list written to {args.species_list_out}", file=sys.stderr)

    # Step 2
    t1 = time.time()
    print(f"[2/4] Loading species/gene/representative map from "
          f"{args.species_gene_rep_map} ...", file=sys.stderr)
    species_to_reps, n_rep_rows, n_rep_kept = load_species_to_valid_reps(
        args.species_gene_rep_map, species_set
    )
    n_species_with_reps = len(species_to_reps)
    n_total_reps = sum(len(v) for v in species_to_reps.values())
    print(f"      -> {n_rep_rows:,} rows scanned, {n_rep_kept:,} kept for species of "
          f"interest ({n_species_with_reps:,}/{len(species_set):,} species matched, "
          f"{n_total_reps:,} unique representative IDs) in {time.time() - t1:.1f}s",
          file=sys.stderr)
    missing_species = species_set - set(species_to_reps.keys())
    if missing_species:
        print(f"      WARNING: {len(missing_species):,} species of interest had NO "
              f"representatives in {args.species_gene_rep_map} - any of their samples "
              f"will contribute zero rows below.", file=sys.stderr)

    # Step 3
    t2 = time.time()
    print(f"[3/4] Scanning {args.gene_rep_to_sample} for matching samples + reps ...",
          file=sys.stderr)
    keep_ids, species_output, n_total, n_sample_matched, n_matched = filter_and_write_catalog(
        args.gene_rep_to_sample, sample_to_species_list, species_to_reps, args.filtered_catalog_out
    )
    print(f"      -> {n_total:,} rows scanned, {n_sample_matched:,} matched a sample of "
          f"interest, {n_matched:,} output rows written after per-species rep matching "
          f"({len(keep_ids):,} unique gene catalog reps) in {time.time() - t2:.1f}s",
          file=sys.stderr)
    print(f"      -> {len(species_output):,}/{len(species_set):,} species of interest "
          f"actually appear in the filtered catalog", file=sys.stderr)
    unrepresented = species_set - species_output
    if unrepresented:
        print(f"      -> {len(unrepresented):,} species of interest have ZERO rows in "
              f"the output - e.g.: {sorted(unrepresented)[:10]}", file=sys.stderr)
    print(f"      -> filtered catalog written to {args.filtered_catalog_out}",
          file=sys.stderr)

    # Step 4
    t3 = time.time()
    print(f"[4/4] Scanning {args.fasta_in} for matching sequences ...", file=sys.stderr)
    n_seen, n_written = extract_fasta(args.fasta_in, keep_ids, args.fasta_out)
    print(f"      -> {n_seen:,} FASTA records scanned, {n_written:,} written "
          f"in {time.time() - t3:.1f}s", file=sys.stderr)

    print(f"Done. Filtered catalog: {args.filtered_catalog_out}", file=sys.stderr)
    print(f"Done. Output FASTA: {args.fasta_out}", file=sys.stderr)


if __name__ == "__main__":
    main()