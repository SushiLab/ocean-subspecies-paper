import argparse
import gzip
import os
import shutil
import subprocess
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

PIGZ_THREADS = 4
PIGZ_AVAILABLE = shutil.which("pigz") is not None
if not PIGZ_AVAILABLE:
    logging.warning("pigz not found on PATH; falling back to standard gzip (slower for large files).")


def read_tsv(path: str) -> pd.DataFrame:
    """Read TSV file (gzipped or plain), using pigz for parallel decompression if available."""
    if not path.endswith(".gz"):
        return pd.read_csv(path, sep="\t")

    if PIGZ_AVAILABLE:
        proc = subprocess.Popen(
            ["pigz", "-dc", "-p", str(PIGZ_THREADS), path],
            stdout=subprocess.PIPE,
        )
        try:
            df = pd.read_csv(proc.stdout, sep="\t")
        finally:
            proc.stdout.close()
            ret = proc.wait()
            if ret != 0:
                raise RuntimeError(f"pigz failed to decompress {path} (exit code {ret})")
        return df
    else:
        with gzip.open(path, "rt") as f:
            return pd.read_csv(f, sep="\t")


def write_tsv(df: pd.DataFrame, path: str):
    """Write gzipped TSV file, using pigz for parallel compression if available."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    if PIGZ_AVAILABLE:
        with open(path, "wb") as f_out:
            proc = subprocess.Popen(
                ["pigz", "-c", "-p", str(PIGZ_THREADS)],
                stdin=subprocess.PIPE,
                stdout=f_out,
            )
            try:
                df.to_csv(proc.stdin, sep="\t", index=False)
            finally:
                proc.stdin.close()
                ret = proc.wait()
                if ret != 0:
                    raise RuntimeError(f"pigz failed to compress output to {path} (exit code {ret})")
    else:
        with gzip.open(path, "wt") as f:
            df.to_csv(f, sep="\t", index=False)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--coverage", required=True, help="Merged vertical + horizontal coverage file per species/sample.")
    p.add_argument("--phaseability", required=True, help="File containing phaseability score of each species and sample.")
    p.add_argument("--subspecies_summary", required=True, help="Populations summary file (per species/sample cluster info).")
    p.add_argument("-o", "--output_sample_list", required=True, help="Final Species/Subspecies/Sample table.")

    p.add_argument("--min_vertical_cov", type=float, default=5.0, help="Minimum vertical coverage.")
    p.add_argument("--min_horizontal_cov", type=float, default=0.9, help="Minimum frac_core_detected (horizontal coverage).")
    p.add_argument("--min_samples_vh", type=int, default=20, help="Minimum samples per species passing coverage filter.")
    p.add_argument("--min_phaseability", type=float, default=0.8, help="Minimum phaseability value.")
    p.add_argument("--min_samples_phaseable", type=int, default=15, help="Minimum phaseable samples per species.")
    p.add_argument("--min_umap_trustworthiness", type=float, default=0.6, help="Minimum UMAP trustworthiness per species.")
    p.add_argument("--min_dbcv", type=float, default=0.4, help="Minimum DBCV score per species.")
    p.add_argument("--min_subspecies", type=int, default=2, help="Minimum distinct subspecies (Cluster != -1) per species.")
    return p.parse_args()


def unique_or_first(series: pd.Series, label: str, species: str):
    """Mirrors dplyr's summarize(unique(x)) -- warns if the group isn't actually constant."""
    vals = series.dropna().unique()
    if len(vals) > 1:
        logging.warning(f"Species {species}: expected constant '{label}' but found {len(vals)} distinct values; using first.")
    return vals[0] if len(vals) else pd.NA


def main():
    args = parse_args()
    Path(args.output_sample_list).parent.mkdir(parents=True, exist_ok=True)

    # ============================================================
    # Section 1: filter based on species-sample coverage before SNV calling
    # ============================================================
    coverage_df = read_tsv(args.coverage)

    vh_filt = coverage_df.loc[
        (coverage_df["vertical_coverage"] >= args.min_vertical_cov)
        & (coverage_df["frac_core_detected"] >= args.min_horizontal_cov)
    ].copy()
    logging.info(f"Coverage filter (vertical>={args.min_vertical_cov} & horizontal>={args.min_horizontal_cov}): {len(vh_filt)} species/sample rows retained")

    species_for_snvs = (
        vh_filt.groupby("Species").size().reset_index(name="n_samples_with_vh_coverage")
    )
    species_for_snvs = species_for_snvs[species_for_snvs["n_samples_with_vh_coverage"] >= args.min_samples_vh]
    logging.info(f"Species passing >= {args.min_samples_vh} samples with vh coverage: {len(species_for_snvs)}")

    species_samples_for_snvs = vh_filt[vh_filt["Species"].isin(species_for_snvs["Species"])].copy()

    # ============================================================
    # Section 2: filter based on SNV coverage and phaseability before population pipeline
    # ============================================================
    phaseability_df = read_tsv(args.phaseability).rename(
        columns={"species": "Species", "sample_id": "Sample", "value": "Phaseability"}
    )[["Species", "Sample", "Phaseability"]].drop_duplicates()

    per_sample = species_samples_for_snvs.merge(phaseability_df, on=["Species", "Sample"], how="left")
    per_sample["samples_with_SNV_coverage"] = per_sample["Phaseability"].notna()
    per_sample["samples_phaseable"] = per_sample["samples_with_SNV_coverage"] & (
        per_sample["Phaseability"] >= args.min_phaseability
    )

    per_sample_with_counts = per_sample.merge(species_for_snvs, on="Species", how="left")

    species_for_pop = (
        per_sample_with_counts.groupby("Species")
        .agg(
            n_samples_with_vh_coverage=("n_samples_with_vh_coverage", "first"),
            n_samples_with_SNV_coverage=("samples_with_SNV_coverage", "sum"),
            n_samples_phaseable=("samples_phaseable", "sum"),
        )
        .reset_index()
    )
    species_for_pop = species_for_pop[species_for_pop["n_samples_phaseable"] >= args.min_samples_phaseable]
    logging.info(f"Species passing >= {args.min_samples_phaseable} phaseable samples: {len(species_for_pop)}")

    species_samples_for_pop = per_sample[
        per_sample["samples_with_SNV_coverage"] & per_sample["Species"].isin(species_for_pop["Species"])
    ].copy()

    # ============================================================
    # Section 3: filter based on subspecies assignment & confidence for summary
    # ============================================================
    populations_summary_df = read_tsv(args.subspecies_summary)

    species_sample_with_pop = species_samples_for_pop.merge(
        populations_summary_df, on=["Species", "Sample"], how="left"
    )

    quality_filt = species_sample_with_pop[
        (species_sample_with_pop["umap_trustworthiness"] >= args.min_umap_trustworthiness)
        & (species_sample_with_pop["DBCV"] >= args.min_dbcv)
    ].copy()

    quality_filt_with_counts = quality_filt.merge(species_for_pop, on="Species", how="left")

    def summarize_species(group: pd.DataFrame) -> pd.Series:
        species = group.name
        n_subspecies = group.loc[group["Cluster"] != -1, "Cluster"].nunique()
        return pd.Series({
            "n_samples_with_vh_coverage": unique_or_first(group["n_samples_with_vh_coverage"], "n_samples_with_vh_coverage", species),
            "n_samples_with_SNV_coverage": group["samples_with_SNV_coverage"].sum(),
            "n_samples_phaseable": group["samples_phaseable"].sum(),
            "umap_trustworthiness": unique_or_first(group["umap_trustworthiness"], "umap_trustworthiness", species),
            "DBCV": unique_or_first(group["DBCV"], "DBCV", species),
            "n_subspecies": n_subspecies,
        })

    species_with_subspecies = (
        quality_filt_with_counts.groupby("Species").apply(summarize_species, include_groups=False).reset_index()
    )
    species_with_subspecies = species_with_subspecies[species_with_subspecies["n_subspecies"] > (args.min_subspecies - 1)]
    logging.info(f"Species passing subspecies quality filter (n_subspecies > {args.min_subspecies - 1}): {len(species_with_subspecies)}")

    # Final per-sample table: species_samples_for_pop + populations_summary, filtered to qualifying species
    species_sample_with_subspecies = species_sample_with_pop[
        species_sample_with_pop["Species"].isin(species_with_subspecies["Species"])
    ].copy()

    species_sample_with_subspecies["Subspecies"] = (
        "Subspecies_" + species_sample_with_subspecies["Cluster"].astype("Int64").astype(str)
    )

    final = (
        species_sample_with_subspecies[["Species", "Subspecies", "Sample"]]
        .drop_duplicates()
        .sort_values(["Species", "Subspecies", "Sample"])
    )
    write_tsv(final, args.output_sample_list)
    logging.info(f"Wrote final sample list to {args.output_sample_list}")


if __name__ == "__main__":
    main()