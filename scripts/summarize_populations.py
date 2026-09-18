#!/usr/bin/env python

"""
This script can be used to parse a 'species_clusters' directory and summarize the populations output
it will look for folders that contain both the "*-POPULATIONS.done" marker and the "*.dbcv" files and produce a summary table with the columns:
    species, sample, cluster, probability, phaseability, umap trustworthiness, dbcv
"""

import pandas
import argparse

from pathlib import Path


def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-i', '--input_dir', help = 'the species_clusters directory to summarize.', required = True, type = str)
    parser.add_argument('-o', '--output_table', help = 'Path and prefix.', required = True, type = str)
    # return namespace
    return parser.parse_args()


def summarize_clusters(input_dir, output_table):
    """
    function to iterate over the species_clusters directory
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise Exception(f"Input path {input_path} doesn't exist.")
    output_path = Path(output_table + "-SUMMARY.tsv")
    if output_path.exists():
        raise Exception(f"Output table {output_path} already exists.")

    markers_list = [i for i in input_path.glob("*/*-POPULATIONS.done")]
    species_refg_dict = {i.parent.name: i.name.replace("-POPULATIONS.done", "") for i in markers_list}

    res = pandas.DataFrame()

    print(f"Found a total of {len(markers_list)} species with the populations marker file, will process those.")
    
    for species, refg in species_refg_dict.items():
        dbcv_file = input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-HDBSCAN.dbcv")
        trust_file = input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-UMAP.trustworthiness")
        if trust_file.exists():
            if not dbcv_file.exists():
                raise Exception(f"Found {trust_file} but not {dbcv_file}, that shouldn't be.")
            
            print(f"Processing species {species}...")

            trust = trust_file.read_text().strip()
            dbcv = dbcv_file.read_text().strip()
            phaseable = pandas.read_csv(input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-HDBSCAN-PHASEABLE.tsv"), sep = "\t") 
            phaseable.insert(1, "Phaseable", True)
            unphaseable_tsv = input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-HDBSCAN-UNPHASEABLE.tsv")
            if unphaseable_tsv.exists(): # May not exist if all samples are phaseable
                unphaseable = pandas.read_csv(unphaseable_tsv, sep = "\t") 
                unphaseable.insert(1, "Phaseable", False)
            else:
                unphaseable = pandas.DataFrame()
            species_df = pandas.concat([phaseable, unphaseable]) 
            species_df.insert(0, "Species", species)
            species_df.insert(1, "ref_genome", refg)
            species_df["umap_trustworthiness"] = trust
            species_df["DBCV"] = dbcv
            res = pandas.concat([res, species_df])
    
            res.to_csv(output_path, sep = "\t", index = False)

        else:
            print(f"Skipping species {species} as it seems population clustering was skipped (most likely not enough phaseable samples).")


def summarize_embeddings(input_dir, output_table):
    """
    function to iterate over the species_clusters directory
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise Exception(f"Input path {input_path} doesn't exist.")
    output_path = Path(output_table + "-EMBEDDINGS.tsv")
    if output_path.exists():
        raise Exception(f"Output table {output_path} already exists.")

    markers_list = [i for i in input_path.glob("*/*-POPULATIONS.done")]
    species_refg_dict = {i.parent.name: i.name.replace("-POPULATIONS.done", "") for i in markers_list}

    res = pandas.DataFrame()

    print(f"Found a total of {len(markers_list)} species with the populations marker file, will process those.")
    
    for species, refg in species_refg_dict.items():
        dbcv_file = input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-HDBSCAN.dbcv")
        trust_file = input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-UMAP.trustworthiness")
        if trust_file.exists():
            if not dbcv_file.exists():
                raise Exception(f"Found {trust_file} but not {dbcv_file}, that shouldn't be.")
            
            print(f"Processing species {species}...")

            trust = trust_file.read_text().strip()
            phaseable = pandas.read_csv(input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-UMAP-PHASEABLE.tsv"), sep = "\t") 
            phaseable.insert(1, "Phaseable", True)
            unphaseable_tsv = input_path.joinpath(species, f"{refg}-ANVIO_VARIABILITY-NT-UMAP-UNPHASEABLE.tsv")
            if unphaseable_tsv.exists(): # May not exist if all samples are phaseable
                unphaseable = pandas.read_csv(unphaseable_tsv, sep = "\t") 
                unphaseable.insert(1, "Phaseable", False)
            else:
                unphaseable = pandas.DataFrame()
            species_df = pandas.concat([phaseable, unphaseable]) 
            species_df.insert(0, "Species", species)
            species_df.insert(1, "ref_genome", refg)
            species_df["umap_trustworthiness"] = trust
            res = pandas.concat([res, species_df])
    
            res.to_csv(output_path, sep = "\t", index = False)

        else:
            print(f"Skipping species {species} as it seems population clustering was skipped (most likely not enough phaseable samples).")


if __name__ == '__main__':
    args = get_arguments()
    summarize_clusters(args.input_dir, args.output_table)
    summarize_embeddings(args.input_dir, args.output_table)
