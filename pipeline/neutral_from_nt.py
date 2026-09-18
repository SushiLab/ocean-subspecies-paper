#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script contains function to process the anvio nucleotide variability profile to
1 - produce a between sample neutral distances, which means
    - intergenic regions
    - first and last 10% of genes
"""

import pandas
import argparse
import sklearn.metrics

from umap import validation

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-i', '--input_file', help = 'The anvio nucleotide profile to work with (needs to be produced using the --quince-mode paramter).', required = True, type = str)
    # return namespace
    return parser.parse_args()


def from_lt_table_to_dist(table, output_file):
    """
    Given a table and an output file name,
    compute the nt variab profile distances
    """
    # transform table for distance computation
    table_melt = table.melt(id_vars=["unique_pos_identifier", "sample_id"])

    # reshape the table for distances
    table_unstack = table_melt.set_index(["unique_pos_identifier", "variable", "sample_id"]).unstack()
    table_for_dist = table_unstack.T

    # Get distances and divide by the number of variable positions --> FIXME replace by the genome length ?
    dist_nucl = sklearn.metrics.pairwise.nan_euclidean_distances(table_for_dist)
    print(f"shape of the transformed table for distances: {table_for_dist.shape}")
    unique_positions = set(table.unique_pos_identifier.to_list())
    print(f"number of unique positions {len(unique_positions)}")
    dist_nucl_norm = dist_nucl / len(unique_positions) 

    # Save to tsv
    samples_list = table_for_dist.index.get_level_values(level = "sample_id").tolist()
    pandas.DataFrame(dist_nucl_norm, index = samples_list, columns = samples_list).to_csv(output_file, sep = "\t")


def get_neutral_dist(input_file, min_cov = 3, pos_cov = 3, pos_freq = 0.9, gene_fraction = 0.1):
    """
    This function takes an anvio variability profile computed with the NT engine and using the quince-mode
    """
    print("\n0 - Loading and preparing nucleotide profile.")

    # Load the data
    print(f"Loading {input_file}...")
    table = pandas.read_csv(input_file, sep = "\t")

    # Only keep relevant columns for simplicity
    table_lt = table[["unique_pos_identifier", "sample_id", "codon_order_in_gene", "gene_length", "coverage", "A", "C", "G", "T"]]

    # Check coverage and remove what is below threshold
    # Note that anvio only looked at positions above a given threshold (10X by default)
    # But quince-mode retrieves values for a given position in all samples.
    table_coverage = table_lt[["sample_id", "coverage"]]
    table_coverage_med = table_coverage.groupby(["sample_id"]).median().sort_values(by = "coverage")
    samples_whitelist = table_coverage_med[table_coverage_med.coverage >= min_cov].index.tolist()
    table_lt.set_index("sample_id", inplace = True, drop = False)
    table_lt_filt_s = table_lt.loc[samples_whitelist]
    print(f"Out of the {len(table_coverage_med.index)} samples, {len(samples_whitelist)} went through the {min_cov}X threshold.")

    # Filter for core positions
    table_positions = table_lt_filt_s[table_lt_filt_s.coverage >= pos_cov][['unique_pos_identifier', 'sample_id']].groupby('unique_pos_identifier').count()
    positions_whitelist = table_positions[table_positions.sample_id / len(samples_whitelist) >= pos_freq].index.to_list()
    table_lt_filt_s.set_index("unique_pos_identifier", inplace = True, drop = False)
    table_lt_filt_sp = table_lt_filt_s.loc[positions_whitelist]
    print(f"Out of the {table_positions.shape[0]} variable positions, {len(positions_whitelist)} were found in >={pos_freq} of the {len(samples_whitelist)} samples with coverage >={pos_cov}X and considered core positions.")

    # Get nucleotide freq
    table_lt_filt_sp.loc[:,"A"] = table_lt_filt_sp.loc[:,"A"] / table_lt_filt_sp.coverage
    table_lt_filt_sp.loc[:,"C"] = table_lt_filt_sp.loc[:,"C"] / table_lt_filt_sp.coverage
    table_lt_filt_sp.loc[:,"G"] = table_lt_filt_sp.loc[:,"G"] / table_lt_filt_sp.coverage
    table_lt_filt_sp.loc[:,"T"] = table_lt_filt_sp.loc[:,"T"] / table_lt_filt_sp.coverage

    # filter tables for different neutral subsets
    table_interg = table_lt_filt_sp.loc[table_lt_filt_sp.codon_order_in_gene == -1]

    table_gene_startend = table_lt_filt_sp.loc[table_lt_filt_sp.codon_order_in_gene != -1]
    table_gene_startend.insert(4, "n_codons", [int(i/3) for i in table_gene_startend.gene_length])
    table_gene_startend = table_gene_startend.loc[[x <= gene_fraction*y or x >= (1-gene_fraction)*y for x,y in zip(table_gene_startend.codon_order_in_gene, table_gene_startend.n_codons)]]

    # Compute and save distances
    from_lt_table_to_dist(table_interg.drop(columns = ["codon_order_in_gene", "gene_length", "coverage"]), input_file.replace(".tsv", "-INTERG-DIST.tsv"))
    from_lt_table_to_dist(table_gene_startend.drop(columns = ["codon_order_in_gene", "gene_length", "n_codons", "coverage"]), input_file.replace(".tsv", "-STARTEND-DIST.tsv"))


if __name__ == '__main__':
    args = get_arguments()
    get_neutral_dist(args.input_file)
