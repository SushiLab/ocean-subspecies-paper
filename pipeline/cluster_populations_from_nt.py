#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script contains function to process the anvio nucleotide variability profile to
1 - produce a betweem sample distance matrix
2 - store the sample phaseability
3 - get a UMAP embedding
4 - and HDBSCAN clustering of the data (i.e. populations)
"""

import umap
import pandas
import random
import pickle
import hdbscan
import warnings
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
    # OPTIONAL arguments:
    parser.add_argument("-t", "--threads", help="Number of threads to use.", required=False, default=8, type=int)
    # return namespace
    return parser.parse_args()


def process_anvio_nucleotide_variability(input_file, min_cov = 3, pos_cov = 3, pos_freq = 0.9, pos_phase = 0.9, samples_phase = 0.8, min_phaseable_samples = 15):
    """
    This function takes an anvio variability profile computed with the NT engine and using the quince-mode
    the default coverage paramters were adjusted from metaSNV v1 to capture samples with MAGs, can be as low as 3X.
    the default phaseability parameters are based on metaSNV v2:
    For each species, a ‘discovery subset’ of metagenomes is selected wherein the species is abundant and its population likely contains a single subspecies.
    The latter criterium is satisfied if a metagenome contains minimal internal allele variation relative to the SNV variation across all sampled metagenomes
    (e.g. at least 80% of the dataset-wide species SNVs have the same allele in over 90% of reads in a metagenome).
    This criterium has been previously used (Costea et al., 2017b), and since conceptualized as ‘quasi-phaseability’ (Garud et al., 2019). 
    """
    print("\n0 - Loading and preparing nucleotide profile.")

    # Load the data
    print(f"Loading {input_file}...")
    table = pandas.read_csv(input_file, sep = "\t")

    # Only keep relevant columns for simplicity
    table_lt = table[["unique_pos_identifier", "sample_id", "coverage", "A", "C", "G", "T"]]

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

    # transform table for distance computation
    table_melt = table_lt_filt_sp.drop(columns = ["coverage"]).melt(id_vars=["unique_pos_identifier", "sample_id"])

    # reshape the table for distances
    table_unstack = table_melt.set_index(["unique_pos_identifier", "variable", "sample_id"]).unstack()
    table_for_dist = table_unstack.T

    # Get distances and divide by the number of variable positions --> FIXME replace by the genome length ?
    dist_nucl = sklearn.metrics.pairwise.nan_euclidean_distances(table_for_dist)
    dist_nucl_norm = dist_nucl / len(positions_whitelist)

    # Save to tsv
    samples_list = table_for_dist.index.get_level_values(level = "sample_id").tolist()
    pandas.DataFrame(dist_nucl_norm, index = samples_list, columns = samples_list).to_csv(input_file.replace(".tsv", "-DIST.tsv"), sep = "\t")

    # Add phaseability info
    print(f"Preparing phaseability data with minimum phase frequency of {pos_phase}. Reporting % of positions above the thershold per sample.")
    table_phase = table_melt.groupby(["sample_id", "unique_pos_identifier"]).max("value")
    table_phase_position = table_phase[table_phase.value >= pos_phase].groupby("sample_id").count() / len(positions_whitelist)
    table_phase_position.to_csv(input_file.replace(".tsv", "-PHASEABILITY.tsv"), sep = "\t")

    # Filter profiles for UMAP
    phaseable_samples = table_phase_position.loc[table_phase_position.value >= samples_phase].index.to_list()
    unphaseable_samples = table_phase_position.loc[table_phase_position.value < samples_phase].index.to_list()
    table_for_dist_phaseable = table_for_dist[table_for_dist.index.get_level_values('sample_id').isin(phaseable_samples)]
    table_for_dist_unphaseable = table_for_dist[table_for_dist.index.get_level_values('sample_id').isin(unphaseable_samples)]
    print(f"Using a phaseability threshold of {samples_phase}, {table_for_dist_phaseable.shape[0]} samples are phaseable and {table_for_dist_unphaseable.shape[0]} samples are not.")
    if table_for_dist_phaseable.shape[0] < min_phaseable_samples:
        print(f"Number of phaseable samples ({table_for_dist_phaseable.shape[0]}) is below the threshold ({min_phaseable_samples}), won't run population clutering.")
        return None, None, None

    return table_for_dist_phaseable, table_for_dist_unphaseable, dist_nucl_norm 


def get_umap_embedding(processed_nucl_profile, file_prefix, n_epochs=1000, n_iter=30, samples_frac=0.1, K_weight=3, min_trust=0.9, n_jobs=32):
    """
    This fucntion takes the processed variability profile and gets the umap embedding
    with some hyperparameter space exploration
    """

    def wrapper_umap_on_nt_profile(df, params, K, **kwargs):
        """
        Function to run umap on a nt profile df with a given set of hyper parameters
        returns the umap output and the trustworthiness of the embedding
        Note that as opposed to precomputed distances that can handle NaN values, this option doesn't and NaN are replaced by 0.
        The lower the pos_freq filter, the more problematic this becomes
        """
        # DEBUG
        print(f"Iterating umap with params {params}.")
        umap_out = umap.UMAP(n_neighbors=params["n_neighbors"], min_dist=params["min_dist"], metric="euclidean", **kwargs).fit(df.fillna(0).to_numpy())
        umap_trust = sum(validation.trustworthiness_vector(source=df.fillna(0).to_numpy(), embedding=umap_out .embedding_, max_k=K))/K
        return {"params": params, "umap": umap_out, "trust": umap_trust} 


    def wrapper_umap_on_dist(dist_array, params, K, **kwargs):
        """
        Function to run umap on a nt profile df with a given set of hyper parameters
        returns the umap output and the trustworthiness of the embedding
        left as legacy, but the other option should be best
        """
        umap_out = umap.UMAP(n_neighbors=params["n_neighbors"], min_dist=params["min_dist"], metric="precomputed", **kwargs).fit(dist_array)
        umap_trust = sum(validation.trustworthiness_vector(source=dist_array, embedding=umap_out .embedding_, max_k=K))/K
        return {"params": params, "umap": umap_out, "trust": umap_trust} 


    # Get embedding
    print("\n1 - Embedding the data")
    # Init (hyper)parameters
    init_n_neighbors = round(samples_frac*processed_nucl_profile.shape[0])
    K = K_weight*init_n_neighbors
    print(f"Using frac:{samples_frac}, n_iter:{n_iter}, n_epochs:{n_epochs} and K:{K}.")

    # Explore hyper parameter space
    umap_params = {'n_neighbors' : list({2, 5, 10, round(init_n_neighbors/2), init_n_neighbors, 2*init_n_neighbors, 3*init_n_neighbors}),
                   'min_dist'    : [0, 10**-5, 10**-3, 10**-2, 10**-1]}  
    # safety, since umap cannot handle 1 neighbor
    if 0 in umap_params["n_neighbors"]: umap_params["n_neighbors"].remove(0)
    if 1 in umap_params["n_neighbors"]: umap_params["n_neighbors"].remove(1)
    print(f"Exploring embedding space with params:\n{umap_params}")

    umap_explo = []
    for i in range(n_iter):
        params_iter = {"n_neighbors" : random.choice(umap_params["n_neighbors"]),
                       "min_dist"    : random.choice(umap_params["min_dist"])}
        umap_explo.append(wrapper_umap_on_nt_profile(processed_nucl_profile, params_iter, K, n_jobs=n_jobs))

    optimal_trust = max([x['trust'] for x in umap_explo])
    if optimal_trust < min_trust:
        #raise Exception(f"Best embedding returns K-average trust of {optimal_trust}, conservatively, this is not good enough.")
        warnings.warn(f"Best embedding returns K-average trust of {optimal_trust}, conservatively, this is not good enough.")
    optimal_umap = [x["umap"] for x in umap_explo if x["trust"] == optimal_trust][0]
    optimal_params = [x["params"] for x in umap_explo if x["trust"] == optimal_trust][0]
    print(f"Found best embedding with K-averaged trustworthiness {optimal_trust} and params {optimal_params}.")

    # Save to file
    with open(f"{file_prefix}-UMAP.pickle", "wb") as handle:
        pickle.dump(optimal_umap, handle)
    with open(f"{file_prefix}-UMAP.trustworthiness", "w") as handle:
        handle.write(f"{optimal_trust}\n")
    umap_df = pandas.DataFrame(optimal_umap.embedding_, columns = ["V1", "V2"])
    umap_df.insert(0, "Sample", processed_nucl_profile.index.get_level_values('sample_id'))
    umap_df.to_csv(f"{file_prefix}-UMAP-PHASEABLE.tsv", sep = "\t", index = False)
        
    return optimal_umap, umap_params 


def get_hdbscan_clustering(optimal_umap, umap_params, file_prefix, samples, min_DBCV = 0.75, n_iter = 50): 
    """
    From the UMAP embedding, cluster the data to identify discrete populations
    """
    print("\n2 - Clustering the final embedding with HDBSCAN*")
    optimal_embedding = optimal_umap.embedding_

    # Random search for parameter tuning
    print("Optimizing HDBSCAN...")
    # specify parameters and distributions to sample from
    hdbscan_params = {'min_samples'      : [None] + umap_params["n_neighbors"], 
                      'min_cluster_size' : umap_params["n_neighbors"]} 
    print(f"Exploring parameters: {hdbscan_params}")

    hdbscan_explo = []
    for i in range(n_iter):
        params_iter = {"min_samples" : random.choice(hdbscan_params["min_samples"]),
                       "min_cluster_size"    : random.choice(hdbscan_params["min_cluster_size"])}
        clusterer_iter = hdbscan.HDBSCAN(gen_min_span_tree=True, prediction_data=True, min_samples=params_iter["min_samples"], min_cluster_size=params_iter["min_cluster_size"]).fit(optimal_embedding)
        hdbscan_explo.append({"params": params_iter, "hdbscan": clusterer_iter, "DBCV": clusterer_iter.relative_validity_})

    optimal_DBCV = max([x['DBCV'] for x in hdbscan_explo])
    if optimal_DBCV < min_DBCV:
        #raise Exception(f"Best HDBSCAN clustering on optimal embedding returns DBCV of {optimal_DBCV}, conservatively, this is not good enough.")
        warnings.warn(f"Best HDBSCAN clustering on optimal embedding returns DBCV of {optimal_DBCV}, conservatively, this is not good enough.")

    optimal_clusterer = [x["hdbscan"] for x in hdbscan_explo if x["DBCV"] == optimal_DBCV][0]
    optimal_params = [x["params"] for x in hdbscan_explo if x["DBCV"] == optimal_DBCV][0]
    print(f"Found best clusterer with DBCV of {optimal_DBCV} and params {optimal_params}.")
    print(f"It identified {max(optimal_clusterer.labels_) + 1} clusters with a coverage of {sum(optimal_clusterer.labels_ >= 0) / optimal_embedding.shape[0]}.")

    # Save to file
    with open(f"{file_prefix}-HDBSCAN.pickle", "wb") as handle:
        pickle.dump(optimal_clusterer, handle)
    with open(f"{file_prefix}-HDBSCAN.dbcv", "w") as handle:
        handle.write(f"{optimal_DBCV}\n")
    hdbscan_df = pandas.DataFrame({"Sample": samples, "Cluster": optimal_clusterer.labels_, "Probability": optimal_clusterer.probabilities_})
    hdbscan_df.to_csv(f"{file_prefix}-HDBSCAN-PHASEABLE.tsv", sep = "\t", index = False)

    return optimal_clusterer


def predict_new_samples(df, optimal_umap, optimal_clusterer, file_prefix):
    """
    This function takes new nucleotide profile and runs them through the learned
    embedding and clustering to predict populations
    """
    print("\n3 - Predicting the clustering of additional (e.g. Non-phaseable) samples.")
    new_samples_embedded = optimal_umap.transform(df.fillna(0))
    new_samples_labels, new_samples_strengths = hdbscan.approximate_predict(optimal_clusterer, new_samples_embedded)
    print(f"Clustered the new samples in {new_samples_labels}.")

    # Save to file
    umap_df = pandas.DataFrame(new_samples_embedded, columns = ["V1", "V2"])
    umap_df.insert(0, "Sample", df.index.get_level_values('sample_id'))
    umap_df.to_csv(f"{file_prefix}-UMAP-UNPHASEABLE.tsv", sep = "\t", index = False)
    hdbscan_df = pandas.DataFrame({"Sample": df.index.get_level_values('sample_id'), "Cluster": new_samples_labels, "Probability": new_samples_strengths})
    hdbscan_df.to_csv(f"{file_prefix}-HDBSCAN-UNPHASEABLE.tsv", sep = "\t", index = False)

    return new_samples_labels, new_samples_strengths


def get_umap_hdbscan_from_pickle(umap_file, hdbscan_file):
    """
    Utils function to load models from pickle files
    """
    with open(umap_file, 'rb') as handle:
        trained_umap = pickle.load(handle)
    with open(hdbscan_file, 'rb') as handle:
        trained_hdbscan = pickle.load(handle)
    return trained_umap, trained_hdbscan


if __name__ == '__main__':
    args = get_arguments()
    phaseable, unphaseable, dist = process_anvio_nucleotide_variability(args.input_file)
    if phaseable is not None:
        optimal_umap, umap_params = get_umap_embedding(phaseable, args.input_file.replace(".tsv", ""), n_jobs=args.threads)
        optimal_clusterer = get_hdbscan_clustering(optimal_umap, umap_params, args.input_file.replace(".tsv", ""), phaseable.index.get_level_values('sample_id'))
        if unphaseable.shape[0] > 0:
            new_samples_labels, new_samples_strengths = predict_new_samples(unphaseable, optimal_umap, optimal_clusterer, args.input_file.replace(".tsv", ""))
