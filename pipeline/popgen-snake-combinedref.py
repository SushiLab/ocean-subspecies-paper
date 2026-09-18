#!usr/bin/env python
# -*- coding: utf-8

from pathlib import Path

import Bio.SeqIO.FastaIO as FastaIO

# Variables ====================================================================

working_dir = "/cluster/project/sunagawa/paolil/ocean-popgen"
working_path = Path(working_dir)
genomes_path = Path("/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/integrated/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-genomes")
samples_map = working_path.joinpath('scripts/samples.map').read_text().splitlines()
seqstor_index = working_path.joinpath('scripts/seqstor.index').read_text().splitlines()

input_data = working_path.joinpath('scripts/popgen-processed_input.tsv').read_text().splitlines() 

# Setup ========================================================================

# prepare a dict between old and new sample name dict[old] = new
samples_old2new = {}
samples_new2old = {}
for sample_line in samples_map:
    if not "METAT" in sample_line:
        sample_split = sample_line.split('\t')
        samples_old2new[sample_split[0]] = sample_split[1]
        samples_new2old[sample_split[1]] = sample_split[0]

# prepare a dict between samples and their read files in the seqstor
seqstor_dict = {}
for seqstor_line in seqstor_index:
    seqstor_file = seqstor_line.split('/')[-1]
    if not seqstor_file in seqstor_dict.keys():
        seqstor_dict[seqstor_file] = seqstor_line

# prepare the main input dict with the species, a list of samples (oldnames) and the repr genome
input_dict = {}
flag = 0
for input_line in input_data:
    """
    Expects a table with the format
    species representative_genome   motu    type    n5x     samples
    """
    input_split = input_line.split('\t')
    if input_line.startswith("species"):
        flag += 1
        if input_split[5] != "samples" or input_split[1] != "representative_genome":
            raise Exception("Looks like input table is not in the right format.")

    else:
        species = input_split[0]
        species_samples = input_split[5].split(';')
        species_refgenome = input_split[1].replace(".fa", "")
        input_dict[species] = {'samples': species_samples, 'refgenome': species_refgenome} 
        species_path = working_path.joinpath("species_clusters", species)
        species_path.mkdir(exist_ok = True)
        contigs_list = species_path.joinpath(species_refgenome + "-contigs.txt")
        species_refgenome_path = species_path.joinpath(species_refgenome + ".fa")
        if not species_refgenome_path.exists():
            species_refgenome_path.symlink_to(genomes_path.joinpath(species_refgenome + ".fa"))

        if not contigs_list.exists():
            with open(contigs_list, 'w') as target:
                with open(species_refgenome_path) as handle:
                    for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                        target.write(header + '\n')


if not flag == 1:
    raise Exception("Looks like input table is not in the right format.")

# Define the targets ===========================================================

TARGETS = []

# BAMFILES
#for s in samples_new2old.keys():
#    TARGETS.append(f"{working_dir}/refgenomes/bams/{s}/{s}-vs-COMBINED_REF_GENOMES_FOR_POPGEN.done")

# PROFILES
input_anvio_merge = {}
for key, value in input_dict.items():
    #input_anvio_merge[key] = [f"{working_dir}/species_clusters/{key}/profiles/{samples_old2new[sample]}-vs-{value['refgenome']}.anvio.done" for sample in value['samples']]
    #TARGETS.append(f"{working_dir}/species_clusters/{key}/{value['refgenome']}-MERGED_PROFILE.done")
    TARGETS.append(f"{working_dir}/species_clusters/{key}/{value['refgenome']}-GENES_CONSENSUS.done")
    #TARGETS.append(f"{working_dir}/species_clusters/{key}/{value['refgenome']}-ANVIO_VARIABILITY.done")
    #TARGETS.append(f"{working_dir}/species_clusters/{key}/{value['refgenome']}-ANVIO_VARIABILITY-CDN-PNPS.done")
    #TARGETS.append(f"{working_dir}/species_clusters/{key}/{value['refgenome']}-POPULATIONS.done")
    #TARGETS.append(f"{working_dir}/species_clusters/{key}/{value['refgenome']}-NEUTRAL.done")

#DEBUG
#TARGETS = ['/cluster/work/sunagawa/paolil/ocean-popgen/refgenomes/bams/TARA_SAMEA2623868_METAG/TARA_SAMEA2623868_METAG-vs-COMBINED_REF_GENOMES_FOR_POPGEN.done']

# One rule to rule them all ====================================================

print("Starting snakemake")
rule all:
    input: TARGETS

# Include the required rules ===================================================

include: 'popgen-rules-combinedref.smk' 

