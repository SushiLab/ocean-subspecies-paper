# ===== Script to produce the full pangenome and gene catalog membership =====

log_msg <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
fmt_n   <- function(x) format(x, big.mark = ",", scientific = FALSE)

# Libraries ======
library(tidyverse)
library(data.table)

# Loading data ======

# Load the genomes summary 
genome_summary = read_csv("https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-summary.csv.gz") %>% dplyr::rename(Species = `dRep Species Cluster`)
species_summary = genome_summary %>% group_by(Species) %>% summarize(n = n(), motus_assignment = any(!is.na(`mOTUs Species Cluster`))) %>% filter(n >= 10 & motus_assignment) %>% select(Species)
genome_summary = genome_summary %>% filter(Species %in% species_summary$Species) %>% select(Genome, Species, `Mean Completeness`)
species_summary = data.table(species_summary)
genome_summary = data.table(genome_summary)

log_msg("Identified %s species with at least 10 genomes & mOTUs assignmnent.",
        fmt_n(nrow(species_summary)))

# Load the genomes scaffold information
genome_scaffold_membership = fread(
  "https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-scaffolds-membership.tsv.gz",
  header = FALSE, col.names = c("Scaffold", "Genome")
)[Genome %in% genome_summary$Genome]

# Load the gene catalog information
gene_catalog_membership = fread("https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/gene-catalog-membership.tsv.gz")

# This file only ships `gene` and `representative` -- no separate Scaffold column.
# The scaffold is embedded in the gene ID itself (SAMPLE-scaffold_N-gene_M), and
# stripping the trailing "-gene_M" reproduces genome_scaffold_membership's own
# Scaffold naming exactly, so we derive it here before joining.
gene_catalog_membership[, Scaffold := sub("-gene_[0-9]+$", "", gene)]

# ---- Join, restricting to species of interest as early as possible ----
# The original code did three left_join()s and only implicitly relied on filtering
# later, so every downstream step -- especially the per-species loop -- scanned the
# FULL gene catalog, including genes from every species NOT in species_summary.
# genome_scaffold_membership is already restricted to genomes of interest, so joining
# on Scaffold with nomatch = NULL (an inner join) drops everything irrelevant right
# away, before the two smaller joins even run.
gene_catalog_membership = genome_scaffold_membership[gene_catalog_membership, on = "Scaffold", nomatch = NULL]
gene_catalog_membership = genome_summary[gene_catalog_membership, on = "Genome", nomatch = NULL]
gene_catalog_membership = species_summary[gene_catalog_membership, on = "Species", nomatch = NULL]

# Compute pangenome =====

gene_catalog_membership[, n_genomes := uniqueN(Genome), by = Species]

# correct_for_incompleteness(), vectorized: for a given (Species, representative) group,
# "incompleteness" is the total incompleteness of every genome of that species that does
# NOT carry this representative gene. Instead of looping per species and re-filtering the
# full table 271 times (the actual bottleneck, compounded by a growing rbind()), compute:
#   incompleteness = total_incompleteness_of_species - incompleteness_of_genomes_present
# both terms obtained with a single grouped pass over the (now much smaller) table.
genome_summary[, incompl := (100 - `Mean Completeness`) / 100]
total_incompl_by_species = genome_summary[, .(total_incompl = sum(incompl)), by = Species]

# One row per (Species, representative, Genome) -- collapses duplicate gene-level rows,
# since only genome *presence* matters for n and incompleteness from here on
presence_tbl = unique(gene_catalog_membership[, .(Species, representative, Genome, n_genomes)])
presence_tbl = genome_summary[, .(Species, Genome, incompl)][presence_tbl, on = c("Species", "Genome")]

summary_table = presence_tbl[, .(
  n               = .N,                # already unique Genomes at this point
  incompl_present = sum(incompl),
  n_genomes       = unique(n_genomes)
), by = .(Species, representative)]

summary_table = total_incompl_by_species[summary_table, on = "Species"]
summary_table[, incompleteness := total_incompl - incompl_present]
summary_table = summary_table[, .(Species, representative, n, incompleteness, n_genomes)]
setDF(summary_table)   # hand back a plain data frame for the tidyverse code below

pangenome_frequency = summary_table %>%
  mutate(gc_frequency = (n + incompleteness)/n_genomes) %>%
  group_by(Species) %>%
  mutate(min_frequency = min(gc_frequency[gc_frequency > 0])) %>%
  ungroup()

pangenome_table = pangenome_frequency %>%
  mutate(pangenome = case_when(
    gc_frequency >= 0.75                                       ~ "core",
    gc_frequency >= 1.15 * min_frequency & gc_frequency < 0.75 ~ "shell",
    gc_frequency < 1.15 * min_frequency                        ~ "cloud"
  )
  )

fwrite(pangenome_table, "PANGENOME-GC-FREQUENCY-FULL.tsv.gz", sep = "\t")

pangenome_summary = pangenome_table %>%
  group_by(Species) %>%
  summarize(
    all_gc = n(),
    all_gc_safe = length(unique(representative)),
    core = sum(pangenome == "core"),
    shell = sum(pangenome == "shell"),
    cloud = sum(pangenome == "cloud")
  )

pangenome_summary %>%
  select(Species, all_gc, core_gc = core, shell_gc = shell, cloud_gc = cloud) %>%
  write_tsv("PANGENOMES-FULL-SUMMARY.tsv")

gene_catalog_membership_pangenome = gene_catalog_membership %>%
  select(Species, gene, representative) %>%
  left_join(pangenome_table %>% select(Species, representative, pangenome))

fwrite(gene_catalog_membership_pangenome, "PANGENOME-GENE-CATEGORY-FULL.tsv.gz", sep = "\t")
