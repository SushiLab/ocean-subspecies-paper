
rm(list = ls())

# We load the genome summary from the OMD =====

genome_summary = read_csv("https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-summary.csv.gz") %>% dplyr::rename(Species = `dRep Species Cluster`)

# Add taxonomy to the reference genome per species =====

species_representatives = googlesheets4::read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_representatives")
species_representatives = species_representatives %>% left_join(genome_summary %>% select(representative_genome = Genome, `GTDB Taxonomy`))
species_representatives %>% googlesheets4::write_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_representatives")

# We then select species with at least 10 genomes and mOTUs-level assignments =====

species_summary = genome_summary %>% group_by(Species) %>% summarize(n = n(), motus_assignment = any(!is.na(`mOTUs Species Cluster`))) %>% filter(n >= 10 & motus_assignment)

# Filter based on species-samples coverage before SNV calling =====

coverage = read_tsv("/Users/paolil/polybox/PhD/Exploratorium/popgen/data/merged-coverage-per-sample-per-species.tsv.gz")

species_for_SNVs = coverage %>%
  filter((vertical_coverage >= 5) & frac_core_detected >= 0.9) %>%
  group_by(Species) %>%
  summarize(n_samples_with_vh_coverage = n()) %>%
  filter(n_samples_with_vh_coverage >= 20)

species_samples_for_SNVs = coverage %>%
  filter((vertical_coverage >= 5) & frac_core_detected >= 0.9) %>%
  filter(Species %in% species_for_SNVs$Species)

species_for_SNVs %>%
  googlesheets4::write_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_for_SNVs")

# Filter based on SNV coverage and phaseability before population pipeline ===== 

phaseability = read_tsv("/Users/paolil/polybox/PhD/Exploratorium/popgen/data/phaseability.all.merged.tsv")

species_for_pop = species_samples_for_SNVs %>%
  left_join(phaseability %>% select(Species = species, Sample = sample_id, Phaseability = value) %>% unique) %>%
  mutate(samples_with_SNV_coverage = !is.na(Phaseability),
         samples_phaseable = samples_with_SNV_coverage & Phaseability >= 0.8) %>%
  left_join(species_for_SNVs) %>%
  group_by(Species) %>%
  summarize(n_samples_with_vh_coverage = unique(n_samples_with_vh_coverage),
            n_samples_with_SNV_coverage = sum(samples_with_SNV_coverage),
            n_samples_phaseable = sum(samples_phaseable)) %>%
  filter(n_samples_phaseable >= 15)

species_samples_for_pop = species_samples_for_SNVs %>%
  left_join(phaseability %>% select(Species = species, Sample = sample_id, Phaseability = value) %>% unique) %>%
  mutate(samples_with_SNV_coverage = !is.na(Phaseability),
         samples_phaseable = samples_with_SNV_coverage & Phaseability >= 0.8) %>%
  filter(samples_with_SNV_coverage) %>%
  filter(Species %in% species_for_pop$Species)
  
species_for_pop %>%
  googlesheets4::write_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_for_pop")

# Filter based on subspecies assignment & confidence for summary =====

populations_summary = read_tsv("/Users/paolil/polybox/PhD/Exploratorium/popgen/data/POPULATIONS-SUMMARY.tsv")

species_with_subspecies = species_samples_for_pop %>%
  left_join(populations_summary) %>%
  filter(umap_trustworthiness >= 0.6 & DBCV >= 0.4) %>%
  left_join(species_for_pop) %>%
  group_by(Species) %>%
  summarize(n_samples_with_vh_coverage = unique(n_samples_with_vh_coverage),
            n_samples_with_SNV_coverage = sum(samples_with_SNV_coverage),
            n_samples_phaseable = sum(samples_phaseable),
            umap_trustworthiness = unique(umap_trustworthiness),
            DBCV = unique(DBCV),
            n_subspecies = length(unique(Cluster[Cluster != -1]))) %>%
  filter(n_subspecies > 1) # Actually nothing is filtered here

species_sample_with_subspecies = species_samples_for_pop %>%
  left_join(populations_summary) %>%
  filter(Species %in% species_with_subspecies$Species)

species_with_subspecies %>%
  googlesheets4::write_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_with_subspecies")

species_sample_with_subspecies %>%
  googlesheets4::write_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "samples_subspecies_assignment")

species_with_subspecies %>%
  left_join(species_representatives) %>%
  mutate(phylum = gsub(";c__.*", "", `GTDB Taxonomy`)) %>%
  pull(phylum) %>%
  table() %>% 
  length
