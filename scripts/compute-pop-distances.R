# This script produces the summary distance table based on individually exported distances 

library(tidyverse)
library(googlesheets4)

gs4_deauth()

populations = read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "samples_subspecies_assignment")

rep_genomes = read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_representatives")

prefix = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/OMDV1_POPGEN_EAN/scratch/processed/export_distances/"

out = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/OMDV1_POPGEN_EAN/data/processed/summaries/POPULATIONS-DIST-SUMMARY.tsv.gz"

# pre-define variables
dist_types = cols(.default = col_double(), .delim = "\t")
res = vector("list", length(unique(populations$Species)))
names(res) = unique(populations$Species)

dist_files = tibble(species = unique(populations$Species)) %>%
    left_join(rep_genomes %>% select(species = Species, refg = representative_genome), by = "species") %>%
    mutate(dist_file = paste0(prefix, refg, "-ANVIO_VARIABILITY-NT-DIST.tsv"),
           exists    = file.exists(dist_file))

assertthat::assert_that(!any(is.na(dist_files$refg)))

missing_dist_files = dist_files %>% filter(!exists)
if (nrow(missing_dist_files) > 0){
    message(nrow(missing_dist_files), " of ", nrow(dist_files), " distance files missing:\n",
         paste0("  - ", missing_dist_files$species, " -> ", missing_dist_files$dist_file, collapse = "\n"))
    stop()
}
# If it fails: run export-distances-fix.sh in data/processed/species_clusters

message("All ", nrow(dist_files), " distance files found.")

for (i in unique(populations$Species)){
    refg = rep_genomes %>% filter(Species == i) %>% pull(representative_genome)
    dist_file = paste0(prefix, refg, "-ANVIO_VARIABILITY-NT-DIST.tsv")
    message(paste0("Processing species ", i, " and distance file ", dist_file, "..."))

    tmp_table = populations %>% filter(Species == i, Cluster != -1)
    tmp_samples = tmp_table$Sample

    tmp_dist = read_tsv(dist_file,
                        col_select = c(samples = 1, all_of(tmp_samples)),
                        col_types  = cols(`...1` = col_character(),
                                          .default = col_double()),
                        name_repair = "unique_quiet",
                        progress   = FALSE) %>%
               filter(samples %in% tmp_samples) %>%
               select(-samples)

    tmp_dist[upper.tri(tmp_dist, diag = T)] = NA
    tmp_dist_proc = tmp_dist %>%
        mutate(species = i,
               refg = refg) %>%
        mutate(sample_1 = names(tmp_dist)) %>%
        gather(key = sample_2, value = dist, -species, -refg, -sample_1) %>%
        filter(!is.na(dist))

    res[[i]] = tibble(tmp_dist_proc)
}

species_dist = bind_rows(res)

species_dist_proc = species_dist %>%
    left_join(populations %>% select(species = Species, refg = ref_genome, sample_1 = Sample, cluster_1 = Cluster)) %>%
    left_join(populations  %>% select(species = Species, refg = ref_genome, sample_2 = Sample, cluster_2 = Cluster)) %>%
    mutate(comp = ifelse(cluster_1 == cluster_2, "Same population", "Different populations"))

message(paste0("Writing combined output to ", out, "..."))

species_dist_proc %>% write_tsv(out)
