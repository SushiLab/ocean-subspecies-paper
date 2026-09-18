# Format metadata

library(tidyverse)

# Load samples map =====

samples_map = read_tsv("~/polybox/PhD/Exploratorium/popgen/metadata/samples-map-final-1669.tsv") %>% rename(barcode = `Internal Sample Name`) %>% mutate(barcode = gsub("_METAG", "", barcode))

# Load metadata =====

popgen_metadata_raw = read_tsv("~/polybox/PhD/Exploratorium/popgen/data/popgen-metadata.tsv") # In situ measurements

metadata_models = read_tsv("~/polybox/PhD/Exploratorium/popgen/metadata/MatchedParameters_station_metadata_for_dominic.tsv") %>% select(-`...1`) # model metadata

# Connectivity =====

# We need the connectivity eigen values:

connectivity_raw = read_tsv("/Users/paolil/polybox/PhD/Exploratorium/Connectivity/Connectivity-processed.tsv", guess_max = 10000) %>%
  gather(key = Sample2, value = connectivity, -Sample) %>%
  dplyr::rename(Sample1 = Sample)

missing_connect = connectivity_raw %>%
  group_by(Sample1) %>%
  summarize(na = sum(is.na(connectivity))) %>%
  filter(na > 1000) %>%
  pull(Sample1)

connectivity_tbl = connectivity_raw %>%
  filter(!Sample1 %in% missing_connect & !Sample2 %in% missing_connect)

connectivity_tbl %>%
  ggplot() +
  geom_density(aes(x = connectivity))

connectivity_tbl = connectivity_tbl %>%
  mutate(connectivity = ifelse(is.na(connectivity) | connectivity > 4000, 2*4000, connectivity))

connectivity_dist = connectivity_tbl %>%
  spread(Sample2, connectivity) %>%
  select(-Sample1) %>%
  as.dist()

connectivity_pcoa = ape::pcoa(connectivity_dist)
eigen = connectivity_pcoa$values$Eigenvalues
eigen[eigen < 0] = 0
eigen = round(eigen/sum(eigen)*100)
eigen[eigen > 0]

connectivity_vectors = as_tibble(connectivity_pcoa$vectors) %>%
  add_column(Sample = labels(connectivity_dist), .before = 1)

# Combined metadata =====

samples_metadata_for_popgen = samples_map %>%
  left_join(popgen_metadata_raw) %>%
  mutate(depth_num = as.numeric(depth)) %>%
  select(Sample, barcode, ocean_province, station, latitude, longitude, size_fraction, depth_num, depth_layer) %>%
  left_join(connectivity_vectors %>% select(Sample, connectivity_axis_1 = `Axis.1`, connectivity_axis_2 = `Axis.2`)) %>%
  left_join(metadata_models %>% select(station, depth_num,
                                       temperature_model = `temperature_degreeCeslius`,
                                       oxygen_model = `oxygen_micromol_kg-1`,
                                       iron_model = `iron_nanomol_L-1`,
                                       phosphate_model = `phosphate_micromol_kg-1`,
                                       nitrate_model = `nitrate_micromol_kg-1`,
                                       salinity_model = salinity)) %>%
  left_join(popgen_metadata_raw %>% select(Sample,
                                           temperature_situ = temperature,
                                           oxygen_situ = oxygen,
                                           iron_situ = iron,
                                           phosphate_situ = phosphate,
                                           nitrate_situ = nitrate,
                                           salinity_situ = salinity 
  ))

samples_metadata_for_popgen %>% googlesheets4::sheet_write(ss = "1Npt0J8dbFCfUi3Lugw1KQIYTv4G3uRih3Gdx-pibxIw", sheet = "Samples")
