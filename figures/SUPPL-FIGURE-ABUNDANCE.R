rm(list = ls())

# Libraries ------------------------------------------------------------------------------

library(tidyverse)

# Load data ------------------------------------------------------------------------------

motus_profile = read_tsv("https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/motus-profiles.tsv.gz", skip = 2)

subspecies = googlesheets4::read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_with_subspecies")

motus_map = googlesheets4::read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_for_coverage")

samples = googlesheets4::read_sheet(ss = "1Npt0J8dbFCfUi3Lugw1KQIYTv4G3uRih3Gdx-pibxIw", sheet = "Samples")

# Process data ---------------------------------------------------------------------------

# Identify the relevant motus
popgen_motus = subspecies %>% 
  left_join(motus_map) %>% View()
  separate_rows(motu_match, sep = ";") %>%
  filter(!is.na(motu_match)) %>%
  pull(`motu_match`) %>% unique()

# Get the profiles
motus_profile_processed = motus_profile %>%
  mutate(motu = gsub(".*\\[|\\]", "", `#consensus_taxonomy`)) %>%
  select(-`#consensus_taxonomy`) %>%
  gather(key = Sample, value = abd, -motu) %>%
  filter(Sample %in% samples$Sample) %>%
  group_by(Sample) %>%
  mutate(relabd = abd / sum(abd)) %>%
  ungroup()

motus_profile_popgen = motus_profile_processed %>%
  filter(motu %in% popgen_motus) %>%
  group_by(Sample) %>%
  summarize(relabd = sum(relabd)) %>%
  left_join(samples)

summary(motus_profile_popgen$relabd)

motus_profile_popgen %>%
  ggplot() +
  geom_jitter(aes(y = "dummy", x = relabd), height = .35, width = 0, size = .3, color = "lightgrey") +
  geom_boxplot(aes(y = "dummy", x = relabd), outlier.shape = NA, alpha = .6, fill = NA) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 6) +
  xlab("Cumulative relative abundance") +
  ylab("Size fraction") +
  theme_minimal() +
  theme(text = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(2,2,2,2, 'mm'),
        legend.position = 'none')

ggsave("~/polybox/PhD/Exploratorium/popgen/publication/suppl_figures/popgen-ecol-relevance.pdf", width = 175, height = 30, units = 'mm')
ggsave("~/polybox/PhD/Exploratorium/popgen/publication/suppl_figures/popgen-ecol-relevance.png", width = 175, height = 30, units = 'mm')
