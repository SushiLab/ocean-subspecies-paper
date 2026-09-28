# This script plots the outer layer of the figure with the phylogenetic distribution of 
# subspecies

rm(list = ls())

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(treeio) # YuLab-SMU/treeio (devtools::install_github("YuLab-SMU/treeio"))
library(ggtree) # YuLab-SMU/ggtree (devtools::install_github("YuLab-SMU/ggtree"))
library(tidytree)

# Prepare data ---------------------------------------------------------------------------

artree = "https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/gtdb-phylogeny/gtdbtk.ar122.classify.tree"
bctree = "https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/gtdb-phylogeny/gtdbtk.bac120.classify.tree"
genomes_summary = read_csv("https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-summary.csv.gz")

subspecies = googlesheets4::read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "species_with_subspecies")

# Archaea ================================================================================

# Inner plot -----------------------------------------------------------------------------

artrnw <- drop.tip(read.newick(artree),
                   genomes_summary %>%
                     dplyr::filter(!`dRep Representative Genome` |  `GTDB Taxonomy` == "N/A") %>% pull(Genome)
)

artrnw_metadata = 
  as_tibble(artrnw) %>% 
  left_join(subspecies %>% left_join(genomes_summary %>% filter(`dRep Representative Genome`) %>% select(Species = `dRep Species Cluster`, label = Genome)) %>% mutate(user_genome = TRUE)) %>%
  mutate(user_genome = ifelse(is.na(user_genome), FALSE, user_genome))

artree_plot = ggtree(artrnw, layout="fan", open.angle=285, branch.length = "none", size = 0.35, alpha = 1, ) %<+%
  #artree_plot = ggtree(artrnw, layout="fan", open.angle=270, size = 0.35, alpha = 1) %<+%
  artrnw_metadata + # add metadata to the tree data
  scale_color_identity() +
  coord_polar(theta = 'y', start = 1.5*pi, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        line = element_line(size = 0.2))
artree_plot

# Prepare outer barplots -----------------------------------------------------------------

backbone_ar = artree_plot$data %>% filter(isTip)
circle_ratio_ar = 75/360
total_elt_ar = nrow(backbone_ar) * (1/circle_ratio_ar)

max(artree_plot$data$y, na.rm = T)/250

# Pop Barplot ----------------------------------------------------------------------------

table_figure_pop_melt_ar = artree_plot$data %>%
  filter(isTip) %>%
  mutate(y_r = DescTools::RoundTo(y, 40)) %>%
  group_by(y_r) %>%
  summarize(x = 0,
            y = mean(y),
            subsp = any(!is.na(n_subspecies)))

outer_layer_pop_ar = table_figure_pop_melt_ar  %>%
  ggplot() +
  geom_col(aes(x = y, y = subsp), fill = "#287C8E", width = 21) +
  geom_vline(xintercept = 0, size = 0.12) +
  geom_vline(xintercept = max(table_figure_pop_melt_ar$y, na.rm=T) +  1, size = 0.12) +
  theme_minimal() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 0, 0, 0, 'mm'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.12),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_polar(start = 1.5*pi, direction = -1) +
  xlim(0, total_elt_ar) +
  scale_y_continuous(breaks = c(0, 2, 5, 10, 20, 30), limits = c(-110,35))
outer_layer_pop_ar

ggsave("/Users/paolil/polybox/PhD/Exploratorium/popgen/publication/figures/Figure-tree-material/figure-tree-archaeal_outer_pop.pdf", outer_layer_pop_ar, width = 140, height = 140, units = 'mm')

# Bacteria ================================================================================

# Inner plot -----------------------------------------------------------------------------

bctrnw <- drop.tip(read.newick(bctree),
                   genomes_summary %>%
                     dplyr::filter(!`dRep Representative Genome` |  `GTDB Taxonomy` == "N/A") %>% pull(Genome)
)

bctrnw_metadata = 
  as_tibble(bctrnw) %>% 
  left_join(subspecies %>% left_join(genomes_summary %>% filter(`dRep Representative Genome`) %>% select(Species = `dRep Species Cluster`, label = Genome)) %>% mutate(user_genome = TRUE)) %>%
  mutate(user_genome = ifelse(is.na(user_genome), FALSE, user_genome))

bctree_plot = ggtree(bctrnw, layout="fan", open.angle=285, branch.length = "none", size = 0.35, alpha = 1, ) %<+%
  bctrnw_metadata + # add metadata to the tree data
  scale_color_identity() +
  coord_polar(theta = 'y', start = 1.5*pi, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        line = element_line(size = 0.2))
bctree_plot

# Prepare outer barplots -----------------------------------------------------------------

backbone_bc = bctree_plot$data %>% filter(isTip)
circle_ratio_bc = 3/4
total_elt_bc = nrow(backbone_bc) * (1/circle_ratio_bc)

# Pop Barplot ----------------------------------------------------------------------------

table_figure_pop_melt_bc = bctree_plot$data %>%
  filter(isTip) %>%
  mutate(y_r = DescTools::RoundTo(y, 180)) %>%
  group_by(y_r) %>%
  summarize(x = 0,
            y = mean(y),
            subsp = any(!is.na(n_subspecies)))

outer_layer_pop_bc = table_figure_pop_melt_bc  %>%
  ggplot() +
  geom_col(aes(x = y, y = subsp), fill = "#287C8E", width = 90) +
    geom_vline(xintercept = 0, size = 0.12) +
  geom_vline(xintercept = max(table_figure_pop_melt_bc$y, na.rm=T) +  1, size = 0.12) +
  theme_minimal() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 0, 0, 0, 'mm'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.12),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt_bc) +
  scale_y_continuous(breaks = c(0, 2, 5, 10, 20, 30), limits = c(-110,35))
outer_layer_pop_bc

ggsave("/Users/paolil/polybox/PhD/Exploratorium/popgen/publication/figures/Figure-tree-material/figure-tree-bacterial_outer_pop.pdf", outer_layer_pop_bc, width = 140, height = 140, units = 'mm')

