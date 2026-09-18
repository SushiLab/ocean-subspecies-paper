# Ocean Populations ======================================================================

remove(list = ls())

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)

# Functions ------------------------------------------------------------------------------

run_subsp_multinomial <- function(tbl, sp, min_samples = 7, max_it = 10**4){
  require(nnet)
  require(caret)
  
  assertthat::assert_that(all(c("Species", "Cluster", "connectivity_axis_1", "connectivity_axis_2", "temperature", "oxygen", "iron", "phosphate", "nitrate", "salinity") %in% names(tbl)))
  
  # Remove missing values and rescale variables
  p_reg = na.omit(tbl) %>%
    filter(Species == sp) %>%
    select(-Species) %>%
    mutate(connectivity_axis_1 = scale(connectivity_axis_1)[,1],
           connectivity_axis_2 = scale(connectivity_axis_2)[,1],
           temperature = scale(temperature)[,1],
           oxygen = scale(oxygen)[,1],
           iron = scale(iron)[,1],
           phosphate = scale(phosphate)[,1],
           nitrate = scale(nitrate)[,1],
           salinity = scale(salinity)[,1])
  # mutate(connectivity_axis_1 = scales::rescale(connectivity_axis_1, to = c(0, 1)),
  #    connectivity_axis_2 = scales::rescale(connectivity_axis_2, to = c(0, 1)),
  #    temperature = scales::rescale(temperature, to = c(0, 1)),
  #    oxygen = scales::rescale(oxygen, to = c(0, 1)),
  #    iron = scales::rescale(iron, to = c(0, 1)),
  #    phosphate = scales::rescale(phosphate, to = c(0, 1)),
  #    nitrate = scales::rescale(nitrate, to = c(0, 1)),
  #    salinity = scales::rescale(salinity, to = c(0, 1)))
  
  # Filter populations with enough data
  p_reg_clusters = p_reg %>%
    group_by(Cluster) %>%
    summarize(n = n()) %>%
    filter(n >= min_samples)
  
  if (nrow(p_reg_clusters) < 2) {
    
    # Skip if not >1 cluster (=Cluster) with enough data
    print(paste0("Not enough data for specices ", sp, "."))
    return(NULL)
    
  } else {
    
    p_reg = p_reg %>% filter(Cluster %in% p_reg_clusters$Cluster) %>% mutate(Cluster = as.factor(Cluster))
    
    # Setting the reference
    p_reg$Cluster <- relevel(p_reg$Cluster, ref = min(levels(p_reg$Cluster)))
    
    # Accounting for co-linearity
    cor_matrix = WGCNA::cor(p_reg %>% select(-Cluster))
    cor_matrix[lower.tri(cor_matrix, diag = T)] = NA
    
    cor_tbl = cor_matrix %>%
      as_tibble() %>%
      mutate(env1 = row.names(cor_matrix)) %>%
      gather(key = env2, value = cor, -env1) %>%
      filter(!is.na(cor) & abs(cor) > 0.75)
    
    if (nrow(cor_tbl) > 0){message(paste0(c("Found co-linearity, removing:", cor_tbl$env2), collapse = " "))}
    
    p_reg = p_reg %>% select(-cor_tbl$env2)
    
    # Training the full multinomial model
    multinom_full <- multinom(Cluster ~ ., data = p_reg, na.action = na.omit, maxit = max_it)
    
    # Checking the model
    multinom_full$convergence
    summary(multinom_full)
    coef(multinom_full)
    exp(coef(multinom_full))
    
    # Training subsets
    submodel <- function(data, variables, var_name = NA, full_model = multinom_full){
      if (any(variables %in% names(data))){
        variables = variables[variables %in% names(data)] # account for correlation between connectivity_axis_1 and connectivity_axis_2
        multinom_submodel <- multinom(Cluster ~ ., data = data %>% select(-all_of(variables)), na.action = na.omit, maxit = max_it)
        anova_submodel = anova(multinom_submodel, full_model, test = "Chisq")
        tmp = tibble(variable = ifelse(is.na(var_name), variables, var_name),
                     LR = anova_submodel$`LR stat.`[2],
                     pval = anova_submodel$`Pr(Chi)`[2])
      } else {
        tmp = tibble(variable = ifelse(is.na(var_name), variables, var_name),
                     LR = NA,
                     pval = NA)
      }
      return(tmp)
    }
    
    anova_no_connec = submodel(p_reg, c("connectivity_axis_1", "connectivity_axis_2"), var_name = "connectivity")
    anova_no_temp = submodel(p_reg, c("temperature"))
    anova_no_oxy = submodel(p_reg, c("oxygen"))
    anova_no_iron = submodel(p_reg, c("iron"))
    anova_no_phos = submodel(p_reg, c("phosphate"))
    anova_no_nitrate = submodel(p_reg, c("nitrate"))
    anova_no_salinity = submodel(p_reg, c("salinity"))
    
    res = rbind(
      anova_no_connec,
      anova_no_temp,
      anova_no_oxy,
      anova_no_iron,
      anova_no_phos,
      anova_no_nitrate,
      anova_no_salinity
    ) %>%
      add_column(Species = sp, .before = 1)
    
    return(res)
  }
}

run_depth_multinomial <- function(tbl, sp, min_samples = 7){
  require(nnet)
  require(caret)
  
  assertthat::assert_that(all(c("Species", "Cluster", "depth_num") %in% names(tbl)))
  
  # Remove missing values and rescale variables
  p_reg = na.omit(tbl) %>%
    filter(Species == sp) %>%
    select(-Species) %>%
    mutate(depth_num = scale(depth_num)[,1])
  
  # Filter populations with enough data
  p_reg_clusters = p_reg %>%
    group_by(Cluster) %>%
    summarize(n = n()) %>%
    filter(n >= min_samples)
  
  if (nrow(p_reg_clusters) < 2) {
    # Skip if not >1 cluster (=Cluster) with enough data
    print(paste0("Not enough data for specices ", sp, "."))
    return(NULL)
    
  } else {
    p_reg = p_reg %>% filter(Cluster %in% p_reg_clusters$Cluster) %>% mutate(Cluster = as.factor(Cluster))
    
    # Setting the reference
    p_reg$Cluster <- relevel(p_reg$Cluster, ref = min(levels(p_reg$Cluster)))
    
    # Training the full multinomial model
    multinom_null = multinom(Cluster ~ 1, data = p_reg, na.action = na.omit, maxit = 10000)
    multinom_depth = multinom(Cluster ~ ., data = p_reg, na.action = na.omit, maxit = 10000)
    
    anova_depth = anova(multinom_null, multinom_depth, test = "Chisq")
    tmp = tibble(variable = "depth_num",
                 LR = anova_depth$`LR stat.`[2],
                 pval = anova_depth$`Pr(Chi)`[2]) %>%
      add_column(Species = sp, .before = 1)
    
    return(tmp)
  }
}

run_sizefrac_multinomial <- function(tbl, sp, min_samples = 7){
  require(nnet)
  require(caret)
  
  assertthat::assert_that(all(c("Species", "Cluster", "dataset", "ocean_province", "depth_layer", "size_fraction") %in% names(tbl)))
  
  # Remove missing values and rescale variables
  p_reg = na.omit(tbl) %>%
    filter(Species == sp) %>%
    select(-Species)
  
  # Filter populations with enough data
  p_reg_clusters = p_reg %>%
    group_by(Cluster) %>%
    summarize(n = n()) %>%
    filter(n >= min_samples)
  
  if (nrow(p_reg_clusters) < 2 | length(unique(p_reg$size_fraction)) <= 1) {
    # Skip if not >1 cluster (=Cluster) with enough data
    print(paste0("Not enough data for specices ", sp, "."))
    return(NULL)
    
  } else {
    p_reg = p_reg %>% filter(Cluster %in% p_reg_clusters$Cluster) %>% mutate(Cluster = as.factor(Cluster))
    
    # Setting the reference
    p_reg$Cluster <- relevel(p_reg$Cluster, ref = min(levels(p_reg$Cluster)))
    
    # removing variables if no variation
    if (length(unique(p_reg$dataset)) == 1){p_reg = p_reg %>% select(-dataset)}
    if (length(unique(p_reg$ocean_province)) == 1){p_reg = p_reg %>% select(-ocean_province)}
    if (length(unique(p_reg$depth_layer)) == 1){p_reg = p_reg %>% select(-depth_layer)}
    
    # Training the full multinomial model and null model
    if (any(c("dataset", "depth_layer", "ocean_province") %in% names(p_reg))){
      multinom_null = multinom(Cluster ~ ., data = p_reg %>% select(-size_fraction), na.action = na.omit, maxit = 10000)
      multinom_sizefrac = multinom(Cluster ~ ., data = p_reg, na.action = na.omit, maxit = 10000)
    } else {
      multinom_null = multinom(Cluster ~ 1, data = p_reg, na.action = na.omit, maxit = 10000)
      multinom_sizefrac = multinom(Cluster ~ ., data = p_reg, na.action = na.omit, maxit = 10000)
    }
    
    # compute significance
    anova_sizefrac = anova(multinom_null, multinom_sizefrac, test = "Chisq")
    tmp = tibble(variable = "sizefrac",
                 LR = anova_sizefrac$`LR stat.`[2],
                 pval = anova_sizefrac$`Pr(Chi)`[2]) %>%
      add_column(Species = sp, .before = 1)
    
    return(tmp)
  }
}

# Load & prep data -----------------------------------------------------------------------

subspecies = googlesheets4::read_sheet(ss = "1fXzHlT3NDnBKruxyM_VWqkMzMeXD3PCs2RAN5dV61pc", sheet = "samples_subspecies_assignment")
metadata = googlesheets4::read_sheet(ss = "1Npt0J8dbFCfUi3Lugw1KQIYTv4G3uRih3Gdx-pibxIw", sheet = "Samples")

subspecies_summary_metadata = subspecies %>%
  left_join(metadata) %>%
  filter(Cluster != -1) %>% # Remove noise points, i.e. not belonging to any cluster 
  mutate(Cluster = ifelse(grepl("Cluster", Cluster), Cluster, paste("Cluster", Cluster)))

# Defined max depth for analysis =========================================================

subspecies_summary_metadata %>%
  select(Sample, depth_num) %>%
  mutate(depth_num = ifelse(depth_num > 200, 200, depth_num)) %>%
  ggplot() +
  geom_histogram(aes(x = depth_num), binwidth = 10)

depth_test = NULL
for (d in 1:200){
  subspecies_summary_metadata_epi_tmp = subspecies_summary_metadata %>%
    filter(depth_layer == "EPI") %>%
    filter(depth_num <= d)
  depth_test = rbind(
    depth_test,
    tibble(depth = d,
           corr = cor.test(subspecies_summary_metadata_epi_tmp$temperature_model, subspecies_summary_metadata_epi_tmp$temperature_situ)$estimate)
  )
}

depth_test %>%
  ggplot() +
  geom_line(aes(x = depth, y = corr)) +
  theme_minimal()

subspecies_summary_metadata_epi = subspecies_summary_metadata %>%
  filter(depth_layer == "EPI") %>%
  filter(depth_num <= 100)

subspecies_summary_metadata_epi %>%
  select(Sample, temperature_situ, oxygen_situ, iron_situ, phosphate_situ, nitrate_situ, salinity_situ, temperature_model, oxygen_model, iron_model, phosphate_model, nitrate_model, salinity_model) %>%
  unique %>%
  gather(key = variable_situ, value = estimate_situ, -all_of(contains("_model")), -Sample) %>%
  gather(key = variable_model, value = estimate_model, - Sample, -all_of(contains("_situ"))) %>%
  mutate(variable_situ = gsub("_situ", "", variable_situ),
         variable_model = gsub("_model", "", variable_model)) %>%
  filter(variable_situ == variable_model) %>%
  ggplot() +
  geom_point(aes(x = estimate_situ, y = estimate_model), size = 0.5) +
  facet_wrap(~variable_situ, scales = "free") +
  xlab("In situ measurements") +
  ylab("Modelled variables") +
  theme_minimal() +
  theme(text = element_text(size = 6))

ggsave("~/polybox/PhD/Exploratorium/popgen/publication/suppl_figures/multinomial-variables-model-vs-situ.jpg", width = 130, height = 120, units= 'mm')

# Multinomial regressions on surface samples based on models =============================

subspecies_summary_metadata_epi_pred_model = subspecies_summary_metadata_epi %>%
  select(Species,
         Cluster,
         temperature = temperature_model,
         oxygen = oxygen_model,
         iron = iron_model,
         phosphate = phosphate_model,
         nitrate = nitrate_model,
         salinity = salinity_model,
         connectivity_axis_1,
         connectivity_axis_2)

res_model = NULL
for (sp in subspecies_summary_metadata_epi_pred_model %>% pull(Species) %>% unique){
  res_model = rbind(res_model, run_subsp_multinomial(subspecies_summary_metadata_epi_pred_model, sp))
}

res_model = res_model %>% mutate(padj = p.adjust(pval, method = 'fdr'))

res_model %>% googlesheets4::sheet_write(ss = "1-FkH8B4muzFsoCP7MxAf04w8E9w6TyjqMu8RHMPaIkU", sheet = "Env model")

length(unique(res_model$Species))

res_model %>% filter(padj <= 0.05) %>% pull(Species) %>% unique() %>% length()
res_model %>% filter(variable != "connectivity") %>% filter(padj <= 0.05) %>% pull(Species) %>% unique() %>% length()
res_model %>% filter(padj <= 0.05) %>% pull(variable) %>% table()

p1 = res_model %>%
  group_by(variable) %>%
  summarize(
    significant = sum(padj <= 0.05, na.rm = T),
    non_significant = sum(padj > 0.05, na.rm = T),
    not_tested = sum(is.na(padj)),
    n = n()
  ) %>% #View()
  gather(key = type, value = number, -n, -variable) %>%
  mutate(type = factor(type, levels = rev(c("significant", "non_significant", "not_tested"))),
         variable = factor(variable, levels = rev(c("temperature", "connectivity", "phosphate", "salinity", "iron", "nitrate", "oxygen")))) %>%
  ggplot() +
  geom_col(aes(x = number, y = variable, fill = type)) +
  scale_fill_manual(values = c("significant" = "#0F3C5F", "non_significant" = "grey60", "not_tested" = "grey80")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

p2 = res_model %>%
  mutate(variable = factor(variable, levels = rev(c("temperature", "connectivity", "phosphate", "salinity", "iron", "nitrate", "oxygen")))) %>%
  filter(padj <= 0.05) %>%
  ggplot() +
  geom_boxplot(aes(y = variable, x = log2(LR)), outliers = F) +
  geom_jitter(aes(y = variable, x = log2(LR)), size = .4) +
  xlim(2, 8) +
  scale_y_discrete(drop = FALSE) +
  xlab("log2(likelihood ratios)") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 6))

p1 | p2

ggsave("~/polybox/PhD/Exploratorium/popgen/publication/figures/Figure-environment-material/multinomial-variables-model.pdf", width = 120, height = 70, units= 'mm')

# Multinomial regressions on surface samples based on in situ data =======================

subspecies_summary_metadata_epi_pred_situ = subspecies_summary_metadata_epi %>%
  filter(depth_layer == "EPI") %>%
  select(Species,
         Cluster,
         temperature = temperature_situ,
         oxygen = oxygen_situ,
         iron = iron_situ,
         phosphate = phosphate_situ,
         nitrate = nitrate_situ,
         salinity = salinity_situ,
         connectivity_axis_1,
         connectivity_axis_2)


res_situ = NULL
for (sp in subspecies_summary_metadata_epi_pred_situ %>% pull(Species) %>% unique){
  res_situ = rbind(res_situ, run_subsp_multinomial(subspecies_summary_metadata_epi_pred_situ, sp))
}

res_situ = res_situ %>% mutate(padj = p.adjust(pval, method = 'fdr'))

res_situ %>% googlesheets4::sheet_write(ss = "1-FkH8B4muzFsoCP7MxAf04w8E9w6TyjqMu8RHMPaIkU", sheet = "Env situ")

length(unique(res_situ$Species))
res_situ %>% filter(padj <= 0.05) %>% pull(variable) %>% table()

# Multinomial regressions on depth =======================================================

res_depth = NULL
for (sp in subspecies_summary_metadata %>% pull(Species) %>% unique){
  res_depth = rbind(res_depth, run_depth_multinomial(subspecies_summary_metadata %>% select(Species, Cluster, depth_num), sp))
}

res_depth = res_depth %>% mutate(padj = p.adjust(pval, method = 'fdr'))

res_depth %>% googlesheets4::sheet_write(ss = "1-FkH8B4muzFsoCP7MxAf04w8E9w6TyjqMu8RHMPaIkU", sheet = "Depth")

length(unique(res_depth$Species))
res_depth %>% filter(padj <= 0.05) %>% pull(Species) %>% unique() %>% length()

p1_depth = res_depth %>%
  group_by(variable) %>%
  summarize(
    significant = sum(padj <= 0.05, na.rm = T),
    non_significant = sum(padj > 0.05, na.rm = T),
    not_tested = sum(is.na(padj)),
    n = n()
  ) %>%
  gather(key = type, value = number, -n, -variable) %>%
  mutate(type = factor(type, levels = rev(c("significant", "non_significant", "not_tested")))) %>%
  ggplot() +
  geom_col(aes(x = number, y = variable, fill = type)) +
  scale_fill_manual(values = c("significant" = "#0F3C5F", "non_significant" = "grey60", "not_tested" = "grey80")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

p2_depth = res_depth %>%
  filter(padj <= 0.05) %>%
  ggplot() +
  geom_boxplot(aes(y = variable, x = log2(LR)), outliers = F) +
  geom_jitter(aes(y = variable, x = log2(LR)), size = .4) +
  xlim(2, 8) +
  scale_y_discrete(drop = FALSE) +
  xlab("log2(likelihood ratios)") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 6))

p1_depth | p2_depth

ggsave("~/polybox/PhD/Exploratorium/popgen/publication/figures/Figure-environment-material/multinomial-variables-model-depth.pdf", width = 120, height = 30, units= 'mm')

# Combine model and situ data ============================================================

res_model = googlesheets4::read_sheet(ss = "1-FkH8B4muzFsoCP7MxAf04w8E9w6TyjqMu8RHMPaIkU", sheet = "Env model")
res_situ = googlesheets4::read_sheet(ss = "1-FkH8B4muzFsoCP7MxAf04w8E9w6TyjqMu8RHMPaIkU", sheet = "Env situ")
res_depth = googlesheets4::read_sheet(ss = "1-FkH8B4muzFsoCP7MxAf04w8E9w6TyjqMu8RHMPaIkU", sheet = "Depth")

rbind(res_model %>% mutate(type = "model"),
      res_situ %>% mutate(type = "situ"),
      res_depth %>% mutate(type = "depth")) %>%
  mutate(type = ifelse(variable == "connectivity", "connectivity", type)) %>%
  group_by(Species, type) %>%
  summarize(association = any(padj <= 0.05, na.rm = T)) %>%
  ungroup() %>%
  spread(key = type, value = association) %>%
  filter(!is.na(model) | !is.na(situ)) %>%
  mutate(association = case_when(
    model ~ "model",
    situ ~ "situ",
    connectivity ~ "connectivity",
    depth ~ "depth",
    .default = "none"
  )) %>%
  group_by(association) %>%
  summarize(n = n()) %>%
  mutate(association = factor(association, levels = rev(c("model", "situ", "depth", "connectivity", "none")))) %>%
  ggplot() +
  geom_col(aes(x = n, y = "dummy", fill = association)) +
  scale_x_continuous(limits = c(0, 150)) +
  scale_fill_manual(values = c("model" = "#076769", "situ" = "#3EA8A6", "connectivity" = "#5D94CB", "depth" = "#10194D", "none" = "#BFBFBF")) +
  theme_minimal()

ggsave("~/polybox/PhD/Exploratorium/popgen/publication/figures/Figure-environment-material/multinomial-variables-model-combined.pdf", width = 120, height = 30, units= 'mm')

# How many species do we have a env model for?
rbind(res_model %>% mutate(type = "model"), res_situ %>% mutate(type = "situ")) %>%
  filter(!is.na(LR)) %>%
  pull(Species) %>% unique %>% length

rbind(res_model %>% mutate(type = "model"), res_situ %>% mutate(type = "situ")) %>%
  filter(!is.na(LR) & padj <= 0.05) %>%
  pull(Species) %>% unique %>% length

stats = rbind(res_model %>% mutate(type = "model"), res_situ %>% mutate(type = "situ")) %>%
  filter(!is.na(LR) & variable != "connectivity") %>%
  group_by(Species, type) %>%
  summarize(sig = any(padj <= 0.05)) %>%
  spread(key = type, value = sig) %>%
  mutate(sig = any(model, situ, na.rm = T)) %>%
  ungroup()

sum(stats$sig)

pull(res_model, Species) %>% unique %>% length

res_model %>%
  group_by(variable) %>%
  summarize(n_tot = n(),
            n_test = sum(!is.na(pval)),
            n_sig = sum(padj <= 0.05, na.rm = T),
            r_tot = n_sig/n_tot,
            r_test = n_sig/n_test)

pull(res_depth, Species) %>% unique %>% length

res_depth %>%
  group_by(variable) %>%
  summarize(n_tot = sum(!is.na(pval)),
            n_sig = sum(padj <= 0.05, na.rm = T),
            ratio = n_sig/n_tot)
