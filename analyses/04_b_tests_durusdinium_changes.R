#####################################################################
# Last updated 2025-05-25 - Alexander Van Nynatten
# Compares proportion of colonies hosting Durusdinium by expedition
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')

#####################################################################
# Durusdinium results for supplementary figures
#####################################################################

# Identifies which colonies host Durusinium as a binary variable
D_binary_df <- symportal_df %>%
  filter(!expedition_class == '2_during') %>%
  filter(Count_freq > 0.1) %>%
  mutate(Durusdinium_binary = ifelse(Clade == 'D', 1, 0)) %>% # identifies which colonies host Durusdinium
  group_by(sample_name, expedition, expedition_class, site_name, host_genus) %>%
  reframe(Durusdinium_binary = max(Durusdinium_binary)) # For when colonies host more than one symbiont

#####################################################################
# Calculates the number of samples hosting Durusdinium at each site by expedition
#####################################################################

# Montipora samples hosting Durusdinium
D_binary_df %>%
  filter(host_genus == 'Montipora') %>%
  filter(Durusdinium_binary == 1)

# Proportion of Pocillopora colonies hosting Durusdinium at each site surveyed 
prop_D_site <- D_binary_df %>%
  filter(host_genus == 'Pocillopora') %>%
  group_by(host_genus, site_name, expedition_class, expedition) %>%
  mutate(n_samples = n()) %>% # number of colonies hosting Durusdinium
  group_by(host_genus, site_name, expedition, expedition_class, n_samples) %>%
  reframe(n_D = sum(Durusdinium_binary)) %>%
  mutate(prop_D = n_D/n_samples) %>%
  mutate(expedition_plot = recode(expedition_class,
                          "1_before" = "2013-2015b",
                          "3_after" = "2016b-2019",
                          "4_longafter" = "2023b"))

# Violin plots of the results
ggplot() + 
  geom_violin(data = prop_D_site,
              aes(expedition_plot, prop_D),
              fill = 'grey70', scale = 'width', draw_quantiles = c(0.5)) +
  geom_jitter(data = prop_D_site, 
              aes(expedition_plot, prop_D),
              shape = 21, width = 0.1, fill = 'white', size = 2) + 
  theme_test(base_size = 10) + 
  ylab('Proportion of colonies hosting Durusdinium') + 
  theme(axis.title.x = element_blank())
  
ggsave('figures_tables/FigS4_Pocillopora_Durusdinium_incidence.pdf', height = 100, width = 110, units = "mm")

# Subsamples only sites that were surveyed before heatwave and in 2023
before_longafter_sites <- prop_D_site %>%
  distinct(site_name, expedition_class) %>%
  filter(expedition_class %in% c('1_before', '4_longafter')) %>%
  group_by(site_name) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>%
  .$site_name

# Tests for significance between pre- and 2023 Durusdinium incidence
wilcox.test(prop_D~expedition_class, data = prop_D_site %>%
              filter(expedition_class %in% c('1_before', '4_longafter')) %>%
              filter(site_name %in% before_longafter_sites),
            exact = FALSE, paired = FALSE)

#####################################################################

