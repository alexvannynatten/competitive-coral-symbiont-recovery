#####################################################################
# Last updated 2025-05-23 - Alexander Van Nynatten
# Tidies symportal output and incorperates sample and survey metadata
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Sample metadata
sample_metadata <- read_csv('data/sample_metadata.csv')

# Site metadata
survey_metadata <- read_csv('data/survey_metadata.csv') %>%
  select(site_name, pub.name, lat, lon, level, continous.pressure.2km)

# Subsets header rows for SymPortal output
symportal_header <- read_delim('data/20250527T170756.profiles.absolute.abund_and_meta.txt', 
                               delim = '\t', col_names = FALSE) %>%
  slice(c(2,7)) %>%
  summarise(across(everything(), ~ paste(., collapse = "|"))) %>% 
  unlist()

# Reads the rest of the SymPortal output and attaches the modified header rows
symportal_output <- read_delim('data/20250527T170756.profiles.absolute.abund_and_meta.txt', 
                               delim = '\t', skip = 6) %>%
  slice(1:(n() - 2)) %>%
  setNames(symportal_header) %>%
  select(-1) %>%
  rename(library_ID = 1)

# Categorical definitions for expeditions relative heatwave
heatwave_df <- data.frame(expedition = c('2013', '2014', '2015a', '2015b', '2015c', 
                                         '2016a', '2016b', '2017', '2018', '2019', 
                                         '2023b'),
                          expedition_class = c('1_before', '1_before', '1_before', 
                                               '1_before',  '2_during', '2_during', 
                                               '3_after', '3_after', '3_after', 
                                               '3_after', '4_longafter'))

#####################################################################
# Tidies dataset, removes duplicates, and samples with < 1000 reads
#####################################################################

# Simplified version of the SymPortal ASV table for subsequent analyses
symportal_df <- symportal_output %>%
  pivot_longer(!c(library_ID), names_to = 'Profile', values_to = 'Count') %>% # Converts to tidy dataframe
  separate_wider_delim(Profile, delim = "|", names = c("Clade", "Profile")) %>% # Separates Clade and Type-profile classification
  mutate(Count = as.numeric(Count)) %>%
  filter(Count > 0) %>%
  group_by(library_ID) %>%
  mutate(Count_sum = sum(Count)) %>% # Calculates total reads per sample
  ungroup() %>%
  filter(Count_sum > 1000) %>% # Removes any samples with fewer than 1000 reads
  mutate(Count_freq = Count/Count_sum) %>% # Converts read counts to read frequencies
  left_join(sample_metadata) %>%
  left_join(survey_metadata) %>%
  left_join(heatwave_df)

write.csv(symportal_df, 'data/01_symportal_results_df.csv', 
          row.names = FALSE)

#####################################################################
# Some stats on the corals sampled
#####################################################################

# Number of samples
n_samples <- symportal_df %>%
  distinct(coral_tag, expedition, host_genus) %>%
  group_by(host_genus) %>%
  summarize(n_samples = n())

# Number of colonies
n_colonies <- symportal_df %>%
  distinct(coral_tag, host_genus) %>%
  group_by(host_genus) %>%
  summarize(n_colonies = n())

# Number of tracked colonies
n_repeat_samples <- symportal_df %>%
  distinct(coral_tag, expedition, host_genus) %>%
  group_by(coral_tag, host_genus) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>%
  ungroup() %>%
  group_by(host_genus) %>%
  reframe(n_repeat_samples = n())

# Combines above
summary_table <- n_samples %>%
  left_join(n_colonies) %>%
  left_join(n_repeat_samples) %>%
  bind_rows(summarise(., across(where(is.numeric), sum, na.rm = TRUE))) %>%
  mutate(host_genus = replace(host_genus, n(), "Total"))

# Summary table of the number of samples and colonies for each species
summary_table

#####################################################################
# Supplementary table of all samples collected by expedition and site
#####################################################################

Suppl_table <- symportal_df %>%
  distinct(pub.name, sample_name, host_genus, level, expedition) %>%
  group_by(pub.name, host_genus, level, expedition) %>%
  summarize(n_samples = n()) %>%
  ungroup() %>%
  arrange(expedition) %>%
  pivot_wider(names_from = expedition, values_from = n_samples) %>%
  arrange(level) %>%
  select(!level) %>%
  replace(is.na(.), 0)

write.csv(Suppl_table, 'figures_tables/Table_S3.csv', row.names = FALSE)

#####################################################################

