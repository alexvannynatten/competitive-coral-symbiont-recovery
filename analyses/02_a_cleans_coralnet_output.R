################################################################################
# Last updated 2025-05-23 - Alexander Van Nynatten
# Decodes coralNet labels, cleans survey data, and incorporates survey metadata
################################################################################

################################################################################
## Loads required libraries and datasets
################################################################################

library(tidyverse)
coral.raw <- read_csv("data/CoralNet_raw_5July2024_2013-2023b.csv")
coral.codes <- read_csv("data/Coralnet_labels.csv")
survey_data <- read_csv('data/survey_metadata.csv') %>%
  select(site_name, pub.name, pressure.group, level) %>%
  mutate(pressure_level = paste(level, pressure.group)) %>%
  mutate(site_name = as.character(site_name))

# Some higher level classifications of each expedition
heatwave_df <- data.frame(expedition = c('2013', '2014', '2015a', '2015b', '2015c', 
                                         '2016a', '2016b', '2017', '2018', '2019', 
                                         '23b'),
                          expedition_class = c('1_before', '1_before', '1_before', 
                                               '1_before',  '2_during', '2_during', 
                                               '3_after', '3_after', '3_after', 
                                               '3_after', '4_longafter'))

################################################################################
## Some transformations of the dataset
################################################################################

# Renames columns and removes data from photoquads that are not part of this study
coral_cov_df <- coral.raw %>%
  select(Name, Label) %>%
  separate(Name, sep="_", into=c("Field.Season", "Site","Quadrat")) %>%
  mutate(Quadrat = gsub("\\.jpg$|\\.jpeg$", '', Quadrat)) %>%
  filter(!grepl("DEEP", Site)) %>%
  filter(!grepl("MPQ", Quadrat)) %>%
  mutate(site_name = gsub('site', '', Site))
  
# Adds taxonomic and life history information for points and counts number of points per quadrat
# Removes shadows and unclassified points
# Removes quadrats with fewer than 90 points
coral_cov_df_2 <- coral_cov_df %>%
  left_join(coral.codes) %>%
  filter(!is.na(Species)) %>%
  group_by(Field.Season, site_name, Quadrat) %>%
  mutate(Total_points = n()) %>%
  ungroup() %>%
  filter(Total_points > 90)

# Number of quadrats initially
coral_cov_df %>%
  mutate(uniq_id = paste(Field.Season, site_name, Quadrat)) %>%
  distinct(uniq_id) %>%
  nrow(.)

# Appends some important site and survey metadata
coral_cov_df_3 <- coral_cov_df_2 %>%
  mutate(expedition = gsub('KI', '', Field.Season)) %>%
  left_join(heatwave_df) %>%
  left_join(survey_data) %>%
  select(-c(Total_points, Field.Season, Site, Label, coralnet.tag))
  
################################################################################
## Exports the table as a .csv file
################################################################################

write.csv(coral_cov_df_3, 'data/02_coralnet_results_df.csv', row.names = FALSE)

################################################################################
