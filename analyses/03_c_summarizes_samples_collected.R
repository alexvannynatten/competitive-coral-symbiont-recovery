#####################################################################
# Last updated 2025-05-28- Alexander Van Nynatten
# Summarizes the samples collected and surveys conducted by expedition
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

samples_df <- read_csv("data/sample_metadata.csv")
symportal_df <- read_csv("data/01_symportal_results_df.csv")
coral_cover_df <- read_csv("figures_tables/Table_S2.csv")

  
heatwave_df <- data.frame(expedition = c('2013', '2014', '2015a', '2015b', 
                                         '2015c', '2016a', 
                                         '2016b', '2017', '2018', '2019', 
                                         '2023'),
                          expedition_class = c('1_before', '1_before', '1_before', '1_before',  
                                               '2_during', '2_during', 
                                               '3_after', '3_after', '3_after', '3_after', 
                                               '4_longafter'),
                          expedition_date = c('17 July 2013 - 7 August 2013', 
                                              '20 August 2014 - 10 September 2014', 
                                              '21 January 2015 - 4 Febrary 2015', 
                                              '28 April 2015 - 12 May 2015',
                                              '1 July 2015 - 28 July 2015', 
                                              '16 March 2016 - 6 April 2016', 
                                              '9 November 2016 - 23 November 2016', 
                                              '5 July 2017 - 2 August 2017', 
                                              '13 June 2018 - 4 July 2018', 
                                              '17 July 2019 - 14 August 2019', 
                                              '6 September 2023 - 4 october 2023'))


#####################################################################
# Sums of photoquadrats and samples collected for each species by expedition 
#####################################################################

sum_quads <- coral_cover_df %>%
  select(-pub.name) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "expedition", values_to = "total photoquadrats")

sum_dna <- samples_df %>%
  filter(sample_name %in% symportal_df$sample_name) %>%
  mutate(expedition = ifelse(expedition == '2023b', '2023', expedition)) %>%
  group_by(host_genus, expedition) %>%
  reframe(n_samples = n()) %>%
  mutate(host_genus = ifelse(host_genus == 'Pocillopora', 'Pocillopora grandis DNA samples', 'Montipora aequituberculata DNA samples')) %>%
  pivot_wider(names_from = host_genus, values_from = n_samples)

output_data <- heatwave_df %>%
  left_join(sum_quads) %>%
  left_join(sum_dna)
  
write.csv(output_data, 'figures_tables/Table_S1.csv', row.names = FALSE)

#####################################################################