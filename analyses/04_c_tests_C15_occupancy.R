#####################################################################
# Last updated 2025-05-17 - Alexander Van Nynatten
# Detection probability of C15 based on pre-heatwave sampling
#####################################################################

library(tidyverse)

symportal_df <- read_csv('data/01_symportal_results_df.csv')

# Total Montipora samples collected before the heatwave
n_corals <- symportal_df %>%
  filter(host_genus == 'Montipora') %>%
  filter(expedition_class %in% c('1_before')) %>%
  distinct(sample_name) %>%
  reframe(n = n()) %>%
  .$n

n <- n_corals # Number of Montipora samples collected before the heatwave
p <- 0.1 # detection probability based on ITS2 sequencing cutoffs

pc = 1 - (1 - p)^n #Equation for calculating cumulative detection probability

pc

#####################################################################







