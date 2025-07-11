#####################################################################
# Last updated 2025-05-25 - Alexander Van Nynatten
# Tests the recovery of pre-heatwave lineages across disturbance gradient 
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')
profile_groups <- read_csv('data/04_symbiont_lineages.csv')

#####################################################################
# Classifies pre-heatwave profiles
#####################################################################

# List of pre-heatwave symbiont profiles in Pocillopora
pre_heatwave_profiles <- profile_groups %>%
  filter(Profile_cols == '#054a8d')

#####################################################################
# Calculates the probability of a colony hosting pre-heatwave symbiont by yr
#####################################################################

# New dataframe classifying each sample as hosting pre-heatwave profile (1) or not (0)
symbiont_recovery_df <- symportal_df %>%
  filter(Count_freq > 0.1) %>%
  filter(host_genus == 'Pocillopora') %>% # only Pocillopora
  mutate(sqrt.continous.pressure.2km = sqrt(continous.pressure.2km)) %>%
  mutate(Profile_type_binary = ifelse(Profile %in% pre_heatwave_profiles$Profile, 1, 0)) %>%
  group_by(sample_name, expedition, host_genus, coral_tag, sqrt.continous.pressure.2km) %>%
  reframe(Profile_type_binary = max(Profile_type_binary)) %>%
  filter(expedition %in% c('2023b'))
  
#####################################################################
# Plots the logistic binomial regression
#####################################################################

# Plots binomial regression of probability of hosting pre-heatwave profiles by disturbance
ggplot(data = symbiont_recovery_df,
       aes(sqrt.continous.pressure.2km, Profile_type_binary, colour = sqrt.continous.pressure.2km)) +
  geom_jitter(shape = 21, width = 1.5, height = 0.05, size = 3) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, 
              colour = 'grey9', linetype = 2) +
  ylab("Hosting pre-heatwave symbiont") + 
  xlab("Human disturbance gradient") + 
  theme_test(base_size = 10) + 
  xlim(-5,90) +
  scale_colour_gradient2(low = "#0096C7", mid = '#6E9075', high = "#825F48", midpoint = 40, na.value = NA)
  
ggsave('figures_tables/Fig4a_disturbed_recovery.pdf', height = 85, width = 140, units = "mm")

# Statistical analysis using binomial distribution 
Disturbance_output <- glm(Profile_type_binary ~ sqrt.continous.pressure.2km, family = binomial(link = "logit"),
                          data = symbiont_recovery_df)

# Summary of statistical output
summary(Disturbance_output)

#####################################################################

