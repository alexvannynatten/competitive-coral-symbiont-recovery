#####################################################################
# Last updated 2025-05-29 - Alexander Van Nynatten
# Tests the recovery of pre-heatwave symbiont lineages in Pocillopora
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')
profile_groups <- read_csv('data/04_symbiont_lineages.csv')

# Survey dates for comparing recovery over time
Sampling_dates <- data.frame(
  expedition = sort(unique(symportal_df$expedition)),
  Date = as.Date(c('2013-08-01', 
                   '2014-08-01',
                   '2015-01-01',
                   '2015-05-01',
                   '2015-07-01',
                   '2016-03-01',
                   '2016-11-01', 
                   '2017-07-01',
                   '2018-07-01',
                   '2019-07-01',
                   '2023-09-01')
  ))


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
  left_join(Sampling_dates) %>%
  mutate(Profile_type_binary = ifelse(Profile %in% pre_heatwave_profiles$Profile, 1, 0)) %>%
  group_by(sample_name, Date, expedition, expedition_class, coral_tag) %>%
  reframe(Profile_type_binary = max(Profile_type_binary)) %>%
  mutate(days_after = Date - as.Date("2016-04-01")) %>%
  filter(days_after > 0) %>%
  mutate(yrs_after = as.numeric(days_after / 365.25))

# Checks the proportion of samples in 2023 hosting pre-heatwave states
symbiont_recovery_df %>%
  filter(expedition == '2023b') %>%
  group_by(Profile_type_binary) %>%
  summarize(n = n())

# Plots binomial regression of probability of hosting pre-heatwave profiles x year
ggplot(data = symbiont_recovery_df,
       aes(yrs_after, Profile_type_binary)) +
  geom_jitter(shape = 21, width = 0.2, height = 0.04, colour = '#054a8d', size = 3, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), 
              se = TRUE, colour = '#054a8d') +
  xlab("Years after heatwave") +
  ylab("Hosting pre-heatwave symbiont") + 
  scale_x_continuous(breaks = 0:7) + 
  ylim(-0.05,1.05) + 
  theme_test()

ggsave('figures_tables/Fig3a_preheatwave_recovery.pdf', height = 80, width = 80, units = "mm")

  # Statistical analysis using binomial distribution 
Recovery_x_yr <- glm(Profile_type_binary ~ yrs_after, family = binomial(link = "logit"), 
                          data = symbiont_recovery_df)

# Summary of statistical output
summary(Recovery_x_yr)

# For 50% recovery
target_prob <- 0.5
logit_target <- log(target_prob / (1 - target_prob))
(logit_target - coef(Recovery_x_yr)[1]) / coef(Recovery_x_yr)[2]

# For 95% recovery
target_prob <- 0.95
logit_target <- log(target_prob / (1 - target_prob))
(logit_target - coef(Recovery_x_yr)[1]) / coef(Recovery_x_yr)[2]

#####################################################################
# Probability of hosting any Cladocopium each expedition after heatwave
#####################################################################

# New dataframe classifying each sample as hosting Cladocopium (1) or not (0)
cladocopium_recovery_df <- symportal_df %>% 
  filter(Count_freq > 0.1) %>%
  filter(host_genus == 'Pocillopora') %>% # only Pocillopora
  left_join(Sampling_dates) %>%
  mutate(Profile_type_binary = ifelse(Clade == 'C', 1, 0)) %>%
  group_by(sample_name, Date, expedition, coral_tag) %>%
  reframe(Profile_type_binary = max(Profile_type_binary)) %>%
  mutate(days_after = Date - as.Date("2016-04-01")) %>%
  filter(days_after > 0) %>%
  mutate(yrs_after = as.numeric(days_after / 365.25))

# Checks the proportion of samples in 2023 hosting pre-heatwave states
cladocopium_recovery_df %>%
  filter(expedition == '2023b') %>%
  group_by(Profile_type_binary) %>%
  summarize(n = n())

# Plots binomial regression of probability of hosting pre-heatwave profiles x year
ggplot(data = cladocopium_recovery_df,
       aes(yrs_after, Profile_type_binary)) +
  geom_jitter(shape = 21, width = 0.2, height = 0.05, colour = '#4481bd', size = 3, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), 
              se = TRUE, colour = '#4481bd') +
  xlab("Years after heatwave") +
  ylab("Hosting Cladocopium") + 
  scale_x_continuous(breaks = 0:7) + 
  ylim(-0.05,1.05) + 
  theme_test()

ggsave('figures_tables/Fig3b_cladocopium_recovery.pdf', height = 80, width = 80, units = "mm")

# Statistical analysis using binomial distribution 
Recovery_C_yr <- glm(Profile_type_binary ~ yrs_after, family = binomial(link = "logit"), 
                     data = cladocopium_recovery_df)

# Summary of statistical output
summary(Recovery_C_yr)

# For 50% recovery
target_prob <- 0.5
logit_target <- log(target_prob / (1 - target_prob))
(logit_target - coef(Recovery_C_yr)[1]) / coef(Recovery_C_yr)[2]

# For 90% recovery
target_prob <- 0.95
logit_target <- log(target_prob / (1 - target_prob))
(logit_target - coef(Recovery_C_yr)[1]) / coef(Recovery_C_yr)[2]

#####################################################################
# Simulations of Pocillopora samples equal in number of Montipora colonies
#####################################################################

# Number of Montipora colonies sampled collected after heatwave
symportal_df %>%
  filter(!expedition_class %in% c('1_before', '2_during')) %>%
  filter(host_genus == 'Montipora') %>%
  distinct(sample_name) %>%
  nrow()
  
# Loops through random subsamples of 2023 Pocillopora samples on n = Montipora samples
preh_detected <- numeric()
for(i in 1:1000){
  
  # Random subsamples
  sampled_colonies <- symportal_df %>%
    filter(expedition_class %in% c('3_after', '4_longafter')) %>%
    filter(host_genus == 'Pocillopora') %>%
    distinct(coral_tag, .keep_all = TRUE) %>%
    slice_sample(n = 8) # equal to the number of Montipora samples
  
  # Tests how many of the subsampled colonies host pre-heatwave profiles
  n_preheatwave <- symbiont_recovery_df %>%
    filter(expedition_class %in% c('3_after', '4_longafter')) %>%
    filter(sample_name %in% sampled_colonies$sample_name) %>%
    filter(Profile_type_binary == 1) %>%
    nrow()
  
  preh_detected <- c(preh_detected, n_preheatwave) 
  
}

# Median number of samples with pre-heatwave profiles if fewer Pocillopora corals sampled
median(preh_detected) 

# Statistical comparison of rarefied symbiont sampling to Montipora observations (n = 0 pre-heatwave profiles)
t.test(preh_detected, mu = 0)

#####################################################################

