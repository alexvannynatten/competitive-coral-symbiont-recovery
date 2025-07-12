#####################################################################
# Last updated 2025-05-25 - Alexander Van Nynatten
# Tests the recovery of pre-heatwave lineages in tracked colonies
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')
profile_groups <- read_csv('data/04_symbiont_lineages.csv')

#####################################################################
# Identifies when colonies were first sampled for tagged colonies
#####################################################################

# Expedition each colony was first sampled
First_sampled <- symportal_df %>%
  mutate(First_sampled = expedition) %>%
  arrange(First_sampled) %>%
  group_by(coral_tag) %>%
  select(First_sampled, coral_tag) %>%
  slice(1) %>%
  ungroup()

#####################################################################
# Makes a database to show symbiont states in re-sampled colonies
#####################################################################

# Counts of each symbiont profile each expedition (by first sampled)
coral_symportal_barplot <- symportal_df %>%
  ungroup() %>%
  left_join(First_sampled) %>%
  left_join(profile_groups) %>%
  filter(!Profile_cols %in% c('tan')) %>% # removes fragmented profiles
  mutate(Profile_cols = ifelse(is.na(Profile_cols), 'grey', Profile_cols)) %>%
  mutate(Profile_cols = ifelse(Clade == 'D', "#EC3F3A", Profile_cols)) %>%
  group_by(expedition, host_genus, Profile_cols, First_sampled) %>%
  summarize(rel_profile_ab = sum(Count_freq)) %>%
  ungroup() %>%
  complete(expedition, host_genus, First_sampled)

# Labels (total counts by expedition) for the bars in each plot
coral_symportal_labels <- symportal_df %>%
  left_join(First_sampled) %>%
  distinct(expedition, host_genus, First_sampled, coral_tag) %>%
  group_by(expedition, host_genus, First_sampled) %>%
  summarize(n_sum = n())

# Plots sequence results by expedition and expedition each was first sampled
ggplot() + 
  geom_bar(data = coral_symportal_barplot, 
           aes(x = rel_profile_ab, y = First_sampled, fill = Profile_cols),
           stat = 'identity', colour = NA, position = 'fill') +
  geom_text(data = coral_symportal_labels,
            aes(x = 0.5, y = First_sampled, label = n_sum),
            colour = 'grey90', size = 2) + 
  facet_grid(host_genus ~ expedition) + 
  theme_test() + 
  scale_y_discrete(limits=rev) + 
  scale_x_continuous(n.breaks = 3) +
  theme(text = element_text(size = 6)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) %>%
  scale_fill_identity() + 
  xlab('Total RRA') + 
  ylab('First sampled')

ggsave('figures_tables/Fig3c_tracked_colonies.pdf', height = 130, width = 180, units = "mm")


#####################################################################

symportal_df %>%
  left_join(profile_groups) %>%
  left_join(First_sampled) %>%
  filter(First_sampled == '2023b') %>%
  group_by(Profile_cols) %>%
  summarize(n = n())

