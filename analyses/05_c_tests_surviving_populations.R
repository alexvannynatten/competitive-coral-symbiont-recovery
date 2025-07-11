#####################################################################
# Last updated 2025-05-25 - Alexander Van Nynatten
# Compares symbiont assemblages in surviving colonies with others
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)
library(vegan)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')
profile_groups <- read_csv('data/04_symbiont_lineages.csv')

#####################################################################

symportal_df <- symportal_df %>%
  left_join(profile_groups) %>%
  filter(!Profile_cols %in% c('tan')) %>% # removes fragmented profiles
  mutate(Profile_cols = ifelse(is.na(Profile_cols), 'grey', Profile_cols)) %>%
  mutate(Profile_cols = ifelse(Clade == 'D', "#EC3F3A", Profile_cols))

# For plotting by date instead of expedition
sampling_df <- data.frame(
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

# Loads coral size data
coral_size <- read_csv('data/coral_size.csv')

################################################################################
# Makes dataframe and plots coral size by expedition
################################################################################

# Generates a dataframe of the coral size based on when each colony was first sampled
coral_size_plot <- coral_size %>%
  group_by(coral_tag, expedition) %>%
  slice_head(n = 1) %>%
  left_join(sampling_df) %>%
  mutate(colony_width1 = as.numeric(max_colony_width)) %>%
  ungroup()

################################################################################
# Test if community composition is differnet in 2023 samples at high disturbance sites
################################################################################

# List of corals that are larger than expected based on published growth rates
Larger_than_expected <- coral_size_plot %>%
  filter(expedition == '2023b') %>%
  filter(!is.na(colony_width1)) %>%
  filter(colony_width1 > 37.5) # 5 * 7.5 April 2016 to October 2023

# Samples sampled earlier in the heatwave 
Sampled_early_heatwave <- coral_size_plot %>%
  filter(!expedition %in% c('2023b'))

# Combines above
Surviving_colonies <- unique(c(Larger_than_expected$coral_tag, Sampled_early_heatwave$coral_tag))

################################################################################
# Test if community composition is differnt in 2023 samples at high disturbance sites
################################################################################

# Line showing the growth rate estimate
fast_growth_line <- data.frame(Date = as.Date(c('2016-04-01', '2017-04-01', '2018-04-01', '2019-04-01', '2023-04-01', '2024-04-01')),
                               estimated_colony_width_cm = c(0, 5, 10, 15, 35, 40))

# Dataframe of the size each colony was when first sampled 
coral_size_plot <- coral_size_plot %>%
  mutate(Colony_grouping = ifelse(coral_tag %in% Surviving_colonies, 'Survivor', 'Recovered population')) %>%
  arrange(expedition) %>%
  group_by(coral_tag) %>%
  slice_head(n = 1) %>%
  mutate(Colony_grouping = ifelse(Date < as.Date('2016-01-01'), 'Pre-heatwave', Colony_grouping))
  
# Plots the size when first sampled by date
ggplot() + 
  geom_line(data = fast_growth_line, 
            aes(x = Date, y = estimated_colony_width_cm),
            colour = 'grey4', linetype = 'dashed') + 
  geom_point(data = coral_size_plot, 
              aes(x = Date, y = colony_width1, colour = Colony_grouping),
              shape = 21) + 
  geom_vline(xintercept = as.Date(c('2015-05-28', '2016-04-13')),
             colour = 'firebrick', linetype = 2) +
  theme_test(base_size = 12) + 
  scale_colour_manual(values = c('grey', 'dodgerblue', '#ee6363')) +
  ylab('Colony width (cm)')

ggsave('figures_tables/FigS6a_coral_size.pdf', height = 100, width = 150, units = "mm")

################################################################################
# Makes matrix of ITS2 profile grouping level abundance data per colony
################################################################################

# Subsamples dataframe to include only samples surveyed multiple times before and after
abund_df <- symportal_df %>%
  filter(host_genus == 'Pocillopora') %>%
  group_by(sample_name, Profile_cols, coral_tag, expedition_class, level) %>%
  reframe(Count_freq = sum(Count_freq)) %>%
  filter(expedition_class %in% c('4_longafter')) %>%
  mutate(Grouping = ifelse(coral_tag %in% Surviving_colonies, 'Survior', 'General')) %>%
  mutate(row_names = paste(Grouping, coral_tag, sep = '_')) %>%
  select(row_names, Profile_cols, Count_freq) %>%
  pivot_wider(names_from = Profile_cols, values_from = Count_freq) %>%
  replace(is.na(.), 0)

# Converts dataframe to matrix
abund_matrix <- as.matrix(abund_df[ ,c(2:ncol(abund_df))])
rownames(abund_matrix) <- abund_df$row_names

#####################################################################
# Calculates distances between samples and tests differences with PCORA
#####################################################################

# Runs Jaccard dissimilarity comparison
bray_mat <- vegdist(abund_matrix, method = "bray")

# Matrix in correct format for adonis
bray_mat_matrix <- as.matrix(bray_mat)

# Dataframe of grouping variables
test_df <- data.frame(row_names = row.names(bray_mat_matrix)) %>%
  separate_wider_delim(row_names, names = c('Before_2023', 'Sample'), delim = '_')

# Runs test
sub_test <- vegan::adonis2(bray_mat_matrix ~ Before_2023, 
                           data=test_df, permutations=9999)

print(sub_test)

# Formats for PCOA visualization
pcoa_result <- cmdscale(bray_mat, eig = TRUE)
pcoa_result$eig[1:2]

pcoa_df <- as.data.frame(scores(pcoa_result)) %>%
  mutate(Sample_code_1 = rownames(.)) %>%
  separate_wider_delim(Sample_code_1, "_", names = c("Grouping", "coral_tag"))

# Calculates the centroid for each grouping variable
centroids <- pcoa_df %>%
  group_by(Grouping) %>%
  summarize(Dim1 = mean(Dim1), Dim2 = mean(Dim2))

ggplot() + 
  geom_jitter(data = pcoa_df,
              aes(Dim1, Dim2, colour = Grouping), 
              width = 0.02, height = 0.02, shape = 21, size = 3) + 
  geom_point(data = centroids,
             aes(Dim1, Dim2, colour = Grouping), 
             shape = 3, size = 3) + 
  theme_test(base_size = 12) +
  stat_ellipse() +  # Draws an ellipse around each group
  scale_colour_manual(values = c('dodgerblue', 'indianred2')) +
  xlab(paste0("PC1 (", round(pcoa_result$eig[1],1),")")) +
  ylab(paste0("PC2 (", round(pcoa_result$eig[2],1),")"))
  
ggsave('figures_tables/FigS6b_survivor_bias.pdf', height = 100, width = 125, units = "mm")

#####################################################################
# Extra plot separating survivors from other colonies for samples first sampled in 2023
#####################################################################

# Colonies first sampled in 2023
First_sampled_2023 <- symportal_df %>%
  mutate(First_sampled = expedition) %>%
  arrange(First_sampled) %>%
  group_by(coral_tag) %>%
  select(First_sampled, coral_tag) %>%
  slice(1) %>%
  ungroup() %>%
  filter(First_sampled == '2023b')

# Classifies colonies first sampled in 2023 as survivors 
survivor_plot_2023 <- symportal_df %>%
  filter(host_genus == 'Pocillopora') %>%
  filter(coral_tag %in% First_sampled_2023$coral_tag) %>%
  filter(expedition_class %in% c('4_longafter')) %>%
  mutate(Before_2023 = ifelse(coral_tag %in% Larger_than_expected$coral_tag, 'Survivor', 'General population')) %>%
  group_by(Before_2023,Profile_cols, Profile) %>%
  reframe(Total_reads = sum(Count_freq)) %>%
  group_by(Before_2023) %>%
  mutate(RRA = Total_reads / sum(Total_reads))

ggplot() +
  geom_bar(data = survivor_plot_2023, 
           aes(x = RRA, y = Before_2023, fill = Profile_cols),
           stat = 'identity', colour = NA, position = 'fill') + 
  scale_fill_identity() + 
  theme_test()
  
ggsave('figures_tables/Fig3c_tracked_colonies_extra.pdf')

#####################################################################
# Some other useful stats
#####################################################################

# Number of colonies exclusively hosting D
all_D <- symportal_df %>%
  mutate(Before_2023 = ifelse(coral_tag %in% Surviving_colonies, 'Survior', 'General')) %>%
  filter(expedition_class %in% c('4_longafter')) %>%
  filter(Clade == 'D') %>%
  filter(Count_freq == 1) %>%
  select(Before_2023, sample_name, Profile, host_genus, expedition, level)

# Survivors
survivor_check <- symportal_df %>%
  filter(host_genus == 'Pocillopora') %>%
  filter(expedition_class %in% c('4_longafter')) %>%
  mutate(Before_2023 = ifelse(coral_tag %in% Larger_than_expected$coral_tag, 'Yes', 'No')) %>%
  filter(Before_2023 == 'Yes')

#####################################################################

