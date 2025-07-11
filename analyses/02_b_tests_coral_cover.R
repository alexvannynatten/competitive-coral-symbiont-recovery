#####################################################################
# Last updated 2025-01-31 - Alexander Van Nynatten
# Tests differences in coral cover across expeditions and plots results
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)
library(biostat)
library(dunn.test)
library(multcompView)

# Coral cover from surveys
coral_cover_df <- read_csv("data/02_coralnet_results_df.csv") %>%
  mutate(expedition = ifelse(expedition == '23b', '2023', expedition))

# Metadata for sites
survey_df <- read_csv("data/survey_metadata.csv")

################################################################################
## Tables number of quadrat images collected by expedition
################################################################################

# Suppl table of number of quadrats by site
supp_table_1 <- coral_cover_df %>%
  filter(!is.na(pub.name)) %>%
  filter(!is.na(expedition_class)) %>%
  distinct(Quadrat, expedition, pub.name, pressure_level) %>%
  group_by(pub.name, expedition, pressure_level) %>%
  summarize(n_quadrats = n()) %>%
  arrange(expedition) %>%
  pivot_wider(names_from = expedition, values_from = n_quadrats) %>%
  arrange(pressure_level) %>%
  select(-pressure_level) %>%
  replace(is.na(.), 0)
  
write.csv(supp_table_1, 'figures_tables/Table_S2.csv', row.names = FALSE)

################################################################################
## Summarizes coral cover at each site by expedition for comparisons of recovery
################################################################################

# Dataframe containing the % cover of each species classified by CoralNet by site
site_x_expd_coral_cov_df <- coral_cover_df  %>%
  group_by(Species, pub.name, expedition, pressure_level, expedition_class) %>%
  summarize(Species_sub_points = n()) %>%
  ungroup() %>%
  group_by(pub.name, expedition, pressure_level, expedition_class) %>%
  mutate(Total_points = sum(Species_sub_points)) %>%
  ungroup() %>%
  mutate(Percent_sp_cover = Species_sub_points / Total_points) %>%
  complete(Species, nesting(expedition, pressure_level, pub.name, expedition_class)) %>%
  replace_na(list(Percent_sp_cover = 0)) %>%
  mutate(Percent_sp_cover = Percent_sp_cover * 100)

# Dataframe for plot 1 summarizing the coral cover by species 
p1 <- site_x_expd_coral_cov_df %>%
  filter(!is.na(pub.name)) %>%
  filter(Species %in% c('Montipora.foliose', 'Pocillopora.eydouxi')) %>%
  filter(expedition_class %in% c('1_before', '2_during', '3_after', '4_longafter')) %>%
  mutate(expedition_class = factor(expedition_class, levels = c('4_longafter', '3_after', '2_during', '1_before')))

################################################################################
## Posthoc analysis of the differences in coral cover over expeditions
################################################################################

# Kruskal Wallis test on coral cover by expedition
p1 %>%
  group_by(Species) %>%
  do(broom::tidy(kruskal.test(x = .$expedition, g = .$Percent_sp_cover)))

# Posthoc tests
letter_df <- data.frame()
for(i in c('Montipora.foliose', 'Pocillopora.eydouxi')){

# Subsamples species
p1_sub <- p1 %>%
    filter(Species == i) 

# For plotting the letters just above the points
p1_sub_max <- p1_sub %>%
  group_by(expedition) %>%
  summarize(max_Percent_sp_cover = max(Percent_sp_cover) + (max(p1_sub$Percent_sp_cover) / 10)) %>%
  mutate(group = expedition) %>%
  select(group, max_Percent_sp_cover)

# Conducts the test
dunn_results <- dunn.test(p1_sub$Percent_sp_cover, p1_sub$expedition, method = "bh")

p_values <- dunn_results$P.adjusted

groups <- unique(p1_sub$expedition)

# Initialize an empty matrix for storing p-values
comparison_matrix <- matrix(NA, nrow = length(groups), ncol = length(groups))
rownames(comparison_matrix) <- colnames(comparison_matrix) <- groups

# Fill the matrix with pairwise p-values
for (j in seq_along(dunn_results$comparisons)) {
  # Extract group pair from comparison label
  pair <- unlist(strsplit(dunn_results$comparisons[j], " - "))
  
  # Assign the adjusted p-value to the appropriate matrix cells
  comparison_matrix[pair[1], pair[2]] <- dunn_results$P.adjusted[j]
  comparison_matrix[pair[2], pair[1]] <- dunn_results$P.adjusted[j]
}

# Convert p-values to significance letters
group_letters <- multcompLetters(comparison_matrix)$Letters

# Makes dataframe from the letters
letter_df_sub <- data.frame(group_letters) %>%
  mutate(Species = i) %>%
  mutate(group = row.names(.)) %>%
  left_join(p1_sub_max)

letter_df <- rbind(letter_df, letter_df_sub)

}

################################################################################
## Plot 1 Coral cover comparisons by year
################################################################################

ggplot() + 
  geom_boxplot(data = p1,
               aes(x = expedition, y = Percent_sp_cover), 
               fill = 'grey80') +
  geom_text(data = letter_df,
               aes(x = group, y = max_Percent_sp_cover, label = group_letters)) +
  geom_vline(xintercept = c(4.5, 6.5), colour = 'firebrick', linetype = 2) + 
  facet_grid(Species ~ ., scale = 'free_y') + 
  theme_test(base_size = 10) + 
  ylab('Substrate cover (%)') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('figures_tables/Fig2a_Coral_cover.pdf', height = 100, width = 110, units = "mm")

################################################################################
## Plot 2 Recovery in 2023 from post-heatwave lows
################################################################################

# Dataframe for figure 2 comparing 2023 coral cover to low point after heatwave
p2 <- p1 %>%
  filter(expedition_class %in% c('4_longafter', '3_after')) %>%
  group_by(Species, pub.name, expedition_class) %>%
  summarize(min_Percent_sp_cover = min(Percent_sp_cover)) %>%
  ungroup() %>%
  group_by(Species, pub.name) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  pivot_wider(names_from = expedition_class, values_from = min_Percent_sp_cover) %>%
  mutate(change_in_prop_cover = (`4_longafter` - `3_after`)) %>%
  filter(change_in_prop_cover > 0) %>%
  mutate(Species = ifelse(Species == 'Montipora.foliose',
                                     'Montipora aequituberculata',
                                     'Pocillopora grandis'))

# Test of significance
t.test(log10(change_in_prop_cover) ~ Species, data = p2)

# Plots the change in coral cover at each site as a boxplot
ggplot() + 
  geom_boxplot(data = p2,
               aes(x = change_in_prop_cover, y = Species), 
                   fill = 'grey80') +
  theme_test(base_size = 10) + 
  xlab('Change in coral cover as a % of total substrate') + 
  theme(axis.title.y = element_blank())

ggsave('figures_tables/FigS5_Cover_change.pdf', height = 100, width = 150, units = "mm")

# How much the percent of substrate covered changes from post-heatwave lows to 2023 by species
p2 %>%
  group_by(Species) %>%
  reframe(Med_change_in_prop_cover = median(change_in_prop_cover))

################################################################################
