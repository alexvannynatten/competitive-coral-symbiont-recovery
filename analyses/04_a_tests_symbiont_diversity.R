#####################################################################
# Last updated 2025-05-25 - Alexander Van Nynatten
# Visualizes UniFrac distances and plots profile incidence by expedition
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)
library(ape)
library(ggtree)
library(gridExtra)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')

# Inter-profile distance matrix for Cladocopium taxa
C_dist_matrix <- as.matrix(read.table("data/20250527T170756_unifrac_profile_distances_C_sqrt.dist", row.names = 1, sep = '\t'))
C_dist_matrix <- C_dist_matrix[, -1]

#####################################################################
# Generates dendrograms from UniFrac data
#####################################################################

# Converts Cladocopium unifrac dataset to a dendrogram
C_dist_object <- as.dist(C_dist_matrix)
C_cluster_dend <- hclust(C_dist_object, method = "average")
C_profile_dend <- as.phylo(C_cluster_dend) # converts to phylo format for plotting

# Clusters the profiles within 0.0025 dissimilarity of eachother
hc <- hclust(C_dist_object, method = "average")
threshold <- 0.0035
groups <- cutree(hc, h = threshold)

# Assign arbitrary names to groups
group_names <- paste0("Group_", as.numeric(factor(groups)))

# Dataframe of all profiles and groups based on unifrac clusters
group_df <- data.frame(
  Profile = labels(C_dist_object),
  Group = group_names
)

# Assigns each group a colour based on unifrac clusters
group_cols <- group_df %>%
  group_by(Group) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Group) %>%
  mutate(Profile_cols = c('tan','#054a8d','#ae61ac','tan',
                          '#9dbffe', 'tan', '#00b778','tan','tan')) %>%
  select(Group, Profile_cols)

# Adds colours to grouping dataframe
group_df <- group_df %>%
  left_join(group_cols)

# Exports database for subsequent analyses of recovery
write.csv(group_df, 'data/04_symbiont_lineages.csv', row.names = FALSE)

#####################################################################
# Summarizes and formats the data for plotting as a dendrogram
#####################################################################

# Calculates incidence of each symportal profile by expedition for line plot
profile_count_df <- symportal_df %>%
  filter(Count_freq > 0.1) %>% # Removes profiles with RRA less than 0.1
  group_by(expedition, host_genus) %>%
  mutate(n_samples = n_distinct(sample_name)) %>%
  ungroup() %>%
  group_by(expedition, Profile, Clade, host_genus, n_samples) %>%
  reframe(n_detections = n()) %>%
  mutate(Profile_rel_abund = n_detections/n_samples)

# The total times each profile is observed before or after the heatwave
Profile_dend_plot_df <- profile_count_df %>%
  group_by(Profile, host_genus) %>%
  reframe(total_detections_x_sp = sum(n_detections))

Profile_dend_plot_df <- Profile_dend_plot_df %>%
  left_join(group_df)

#####################################################################
# Generates dendrogram with shapes beside tree showing profile variation
#####################################################################

# Fragmented profiles to drop
fragmented_profiles <- Profile_dend_plot_df %>%
  filter(Profile_cols == 'tan')

# Lists Cladocopium profiles in unifrac dendrogram not in count table - always less than 10%
tips_to_remove_C <- data.frame(Profile = C_profile_dend$tip.label) %>%
  left_join(Profile_dend_plot_df) %>%
  mutate(host_genus = ifelse(Profile %in% fragmented_profiles$Profile, NA, host_genus)) %>% # removes fragmented profiles
  filter(is.na(host_genus)) %>%
  .$Profile

# Removes tips in the Cladocopium unifrac dendrogram if not in heatmap
C_profile_dend.pruned <- drop.tip(C_profile_dend, tips_to_remove_C) 

# Adds count data to the dendrogram dataframe
C_abund_df <- data.frame(Profile = C_profile_dend.pruned$tip.label) %>%
  left_join(Profile_dend_plot_df) %>%
  left_join(group_df) %>%
  arrange(host_genus, Profile)
          
# Plots the dendrogram of Cladocopium profiles (branches are unifrac distances)
pC <- ggtree(C_profile_dend.pruned) +
  geom_tiplab(size=2) + 
  geom_treescale(width = 0.001, linesize = 0.5, fontsize = 2)

# Plots dendrogram with circles at tips representing the 
pC + geom_facet(panel = "Abundance", data = C_abund_df, geom = geom_point, 
                 mapping=aes(x = 1, colour = Profile_cols, size = total_detections_x_sp), 
                 shape = 16) + 
  scale_colour_identity() +
  scale_size_continuous(range = c(1,7))

ggsave('figures_tables/Fig2b_symbiont_dend.pdf')

#####################################################################
# Line graph of symbiont recovery
#####################################################################

# Makes dataframe of all profiles detected in each expedition
profile_count_complete_df <- symportal_df %>%
  filter(!Profile %in% fragmented_profiles$Profile) %>% # Removes fragmented profiles
  filter(Count_freq > 0.1) %>% # removes profiles with less than 10% relative read frequency per sample (background lineages)
  group_by(expedition, host_genus) %>%
  mutate(n_samples = n_distinct(sample_name)) %>%
  ungroup() %>%
  left_join(group_df) %>%
  mutate(Profile_cols = ifelse(is.na(Profile_cols), 'grey', Profile_cols)) %>%
  mutate(Profile_cols = ifelse(Clade == 'D', "#EC3F3A", Profile_cols)) %>%
  distinct(sample_name, expedition, Profile_cols, Clade, host_genus, n_samples) %>%
  group_by(expedition, Profile_cols, Clade, host_genus, n_samples) %>%
  reframe(n_detections = n()) %>%
  mutate(Profile_rel_abund = n_detections/n_samples) # calculates the relative abundance of each ITS2 type profile by expedition

# Some modifications of the dataframe to include colours for Durusdinium
profile_count_df_plot <- profile_count_complete_df %>%
  complete(expedition, nesting(Profile_cols, Clade, host_genus)) %>%
  replace_na(list(Profile_rel_abund = 0)) %>%
  mutate(zeros = ifelse(Profile_rel_abund == 0, 'yes', 'no')) %>% # Makes any non-detections an open shape to differentiate rare and totally absent profiles
  ungroup()

# Plots the recovery of symbionts as a line plot
ggplot() + 
  geom_text(data = profile_count_df_plot %>%
              filter(!is.na(n_samples)) %>%
              distinct(host_genus, n_samples, expedition), 
            aes(x = expedition, y = 120, label = n_samples), 
            size = 2, colour = 'grey33') +
  geom_line(data = profile_count_df_plot, 
            aes(x = expedition, Profile_rel_abund * 100, group = Profile_cols, colour = Profile_cols),
            linewidth = 1) +
  geom_point(data = profile_count_df_plot,
             aes(x = expedition, Profile_rel_abund * 100, colour = Profile_cols, shape = zeros), 
             size = 4) +
  scale_shape_manual(values = c(16, 32)) + 
  facet_grid(host_genus ~ .) + 
  theme_test() + 
  scale_colour_identity() + 
  geom_vline(xintercept = c(4.5, 6.5), colour = 'firebrick', linetype = 2) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = "none") + 
  ylab('ITS2 profile incidence (%)')

ggsave('figures_tables/Fig2a_symbiont_recovery.pdf', height = 100, width = 110, units = "mm")

#####################################################################
# Supplementary plots of relative read abundance by sample and expedition
#####################################################################

# Updates the colours to differentiate fragmented lineages
group_cols_supp <- group_cols %>%
  mutate(Profile_cols = c('orange','#054a8d','#ae61ac','pink',
                          '#9dbffe', 'pink4', '#00b778','tan','yellow'))

profile_groups_suppl_df <- group_df %>%
  select(!Profile_cols) %>%
  left_join(group_cols_supp) %>%
  select(Profile, Profile_cols)


# Loop to generate the supplementary figures showing Profile RRA by coral colony / expedition
j <- 2 # for supp figure file name
for(i in c('Pocillopora', 'Montipora')){
  
    # Dataframe for plotting symbiont proportions per coral colony
    supp_fig_plot <- symportal_df %>%
      filter(host_genus == i) %>%
      arrange(expedition) %>%
      mutate(expedition = factor(expedition, levels = unique(expedition))) %>%
      arrange(level, pub.name, expedition) %>%
      mutate(level_tag = paste(pub.name , coral_tag)) %>%
      mutate(level_tag = factor(level_tag, levels = rev(unique(level_tag)))) %>%
      left_join(profile_groups_suppl_df) %>%
      mutate(Profile_cols = ifelse(Clade == 'D', "#EC3F3A", Profile_cols)) %>%
      mutate(Profile_cols = ifelse(is.na(Profile_cols), 'grey50', Profile_cols))
    
    # Splits into two dataframes for plotting on one page
    first_half <- supp_fig_plot %>%
      distinct(level_tag) %>%
      slice(1:(n() / 2)) %>%
      .$level_tag
      
    supp_fig_plot_1 <- supp_fig_plot %>%
      filter(level_tag %in% first_half)
      
    supp_fig_plot_2 <- supp_fig_plot%>%
      filter(!level_tag %in% first_half)
      
    p1 <- ggplot() + 
      geom_col(data = supp_fig_plot_1,
               aes(x = Count_freq, level_tag, 
                   fill = Profile_cols),
               colour = 'black') + 
      facet_grid(. ~ expedition, scales = 'free_x', drop = FALSE) + 
      theme_test() + 
      scale_fill_identity() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.key.size = unit(0.2, "cm")
      ) + 
      theme(
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5)
      ) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    p2 <- ggplot() + 
      geom_col(data = supp_fig_plot_2,
               aes(x = Count_freq, level_tag, 
                   fill = Profile_cols),
               colour = 'black') + 
      facet_grid(. ~ expedition, scales = 'free_x', drop = FALSE) + 
      theme_test() + 
      scale_fill_identity() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        legend.key.size = unit(0.2, "cm")
      ) + 
      theme(
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5)
      ) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    p_combined <- grid.arrange(p1, p2, ncol = 2)
    
    file_name <- paste0("FigS", j, "_", i, "_symbiont_RRA.pdf")
    ggsave(paste0('figures_tables/', file_name), 
           plot = p_combined, height = 10, width = 7, limitsize = FALSE)
    
    j <- j+1
    
}

#####################################################################
# Plots legend for supplementary figures
#####################################################################

# Legend for the symbiont profiles in tracked colonies figures
profile_legend <- symportal_df %>%
  left_join(profile_groups_suppl_df) %>%
  distinct(Profile, Profile_cols, Clade) %>%
  mutate(Profile_cols = ifelse(Clade == 'D', "#EC3F3A", Profile_cols)) %>%
  mutate(Profile_cols = ifelse(is.na(Profile_cols), 'grey50', Profile_cols)) %>%
  arrange(Clade, Profile_cols) %>%
  mutate(Profile = factor(Profile, levels = unique(Profile)))

# Plots data
ggplot() +
  geom_tile(data = profile_legend, 
            aes(x = Profile, y = 1, fill = Profile_cols),
            colour = 'black') + 
  scale_fill_identity() + 
  theme_void() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3))

ggsave('figures_tables/legend_for_S2_S3.pdf')

#####################################################################
