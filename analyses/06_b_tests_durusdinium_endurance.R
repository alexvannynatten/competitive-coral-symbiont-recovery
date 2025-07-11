#####################################################################
# Last updated 2025-05-25 - Alexander Van Nynatten
# Tests differences in durusdinium hosted before and after heatwave
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Tidied symportal output
symportal_df <- read_csv('data/01_symportal_results_df.csv')

#####################################################################
# Dataframe of read frequency of Durusdinium
#####################################################################

# Only at sites sampled before and after the heatwave
before_after <- symportal_df %>%
  filter(expedition_class %in% c('1_before', '4_longafter')) %>%
  distinct(expedition_class, site_name) %>%
  group_by(site_name) %>%
  summarize(n = n()) %>%
  filter(n > 1)

# Makes dataframe of the relative read frequency of Durusdinium
D_rel_freq_df <- symportal_df %>%
  filter(site_name %in% before_after$site_name) %>%
  filter(host_genus == 'Pocillopora') %>%
  filter(Clade %in% c('D', 'C')) %>%
  select(sample_name, expedition, Clade, level, Count_freq, expedition_class, host_genus) %>%
  complete(Clade, nesting(sample_name, host_genus, level,expedition_class,  expedition), fill = list(Count_freq = 0)) %>% # Gives any sample without detectable Durusdinium
  filter(Clade == 'D') %>%
  filter(level %in% c(1,3,5))

# Only samples before and long after
D_rel_freq_df_plot <- D_rel_freq_df %>%
  filter(expedition_class %in% c('1_before', '4_longafter')) %>%
  filter(!is.na(expedition_class)) %>%
  mutate(expedition_plot = recode(expedition_class,
                                  "1_before" = "2013-2015b",
                                  "4_longafter" = "2023b")) %>%
  mutate(level_plot = recode(level,
                                  '1' = "Very low",
                                  '3' = "Medium",
                                  '5' = "Very high")) %>%
  mutate(level_plot = factor(level_plot, levels = c("Very low", "Medium", "Very high" )))

#####################################################################
# Tests difference in Durusdinium RRA
#####################################################################

# Calculates p values for the comparisons of Durusdinium relative read frequency before and after heatwave
D_wilcox_t <- data.frame()
for(i in c(1,3,5)) {
  
  # Loops through each disturbance category and conducts the test
  test_out <- wilcox.test(Count_freq~expedition_class,
                          exact = FALSE, paired = FALSE,
                          data = D_rel_freq_df_plot %>%
                            filter(level == i))
  
  # Adds to table of p values
  temp_out_t <- data.frame(Disturbance = i, 
                           pval = test_out$p.value,
                           W = test_out$statistic
  )
  
  D_wilcox_t <- rbind(D_wilcox_t, temp_out_t)
  
}

# Adjusts for multiple testing
D_wilcox_t$pval_adj <- p.adjust(D_wilcox_t$pval, method = 'bonferroni', n = nrow(D_wilcox_t))

D_wilcox_t <- D_wilcox_t %>%
  mutate(pval_adj_plot = paste0('p = ', round(pval_adj, 3))) %>%
  mutate(level_plot = recode(Disturbance,
                             '1' = "Very low",
                             '3' = "Medium",
                             '5' = "Very high")) %>%
  mutate(level_plot = factor(level_plot, levels = c("Very low", "Medium", "Very high" )))

#####################################################################
# Plots comparison of Durusdinium RRA
#####################################################################

# Plots relative read frequency differences 
ggplot() + 
  geom_violin(data = D_rel_freq_df_plot,
              aes(expedition_plot, Count_freq, fill = level_plot),
              scale = 'width', draw_quantiles = c(0.5), adjust = 1) +
  geom_jitter(data = D_rel_freq_df_plot,
              aes(expedition_plot, Count_freq, colour = level_plot),
              shape = 21, width = 0.1, size = 2.5) + 
  geom_text(data = D_wilcox_t,
            aes(x = 1.5, y = 1.2, label = pval_adj_plot)) + 
  theme_test(base_size = 10) + 
  facet_grid(. ~ level_plot) + 
  scale_fill_manual(values = c('#0096C7', '#6E9075', '#825F48')) + 
  scale_colour_manual(values = c('#0096C7', '#6E9075', '#825F48')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('figures_tables/Fig4b_durusdinium_abudnance.pdf', height = 100, width = 120, units = "mm")

#####################################################################




