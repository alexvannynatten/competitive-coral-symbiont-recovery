#####################################################################
# Last updated 2025-05-23 - Alexander Van Nynatten
# Plots SST and DHW at KI over the period surveyed
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(tidyverse)

# Location data for the sampling sites
sites_df <- read_csv("data/survey_metadata.csv") %>%
  distinct(site_name, pub.name, lat, lon)

# SST and DHW data from NOAA CRW
CRW_df <- read_csv("data/dhw_5km_97a5_e43d_e9ef.csv") # From https://pae-paha.pacioos.hawaii.edu/erddap/griddap/dhw_5km.html

# A dataframe for the expedition dates
exped_labels_df <- data.frame(exp_name = c('2013','2014','2015a','2015b',
                                           '2015c','2016a',
                                           '2016b','2017','2018','2019', 
                                           '2023b'),
                              exp_date = c('2013-08-01', '2014-08-01','2015-01-01','2015-05-01',
                                           '2015-07-01','2016-03-01',
                                           '2016-11-01', '2017-07-01','2018-07-01','2019-07-01',
                                           '2023-09-01')) %>%
  mutate(Measure = 'CRW_DHW_avg')

#####################################################################
# Subsets NOAA SST data to include only measures nearest the sampling sites
#####################################################################

CRW_df <- CRW_df[2:nrow(CRW_df), ] %>%
  mutate(latitude = as.numeric(latitude), 
         longitude = as.numeric(longitude),
         CRW_DHW = as.numeric(CRW_DHW),
         CRW_SST = as.numeric(CRW_SST),
         time = ymd_hms(time)
  ) %>%
  mutate(latlong = paste(latitude, longitude)) %>% 
  mutate(date = as.Date(time))

# Function to find the nearest lat/long value
find_nearest <- function(coord, values) {
  abs_diff <- abs(values - coord)
  return(which.min(abs_diff))
}

# Gets the closest locations in the SST to each of the sites sampled
CRW_lat_long_list <- sites_df %>%
  group_by(site_name) %>%
  mutate(latitude = CRW_df[sapply(lat, find_nearest, values = CRW_df$latitude), ]$latitude) %>%
  mutate(longitude = CRW_df[sapply(lon, find_nearest, values = CRW_df$longitude), ]$longitude) %>%
  ungroup() %>%
  distinct(latitude, longitude) %>%
  mutate(latlong = paste(latitude, longitude))

# Calculates daily averages across the subsampled points in the SST data
CRW_df_df_long <- CRW_df %>%
  filter(latlong %in% CRW_lat_long_list$latlong) %>%
  group_by(date) %>%
  summarize(CRW_DHW_avg = mean(CRW_DHW, na.rm = TRUE), 
            CRW_SST_avg = mean(CRW_SST, na.rm = TRUE)) %>%
  pivot_longer(-date, names_to = 'Measure', values_to = 'Temperature')

################################################################################
# Summarizing
################################################################################

# Subsets by DHW but also keeps track of if it is the first or second heatwave to prevent tearing in the ribbon
CRW_df_df_long_plot <- CRW_df_df_long %>%
  filter(Measure == 'CRW_DHW_avg') %>%
  mutate(DHW_fill_color = cut(Temperature, breaks = c(-Inf, 4, 8, 12, 16, 20, Inf), 
                              labels = c(0, 1, 2, 3, 4, 5))) %>%
  mutate(DHW_fill_color_num = as.numeric(DHW_fill_color)) %>%
  mutate(DHW_hump_num = cumsum(DHW_fill_color_num > lag(DHW_fill_color_num, default = first(DHW_fill_color_num)))) %>% # Counts which heatwave
  mutate(DHW_plot_ribbon = paste(DHW_hump_num, DHW_fill_color))

# New dataframe for the plot
fig1c_plot_df <- CRW_df_df_long %>%
  left_join(CRW_df_df_long_plot)

# Plots the SST and DHW as lines 
p <- ggplot(fig1c_plot_df, aes(x = date)) +
  geom_line(aes(y = Temperature, colour = Measure), linewidth = 0.5) +
  scale_y_continuous(name = "Max Degree Heating Weeks") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  theme_minimal(base_size=13) +
  theme(axis.title.x = element_blank(), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour='NA'),
        panel.background = element_rect(fill = "white", colour='NA'),
        axis.line = element_line(color = "darkgrey"),
        axis.ticks = element_line(color = "darkgrey")) +
  facet_grid(Measure ~ ., scales = 'free_y') + 
  scale_colour_manual(values = c('grey2', 'dodgerblue'))

# Plots the expeditions as vertical lines and labels i through x
p <- p + geom_text(data = exped_labels_df,
                   aes(x = as.Date(exp_date), y = 35, label = exp_name),
                   colour = 'grey8', angle = 90) +
  geom_point(data = exped_labels_df, 
             aes(x = as.Date(exp_date), y = 32), 
             fill = "grey3", shape = 25)

# Loops through DHW groupings made in CRW_df_df_long_plot and adds each grouping separately to the plot to avoid tearing
for(i in na.omit(unique((fig1c_plot_df$DHW_plot_ribbon)))) {
  p <- p +
    geom_ribbon(data = CRW_df_df_long_plot %>%
                  filter(DHW_plot_ribbon == i),
                aes(ymin = 0, ymax = Temperature, fill = DHW_fill_color))
}

alert_cols <- c('#f3e79b','#fac484','#f8a07e','#eb7f86','#ce6693','#a059a0')

# Adds the colour scheme for bleaching alert levels and fixes the legends
p + scale_fill_manual(values = alert_cols) +
  geom_hline(yintercept = 28.1366, colour = 'firebrick2') +
  theme(legend.position = "bottom") +  # Move legend to bottom
  guides(fill = guide_legend(nrow = 1)) +  # Keep legend in a single row
  labs(fill = "Bleaching alert level") + 
  scale_x_date(
    limits = as.Date(c('2013-01-01', '2024-01-01')),   # Set x-axis limits
    breaks = seq(as.Date('2012-01-01'), as.Date('2024-01-01'), by = '2 year'),
    date_labels = "%Y"
  ) + 
  theme_test(base_size = 10) +
  theme(legend.position = 'bottom')


ggsave("figures_tables/Fig1c_DHWSST_plot.pdf", width = 190, height = 120, units = "mm")

################################################################################
