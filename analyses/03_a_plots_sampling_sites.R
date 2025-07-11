#####################################################################
# Last updated 2025-05-23 - Alexander Van Nynatten
# Plots sites surveyed and human disturbance gradient across KI
#####################################################################

#####################################################################
# Loads files and the libraries for analysis
#####################################################################

library(sf)
library(tidyverse)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)

# Bounding box for Kiritimati
KI_bbox <- st_bbox(c(xmin = -157.599, xmax = -157.139, ymin = 1.654, ymax = 2.054)) %>%
  st_as_sfc()

# List of sites with coral surveys conducted for expeditions included in study
coral_survey_sites <- read_csv('data/02_coralnet_results_df.csv') %>%
  distinct(site_name) %>%
  mutate(site_name = as.character(site_name)) %>%
  .$site_name

# DNA sample metadata
samples_df <- read_csv("data/sample_metadata.csv")

# Location data for the sampling sites
sites_df <- read_csv("data/survey_metadata.csv") %>%
  distinct(site_name, pub.name, lat, lon, continous.pressure.2km) %>%
  mutate(sqr_pressure = sqrt(continous.pressure.2km)) %>%
  filter(site_name %in% coral_survey_sites) %>%
  mutate(DNA_site = ifelse(site_name %in% samples_df$site_name, 'DNA', 'Survey only')) %>%
  st_as_sf(coords = c("lon", "lat"), 
           crs = st_crs('WGS84'))

# Shape files from GADM version 1.0, in March 2009
KI_ne10m_sf <- read_sf('data/KI_shape_files/diva-gis/KIR_adm0.shp') %>%
  st_crop(KI_bbox)

# Location and populations of villages on KI
villages_df <- data.frame(
  lat = c(1.989386333, 2.022090868, 1.983594483, 1.865048333),
  lon = c(-157.4760637, -157.4884092, -157.3683462, -157.5522183),
  pop = c(1879, 2311, 955, 441),
  village = c("London", "Tabwakea", "Banana", "Poland")
  ) %>%
    st_as_sf(coords = c("lon", "lat"), 
             crs = st_crs('WGS84'))

#####################################################################
# Plots simple map for main text
#####################################################################

ggplot() +
  geom_sf(data = KI_ne10m_sf, 
          fill = 'lightgrey', colour = NA) +	
  geom_sf(data = villages_df, 
          aes(size = pop),
          colour = 'grey3', shape = 21) +
  geom_sf(data = sites_df %>%
            filter(DNA_site == 'DNA'), 
          aes(colour = sqr_pressure),
          size = 1.5, shape = 18) +
  geom_sf(data = sites_df, 
          aes(colour = sqr_pressure),
          size = 2, shape = 22) +
  theme_test(base_size = 10) + 
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box.margin = margin(l = 10, unit = "mm")
  ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +

  scale_size_continuous(range = c(1,6)) + 
  scale_colour_gradient2(low = "#0096C7", mid = '#6E9075', high = "#825F48", midpoint = 40, na.value = NA)

ggsave('figures_tables/Fig1a_KI_map.pdf', width = 100, height = 80, units = 'mm')

 kiritimati <- data.frame(
  lon = -157.36,
  lat = 1.86
)

# Get a world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the map
ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(-180, -100), ylim = c(-20, 50), expand = FALSE) +
  geom_point(data = kiritimati, aes(x = lon, y = lat), 
             color = "grey2", size = 2, shape = 22) + 
  geom_point(data = kiritimati, aes(x = lon, y = lat + 4), 
             fill = "grey2", size = 1, shape = 25) + 
  
  theme_void() + 
  theme(panel.border = element_rect(fill = NA))

ggsave('figures_tables/Fig1a_KI_inset.pdf', width = 30, height = 30, units = 'mm')


#####################################################################
# Plots detailed map for supplementary materials
#####################################################################

ggplot() +
  geom_sf(data = KI_ne10m_sf, 
          fill = 'lightgrey', colour = NA) +	
  geom_sf(data = villages_df, 
          aes(size = pop),
          colour = 'grey3', shape = 21) +
  geom_sf(data = sites_df %>%
            filter(DNA_site == 'DNA'), 
          aes(colour = sqr_pressure),
          size = 3.5, shape = 18) +
  geom_sf(data = sites_df, 
          aes(colour = sqr_pressure),
          size = 4, shape = 22) +
  geom_sf_text(data = sites_df,
          aes(label = pub.name),
          size = 4, colour = 'grey2') +
  theme_test(base_size = 10) + 
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box.margin = margin(l = 10, unit = "mm")
  ) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  scale_size_continuous(range = c(1,6)) + 
  scale_colour_gradient2(low = "#0096C7", mid = '#6E9075', high = "#825F48", midpoint = 40, na.value = NA)

ggsave('figures_tables/FigS1_detailed_map.pdf', width = 180, height = 180, units = 'mm')

#####################################################################