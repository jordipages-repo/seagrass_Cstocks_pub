# # # # # #
# Analysing the data set provided by Hilary Kennedy
# Jordi F. Pagès
# 20-03-2019
# University of Barcelona
# # # # # # 

# # # # # # # # # # # # # # # # # # # # # 
# Loading the data clean and corrected  #
# # # # # # # # # # # # # # # # # # # # # 

source("01_DataImport&Corrections_CarbonReview.R")
library(gridExtra)
library(patchwork)
# library(cowplot)

# # # # # # # # # # # # # # # # # # # # # 
# Data visualisation of Carbon stocks ----
# # # # # # # # # # # # # # # # # # # # # 


# HISTOGRAM OF CSTOCKS ----

# Histogram for 0-20
cstocks %>% 
  ggplot() +
  geom_histogram(aes(Stock_0_20cm), binwidth = 10, fill = "#9FDA3AFF") +
  xlab("Carbon stock 0-20 cm") + 
  theme_bw()
# ggsave("Figs/Cstocks_histogram_0_20cm.pdf")

# Histogram for 20-50
cstocks %>% 
  ggplot() +
  geom_histogram(aes(Stock_20_50cm), binwidth = 10, fill = "#9FDA3AFF") +
  xlab("Carbon stock 20-50 cm") + 
  theme_bw()
# ggsave("Figs/Cstocks_histogram_20_50cm.pdf")

# Histograms all together in the same plot (for that we need to gather (tidy) the data set)
cstocks_tidy <- cstocks %>% 
  select(-Stock_20_50cm_estimate) %>% 
  gather(key = depth, value = cstocks, Stock_0_20cm:Stock_20_50cm)

cstocks_tidy %>% 
  filter(depth != "Stock_0_50cm") %>% 
  ggplot() +
  geom_histogram(aes(cstocks, fill = depth, colour = depth), alpha = 0.5) +
  theme_bw()
# ggsave("Figs/Cstocks_histogram_allDepths_same_plot.pdf")



# BARPLOT PER SPECIES (COUNT) ----
cstocks %>%
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev()) %>% 
  filter(Meadow_type == "monospecific") %>%
  group_by(Species) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  geom_bar(aes(x = Species, y = n), stat = "identity") +
  geom_text(aes(x = Species, y = n, label = str_c("(", n, ")")), nudge_y = 6, size = 3) +
  ylab("count") +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text.align = 0,
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "italic"))
# ggsave("Figs/CoreCount_by_Species.pdf")



# BARPLOT PER SPECIES PERCENT ----
cstocks %>%
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev()) %>%
  filter(Meadow_type == "monospecific") %>%
  group_by(Species) %>% 
  summarise(n = n()) %>% 
  mutate(percent = 100*(n/sum(n))) %>% 
  ggplot() +
  geom_bar(aes(x = Species, y = percent), stat = "identity") +
  geom_text(aes(x = Species, y = percent, label = str_c(round(percent), "%")), nudge_y = 1, size = 3) +
  ylab("%") +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text.align = 0,
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "italic"))
# ggsave("Figs/CoreCountPercent_by_Species.pdf")



# BARPLOT PER COUNTRY ----
cstocks %>% 
  mutate(Country = factor(Country) %>%  fct_infreq() %>% fct_rev()) %>%
  filter(!is.na(Country)) %>% 
  group_by(Country) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  geom_bar(aes(x = Country, y = n), stat = "identity") +
  geom_text(aes(x = Country, y = n, label = str_c("(", n, ")")), nudge_y = 15, size = 3) +
  ylab("count") +
  coord_flip() +
  theme_bw() + 
  theme(legend.title = element_blank(),
      legend.spacing.x = unit(0.2, 'cm'),
      legend.text.align = 0,
      text = element_text(size=14),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank())
# ggsave("Figs/CountryCount.pdf")



# BARPLOT PER COUNTRY PERCENT ----
cstocks %>% 
  mutate(Country = factor(Country) %>%  fct_infreq() %>% fct_rev()) %>%
  filter(!is.na(Country)) %>% 
  group_by(Country) %>% 
  summarise(n = n()) %>% 
  mutate(percent = 100*(n/sum(n))) %>% 
  ggplot() +
  geom_bar(aes(x = Country, y = percent), stat = "identity") +
  geom_text(aes(x = Country, y = percent, label = str_c(round(percent), "%")), nudge_y = 3, size = 3) +
  ylab("%") +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text.align = 0,
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/CountryPercent.pdf")



# MAP OF CORES BY COUNTRY ----
cstocksCountry <- cstocks %>% 
  group_by(Country) %>% 
  filter(!is.na(Country)) %>% 
  summarise(n = n())

unique(cstocksCountry$Country)
map.world <- map_data('world')
unique(map.world$region)

# Checking potential join mismatches
cstocksCountry %>% 
  anti_join(map.world, by = c('Country' = 'region')) %>% 
  print(n = Inf)
# Ok, no mismatches. Useful to see what will not be joined. In this case, nothing. Everything will be joined.

# So we can proceed with the left_join
map.cstocks <- left_join(map.world, cstocksCountry, by = c('region' = 'Country'))

# Map of number of cores per country
ggplot(map.cstocks, aes(x = long, y = lat, group = group )) +
  geom_polygon(aes(fill = n)) +
  # scale_fill_gradientn(colours = brewer.pal(5, "YlOrRd"), trans = "log10", na.value = "#d0d0d0") +
  # scale_fill_gradientn(colours = rev(c("#9FDA3AFF", "#4AC16DFF", "#1FA187FF", "#277F8EFF", "#365C8DFF", 
  # "#46337EFF")), trans = "log10", na.value = "#d0d0d0") +
  # scale_fill_gradientn(colours = brewer.pal(5, "Blues"), trans = "log10") +
  # theme_minimal()
  labs(fill = '', x = NULL, y = NULL) +
  theme(text = element_text(color = '#EEEEEE'), axis.ticks = element_blank(), axis.text = element_blank(), 
        panel.grid = element_blank(), panel.background = element_rect(fill = '#787878'), 
        plot.background = element_rect(fill = '#787878'), legend.position = c(.18,.36) ,
        legend.background = element_blank(), legend.key = element_blank())
# ggsave(filename = "Figs/CountryMap.pdf")




# MAPPING INDIVIDUAL CORES WITH GGMAP  ----
# From https://lucidmanager.org/geocoding-with-ggmap/
library(ggmap)
# Now, to use Google API you have to be registered (include a credit card) and get an API key. See below.
api <- readLines("google.api") # Text file with the API key
register_google(key = api)
getOption("ggmap")
has_google_key()
# To check we don't exceed the day_limit that Google imposes (2500)
geocodeQueryCheck()

#### Successful trial with https://blog.dominodatalab.com/geographic-visualization-with-rs-ggmaps/
# We get the lat lon coordinates in a df
Df_coord <- cstocks %>% 
  select(Latitude, Longitude) %>% 
  group_by(Latitude, Longitude) %>% 
  summarise(Samples = n())

# We build a world map from google
Map <- ggmap(ggmap = get_googlemap(center = c(lon = 10, lat = 0), 
                                   zoom = 1, 
                                   style = c(feature="all",element="labels",visibility="off"),
                                   maptype = "roadmap",
                                   size = c(512, 512), 
                                   scale = 2,
                                   color = "bw"),
             extent = "panel")
Map

# We plot the frequency of samples per coordinate combination as 'bubbles'
Map +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude, size = Samples), col="#365C8DFF", alpha=0.4) +
  scale_size(range = range(sqrt(Df_coord$Samples))) +
  scale_y_continuous(limits = c(-70,70)) +
  theme_bw()
# ggsave("Figs/CoreMapBubbles.pdf")


# Just plotting the points
Map +
  geom_jitter(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col= "#0063ad", alpha=0.2, width = 1, height = 1) +
  scale_y_continuous(limits = c(-55,70), breaks = c(-50, -25, 25, 50, 75)) +
  theme_bw() +
  xlab("Longitude (º)") +
  ylab("Latitude (º)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))
ggsave("Figs/Final_Figures_March2021/GlobalCoreMapPointsNEW_sized.png", width = 180, height = 110, units = "mm", dpi = 350)



# MAPPING INDIVIDUAL CORES IN GGMAP ZOOMING IN EACH REGION ----
# Europe
MapEurope <- ggmap(ggmap = get_googlemap(center = c(lon = 10, lat = 50), 
                                   zoom = 4, 
                                   maptype = "terrain",
                                   size = c(512, 512), 
                                   scale = 2,
                                   color = "bw"),
             extent = "panel")

Cstocks_Europe <- MapEurope +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col="#0063ad", size = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(35,60)) +
  scale_x_continuous(limits = c(-10,30)) +
  theme_bw()

# N.America
MapAmerica <- ggmap(ggmap = get_googlemap(center = c(lon = -100, lat = 40), 
                                         zoom = 3, 
                                         maptype = "terrain",
                                         size = c(512, 512), 
                                         scale = 2,
                                         color = "bw"),
                   extent = "panel")

Cstocks_America <- MapAmerica +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col="#0063ad", size = 2, alpha = 0.4) +
  # scale_y_continuous(limits = c(15,50)) +
  # scale_x_continuous(limits = c(-125,-65)) +
  theme_bw()

# Brasil
MapBrasil <- ggmap(ggmap = get_googlemap(center = c(lon = -50, lat = -20), 
                                         zoom = 4, 
                                         maptype = "terrain",
                                         size = c(512, 512), 
                                         scale = 2,
                                         color = "bw"),
                   extent = "panel")

Cstocks_Brasil <- MapBrasil +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col="#0063ad", size = 2, alpha = 0.4) +
  # scale_y_continuous(limits = c(-30, -5)) +
  # scale_x_continuous(limits = c(-70,-30)) +
  theme_bw()

# Arabia
MapArabia <- ggmap(ggmap = get_googlemap(center = c(lon = 45, lat = 25),
                                         zoom = 5,
                                         maptype = "terrain",
                                         size = c(512, 512),
                                         scale = 2,
                                         color = "bw"),
                   extent = "panel")

Cstocks_Arabia <- MapArabia +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col="#0063ad", size = 2, alpha = 0.4) +
  # scale_y_continuous(limits = c(-30, -5)) +
  # scale_x_continuous(limits = c(-70,-30)) +
  theme_bw()

# Australia
MapAustralia <- ggmap(ggmap = get_googlemap(center = c(lon = 133, lat = -28),
                                         zoom = 4,
                                         maptype = "terrain",
                                         size = c(512, 512),
                                         scale = 2,
                                         color = "bw"),
                   extent = "panel")

Cstocks_Australia <- MapAustralia +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col="#0063ad", size = 2, alpha = 0.4) +
  # scale_y_continuous(limits = c(-30, -5)) +
  # scale_x_continuous(limits = c(-70,-30)) +
  theme_bw()

# Asia
MapAsia <- ggmap(ggmap = get_googlemap(center = c(lon = 120, lat = 25),
                                            zoom = 3,
                                            maptype = "terrain",
                                            size = c(512, 512),
                                            scale = 2,
                                            color = "bw"),
                      extent = "panel")

Cstocks_Asia <- MapAsia +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude), col="#0063ad", size = 2, alpha = 0.4) +
  # scale_y_continuous(limits = c(-30, -5)) +
  # scale_x_continuous(limits = c(-70,-30)) +
  theme_bw()


# General map
# zooms <- grid.arrange(Cstocks_America, 
#                    Cstocks_Europe, 
#                    Cstocks_Asia,
#                    Cstocks_Brasil,
#                    Cstocks_Arabia,
#                    Cstocks_Australia,
#                    nrow = 2, top = "")
patchwork <- Cstocks_America + Cstocks_Europe + Cstocks_Asia + Cstocks_Brasil + Cstocks_Arabia + Cstocks_Australia + plot_annotation(tag_levels = 'A')
patchwork & xlab("Longitude") & ylab("Latitude")
# ggsave("Figs/Core_coordinates_MapZoom_Points.pdf")



# CSTOCKS PER SPECIES FACET PLOT ----
source('reorder_within_function.R')
cstocks_tidy <- cstocks %>% 
  select(-Stock_20_50cm_estimate) %>% 
  gather(key = depth, value = cstocks, Stock_0_20cm:Stock_20_50cm)

# ggplot(iris_gathered, aes(reorder_within(Species, value, metric), value)) +
  #'   geom_boxplot() +
  #'   scale_x_reordered() +
  #'   facet_wrap(~ metric, scales = "free_x")

### WORKING ON THIS!!! NOW HAVE TO REORDER
cstocks_tidy %>% 
  filter(Meadow_type == "monospecific") %>%
  ggplot() +
  geom_boxplot(aes(x = Species, y = cstocks), fill = "#9FDA3AFF") +
  scale_y_reordered() +
  coord_flip() +
  xlab("Species") +
  ylab("Carbon stock") +
  facet_grid(~depth, scales = "free_x") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Species_facet_depths.pdf", width = 12, height = 6)


# CSTOCKS PER SPECIES grid.arrange PLOT ----

p1 <- cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, Stock_0_20cm, FUN = median), y = Stock_0_20cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Species") +
  ylab("Carbon stock 0-20 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))

p2 <- cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(Stock_20_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, Stock_20_50cm, FUN = median), y = Stock_20_50cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("") +
  ylab("Carbon stock 20-50 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))

p3 <- cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(Stock_0_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, Stock_0_50cm, FUN = median), y = Stock_0_50cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("") +
  ylab("Carbon stock 0-50 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))

ml <- grid.arrange(p1, p2, p3, widths = c(2,2,2), nrow = 1, top = "")
# ggsave("Figs/Cstocks_by_Species_grid.arrange_depths.pdf", plot = ml, width = 12, height = 6)



# CSTOCKS PER SPECIES ----
# Boxplots cstocks 0-20 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, Stock_0_20cm, FUN = median), y = Stock_0_20cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Species") +
  ylab("Carbon stock 0-20 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Species_0_20cm.pdf")

# Boxplots cstocks 0-50 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(Stock_0_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, Stock_0_50cm, FUN = median), y = Stock_0_50cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Species") +
  ylab("Carbon stock 0-50 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Species_0_50cm.pdf")


# Boxplots cstocks 20-50 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(Stock_20_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, Stock_20_50cm, FUN = median), y = Stock_20_50cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Species") +
  ylab("Carbon stock 20-50 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Species_20_50cm.pdf")




# CSTOCKS PER GENUS ----

# Boxplots cstocks 0-20 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(Genus, Stock_0_20cm, FUN = median), y = Stock_0_20cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Genus") +
  ylab("Carbon stock 0-20 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Genus_0_20cm.pdf")

# Boxplots cstocks 0-50 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(Stock_0_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Genus, Stock_0_50cm, FUN = median), y = Stock_0_50cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Genus") +
  ylab("Carbon stock 0-50 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Genus_0_50cm.pdf")

# Boxplots cstocks 20-50 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(Stock_20_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Genus, Stock_20_50cm, FUN = median), y = Stock_20_50cm), fill = "#9FDA3AFF") +
  coord_flip() +
  xlab("Genus") +
  ylab("Carbon stock 20-50 cm") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Cstocks_by_Genus_20_50cm.pdf")



# CSTOCKS PER MEADOW TYPE ----
# Boxplots cstocks 0-20 cm
cstocks %>% 
  ggplot() +
  geom_boxplot(aes(x = Meadow_type, y = Stock_0_20cm), fill = "#9FDA3AFF") +
  xlab("Meadow type") +
  ylab("Carbon stock 0-20 cm") +
  coord_flip() +
  theme_bw()
# ggsave("Figs/Cstocks_meadowtype_0_20cm.pdf")

# Boxplots cstocks 0-50 cm
cstocks %>% 
  ggplot() +
  geom_boxplot(aes(x = Meadow_type, y = Stock_0_50cm), fill = "#9FDA3AFF") +
  xlab("Meadow type") +
  ylab("Carbon stock 0-50 cm") +
  coord_flip() +
  theme_bw()
# ggsave("Figs/Cstocks_meadowtype_0_50cm.pdf")

# Boxplots cstocks 20-50 cm
cstocks %>% 
  ggplot() +
  geom_boxplot(aes(x = Meadow_type, y = Stock_20_50cm), fill = "#9FDA3AFF") +
  xlab("Meadow type") +
  ylab("Carbon stock 20-50 cm") +
  coord_flip() +
  theme_bw()
# ggsave("Figs/Cstocks_meadowtype_20_50cm.pdf")



# DELTA13C PER SPECIES ----
# Boxplots delta13C 0-20 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(d13C_0_20cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, d13C_0_20cm, FUN = median), y = d13C_0_20cm), fill = "#1FA187FF") +
  coord_flip() +
  xlab("Species") +
  ylab(bquote(delta*"13C 0-20 cm")) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/delta13C_by_Species_0_20cm.pdf")

# Boxplots delta13C 0-50 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(d13C_0_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, d13C_0_50cm, FUN = median), y = d13C_0_50cm), fill = "#1FA187FF") +
  coord_flip() +
  xlab("Species") +
  ylab(bquote(delta*"13C 0-50 cm")) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/delta13C_by_Species_0_50cm.pdf")


# Boxplots delta13C 20-50 cm
cstocks %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(d13C_20_50cm)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, d13C_20_50cm, FUN = median), y = d13C_20_50cm), fill = "#1FA187FF") +
  coord_flip() +
  xlab("Species") +
  ylab(bquote(delta*"13C 20-50 cm")) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
 # ggsave("Figs/delta13C_by_Species_20_50cm.pdf")




# # # # # # # # # # # # # # # # # # # # # 
# Data visualisation of Plant traits ----
# # # # # # # # # # # # # # # # # # # # # 

# ABOVEGROUND BIOMASS PER SPECIES ----
# Dotplots mean aboveground biomass
cstocks_traits %>% 
  filter(!is.na(Mean_aboveground_biomass)) %>% 
  ggplot() +
  geom_point(aes(x = reorder(Species, Mean_aboveground_biomass, FUN = mean),
                   y = Mean_aboveground_biomass), fill = "#1FA187FF") +
  geom_errorbar(aes(x = reorder(Species, Mean_aboveground_biomass, FUN = mean),
                    ymin=Mean_aboveground_biomass-SE_Mean_aboveground_biomass,
                    ymax=Mean_aboveground_biomass+SE_Mean_aboveground_biomass),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(x = reorder(Species, Mean_aboveground_biomass, FUN = mean),
                y = Mean_aboveground_biomass+SE_Mean_aboveground_biomass, 
                label = str_c("n = ", N_Mean_aboveground_biomass)),
            nudge_y = 50, size = 3) +
  coord_flip() +
  xlab("Species") +
  ylab("Mean aboveground biomass") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/MeanAbovegroundB_by_Species_.pdf")


# Dotplots max aboveground biomass
cstocks_traits %>% 
  filter(!is.na(Max_aboveground_biomass)) %>% 
  ggplot() +
  geom_point(aes(x = reorder(Species, Max_aboveground_biomass, FUN = mean),
                 y = Max_aboveground_biomass), fill = "#1FA187FF") +
  geom_errorbar(aes(x = reorder(Species, Max_aboveground_biomass, FUN = mean),
                    ymin=Max_aboveground_biomass-SE_Max_aboveground_biomass,
                    ymax=Max_aboveground_biomass+SE_Max_aboveground_biomass),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(x = reorder(Species, Max_aboveground_biomass, FUN = mean),
                y = Max_aboveground_biomass+SE_Max_aboveground_biomass, 
                label = str_c("n = ", N_Max_aboveground_biomass)),
            nudge_y = 100, size = 3) +
  coord_flip() +
  xlab("Species") +
  ylab("Max aboveground biomass") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/MaxAbovegroundB_by_Species_.pdf")



# BELOWGROUND BIOMASS PER SPECIES ----
# Dotplots max aboveground biomass
cstocks_traits %>% 
  filter(!is.na(Mean_belowground_biomass)) %>% 
  ggplot() +
  geom_point(aes(x = reorder(Species, Mean_belowground_biomass, FUN = mean),
                 y = Mean_belowground_biomass), fill = "#1FA187FF") +
  geom_errorbar(aes(x = reorder(Species, Mean_belowground_biomass, FUN = mean),
                    ymin=Mean_belowground_biomass-SE_Mean_belowground_biomass,
                    ymax=Mean_belowground_biomass+SE_Mean_belowground_biomass),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(x = reorder(Species, Mean_belowground_biomass, FUN = mean),
                y = Mean_belowground_biomass+SE_Mean_belowground_biomass, 
                label = str_c("n = ", N_Mean_belowground_biomass)),
            nudge_y = 150, size = 3) +
  coord_flip() +
  xlab("Species") +
  ylab("Mean belowground biomass") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/MeanBelowgroundB_by_Species_.pdf")

# Dotplots mean aboveground biomass
cstocks_traits %>% 
  filter(!is.na(Max_belowground_biomass)) %>% 
  ggplot() +
  geom_point(aes(x = reorder(Species, Max_belowground_biomass, FUN = mean),
                 y = Max_belowground_biomass), fill = "#1FA187FF") +
  geom_errorbar(aes(x = reorder(Species, Max_belowground_biomass, FUN = mean),
                    ymin=Max_belowground_biomass-SE_Max_belowground_biomass,
                    ymax=Max_belowground_biomass+SE_Max_belowground_biomass),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(x = reorder(Species, Max_belowground_biomass, FUN = mean),
                y = Max_belowground_biomass, 
                label = str_c("n = ", N_Max_belowground_biomass)),
            nudge_y = 550, size = 3) +
  coord_flip() +
  xlab("Species") +
  ylab("Max belowground biomass") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/MaxBelowgroundB_by_Species_.pdf")



# ROOT:SHOOT RATIO PER SPECIES ----
cstocks_traits %>% 
  filter(!is.na(Root_shoot_ratio)) %>% 
  ggplot() +
  geom_point(aes(x = reorder(Species, Root_shoot_ratio, FUN = mean),
                 y = Root_shoot_ratio), fill = "#1FA187FF") +
  geom_errorbar(aes(x = reorder(Species, Root_shoot_ratio, FUN = mean),
                    ymin=Root_shoot_ratio-SE_Root_shoot_ratio,
                    ymax=Root_shoot_ratio+SE_Root_shoot_ratio),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(x = reorder(Species, Root_shoot_ratio, FUN = mean),
                y = Root_shoot_ratio+SE_Root_shoot_ratio, 
                label = str_c("n = ", N_Max_belowground_biomass)),
            nudge_y = 1, size = 3) +
  coord_flip() +
  xlab("Species") +
  ylab("Root:Shoot ratio") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/RootShootRatio_by_Species_.pdf")


#########################
