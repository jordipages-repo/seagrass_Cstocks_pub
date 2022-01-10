# # # # # #
# Analysing the data set provided by Hilary Kennedy
# Jordi F. Pagès
# 18-07-2019
# University of Barcelona
# # # # # # 

library(RColorBrewer)
library(gridExtra)
library(scales)
# devtools::install_github("dkahle/ggmap")
library(ggmap)
library(RgoogleMaps)

# # # # # # # # # # # # # # # # # # # # # 
# Loading the data clean and corrected  #
# # # # # # # # # # # # # # # # # # # # # 

source("01_DataImport&Corrections_CarbonReview.R")

# show_col(hue_pal()(4))

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# MAP WITH ALL CORES COLOURED BY MEAN CSTOCKS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

Df_coord <- cstocks_tidy %>% 
  select(Latitude, Longitude, cstocks) %>% 
  group_by(Latitude, Longitude) %>% 
  summarise(Samples = n(),
            mean_stock = mean(cstocks, na.rm = T)) %>% 
  filter(!is.na(Latitude))

# We build a world map from google
api <- readLines("google.api") # Text file with the API key
register_google(key = api)
getOption("ggmap")
has_google_key()
# To check we don't exceed the day_limit that Google imposes (2500)
Map <- ggmap(ggmap = get_googlemap(center = c(lon = 10, lat = 0), 
                                   zoom = 1, 
                                   maptype = "roadmap",
                                   size = c(512, 512), 
                                   scale = 2,
                                   color = "bw"),
             extent = "panel")

# We plot the frequency of samples per coordinate combination as 'bubbles'
Map +
  geom_point(data=Df_coord, mapping = aes(x=Longitude, y=Latitude, colour = mean_stock), alpha = 1) +
  scale_colour_continuous(type = "viridis", trans = "sqrt", name = bquote(atop('Carbon stock', '(Mg C' ~ha^-1* ')')), 
                          breaks=c(5,25,75,150),labels=c(5,25,75,150),
                          limits=c(2,150)) +
  scale_y_continuous(limits = c(-70,70)) +
  theme_bw()
# ggsave(filename = "Figs/Core_coordinates_Map_Points_Coloured_by_meanStocks.pdf")


# # # # # # # # # # # # # # # # # # # # # # # # # # #
# MAP WITH ALL CORES COLOURED BY MEADOW TYPE ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

Df_coord2 <- cstocks_tidy %>% 
  select(Latitude, Longitude, Meadow_type) %>% 
  group_by(Latitude, Longitude, Meadow_type) %>% 
  summarise(Samples = n()) %>% 
  filter(!is.na(Latitude))

# We plot the frequency of samples per coordinate combination as 'bubbles'
Map +
  geom_point(data=Df_coord2, mapping = aes(x=Longitude, y=Latitude, shape = Meadow_type), alpha = 1) +
  scale_shape_manual(values = c(3,16), name = "") +
  scale_y_continuous(limits = c(-70,70)) +
  theme_bw() +
  theme(legend.position = c(0.82,0.1),
        legend.background = element_rect(fill = "transparent"))
# ggsave(filename = "Figs/Core_coordinates_Map_Points_Coloured_by_meadowType_shape.pdf")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Frequency distribution ALL DEPTHS ALL TYPE OF MEADOWS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Either a main graph showing all stocks in reasonably broad bins and an insert graph with 
# stocks minus P.oceanica or else what might be better is to have a stacked bar chart with P. oceanica
# stocks in a different fill to the other species, or maybe as a separate data set plotted on the same graph. 

cstocks_tidy %>% 
  ggplot() +
  geom_histogram(aes(cstocks, fill = Posi), binwidth = 0.5, position = "stack", colour = "black", size = 0.2) +
  xlab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  # scale_y_continuous(limits = c(0,270)) +
  # scale_x_continuous(limits = c(-5,270)) +
  scale_fill_manual(labels = c(expression(italic("Posidonia oceanica")), "Monospecific", "Mixed"),
                    values = c("#1f78b4", "#31a354", "#b2df8a")) +
  # facet_wrap(~depth) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.4),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text.align = 0,
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave(filename = "Figs/Cstocks_histogram_Posi_Mixed_Mono.pdf")

######


# # # # # # # # # # # # # # # # # # # #
# C stocks along latitude by depth ----
# # # # # # # # # # # # # # # # # # # #

# Scatter plot version
cstocks_tidy %>% 
  mutate(Latitude = round((Latitude))) %>% 
  group_by(Latitude, cstocks, Posi, depth) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  # geom_point(aes(y = cstocks, x = Latitude, colour = Posi), alpha = 0.5, size = 2) +
  # geom_jitter(aes(y = cstocks, x = Latitude, colour = Posi), alpha = 0.5, size = 2, width = 0.1, height = 0.1) +
  geom_jitter(aes(y = cstocks, x = Latitude, shape = Posi, colour = Posi), alpha = 0.5, size = 2, width = 0.1, height = 0.1) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  # geom_smooth(aes(y = cstocks, x = Latitude, colour = Posi)) +
  scale_color_manual(labels = c(expression(italic("Posidonia oceanica")), "Monospecific", "Mixed"), 
                     values = c("#1f78b4", "#31a354", "#b2df8a")) +
  scale_shape_manual(labels = c(expression(italic("Posidonia oceanica")), "Monospecific", "Mixed"), 
                     values = c(17, 19, 15)) +
  # facet_wrap(~depth) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text.align = 0,
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_flip()
# ggsave(filename = "Figs/Cstocks_by_latitude_mono_mixed_posi.pdf")


# barplot version
cstocks_tidy %>% 
  filter(!is.na(Latitude)) %>%
  mutate(LatitudeRound = round(Latitude),
         LatitudeBinned = cut(Latitude, breaks = c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70))) %>%
  filter(!is.na(cstocks)) %>% 
  group_by(LatitudeBinned, Posi) %>% 
  summarise(n = n(),
            meanStock = mean(cstocks, na.rm = T),
            std.errorStock = std.error(cstocks)) %>% 
  ggplot() +
  geom_bar(aes(y = meanStock, x = LatitudeBinned, fill = Posi), 
               stat = "identity", position = position_dodge(preserve = "single")) +
  geom_errorbar(aes(y = meanStock, x = LatitudeBinned, fill = Posi, 
                    ymin=meanStock-std.errorStock, ymax=meanStock+std.errorStock),
                width=.3, position = position_dodge(preserve = "single", width = 0.9)) +
  geom_text(aes(x = LatitudeBinned, fill = Posi, ymin = meanStock-std.errorStock, ymax = meanStock+std.errorStock,
                y = ifelse(n>1, meanStock+std.errorStock + 0.5, meanStock + 0.5), label = str_c("(", n, ")")), size = 3,
            position = position_dodge(width = 0.9)) +
  geom_vline(xintercept = 5.5, lty = 2, colour = "grey") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1*')')) +
  xlab("Latitude") + 
  scale_fill_manual(labels = c(expression(italic("Posidonia oceanica")), "Monospecific", "Mixed"), 
                     values = c("#1f78b4", "#31a354", "#b2df8a")) +
  # facet_wrap(~depth) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.position = c(0.8, 0.15),
        legend.text.align = 0,
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_flip()
# ggsave(filename = "Figs/Cstocks_by_latitude_mono_mixed_posiBARPLOTwithNsize.pdf")
# ggsave(filename = "Figs/Cstocks_by_latitude_mono_mixed_posiBARPLOT.pdf")

# # # # # # # # # # # # # # # # # # # #
# C stocks overall mean, species mean, etc ----
# # # # # # # # # # # # # # # # # # # #

# We want to get the means. 
means <- cstocks_tidy %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Meadow_type == "monospecific") %>% 
  group_by(depth) %>% 
  summarise(n = n(),
            mean = mean(cstocks),
            median = median(cstocks),
            max = max(cstocks))


cstocks_tidy %>% 
  filter(!is.na(cstocks)) %>% 
  group_by(depth, Meadow_type) %>% 
  summarise(n = n())
# A tibble: 4 x 3
# Groups:   depth [2]
#   depth    Meadow_type       n
# 1 0-20 cm  monospecific    497
# 2 0-20 cm  multispecific    79
# 3 20-50 cm monospecific    193
# 4 20-50 cm multispecific    33


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# C stocks by species by depth with sample size coloured BARPLOTS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

cstocks_to_order <- cstocks_tidy %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(cstocks)) %>%
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev()) %>% 
  group_by(Species, depth) %>% 
  summarise(n = n(),
            mean_cstocks = mean(cstocks),
            error_cstocks = std.error(cstocks)) %>% 
  ungroup() %>% 
  arrange(depth, mean_cstocks) %>% 
  mutate(order = row_number(),
         nscale = ifelse(n<=10, "n≤10", "n>10"))

ggplot(cstocks_to_order, aes(x = order, y = mean_cstocks, ymin = mean_cstocks-error_cstocks, ymax = mean_cstocks+error_cstocks)) +
  geom_bar(aes(fill = nscale), stat = "identity") +
  scale_fill_manual(name = "Sample size", 
                    labels = c("n > 10", bquote("n" <= 10)),
                    values = c("#00BFC4", "#F8766D")) +
  geom_errorbar(width = 0.3) +
  geom_hline(data = means, aes(yintercept = mean), lty = 2, colour = "darkgrey") +
  geom_text(aes(y = ifelse(n>1, mean_cstocks+error_cstocks + 15, mean_cstocks + 15), label = str_c("(", n, ")")), size = 3) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  facet_wrap(~depth, scales = "free") + 
  coord_flip(ylim = c(0,170)) +
  scale_x_continuous(breaks = cstocks_to_order$order,
                     labels = cstocks_to_order$Species,
                     expand = c(0.008,0.008)) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.justification = c(1,0),
        legend.position = c(0.97, 0.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/Cstocks_by_Species_facet_depth_BARPLOT.pdf")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# C stocks by species by depth with sample size coloured BOXPLOTS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# To count the sample size for each species and depth
listSpeciesN <- cstocks_tidy %>%
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(cstocks)) %>% 
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev()) %>% 
  group_by(Species, depth) %>% 
  summarise(n = n(),
            median_stocks = median(cstocks),
            max_stocks = max(cstocks)) %>% 
  ungroup() %>% 
  arrange(depth, median_stocks) %>% 
  mutate(order = row_number(),
         nscale = ifelse(n<=10, "n≤10", "n>10"))

# We join the above data set with the data set 
cstocksN <- cstocks_tidy %>% 
  filter(Meadow_type == "monospecific") %>%
  filter(!is.na(cstocks)) %>% 
  left_join(listSpeciesN, by = c("Species", "depth"))
  
ggplot(cstocksN, aes(x = as.factor(order), y = cstocks)) +
  geom_boxplot(aes(fill = nscale, colour = nscale), fatten = 1) +
  scale_fill_manual(name = "Sample size", 
                    labels = c("n > 10", bquote("n" <= 10)),
                    values = c("#00BFC4", "#F8766D")) +
  scale_colour_manual(name = "Sample size", 
                    labels = c("n > 10", bquote("n" <= 10)),
                    values = c("#00BFC4", "#F8766D")) +
  geom_point(aes(y = median_stocks), shape = 15, size = 1.5) +
  geom_hline(data = means, aes(yintercept = median), lty = 2, colour = "darkgrey") +
  geom_text(data = listSpeciesN, aes(x = as.factor(order), y = max_stocks + 45, label = str_c("(", n, ")")), size = 3) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  facet_wrap(~depth, scales = "free") + 
  coord_flip(ylim = c(0, 350)) +
  scale_x_discrete(breaks = listSpeciesN$order,
                     labels = listSpeciesN$Species,
                     expand = c(0.025,0.025)) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.justification = c(1,0),
        legend.position = c(0.97, 0.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/Cstocks_by_Species_facet_depth_BOXPLOT.pdf")


# # # # # # # # # # # # # # # # # # # # #
# C stocks by meadow type            ----
# but only comparing multispecific,     #
# with monospecific containing one of   #
# the species in multispecific meadows  #
# # # # # # # # # # # # # # # # # # # # #

# We separate the variable Species into 5 species columns, to be able to know the different species present in multispecific meadows.
cstocksSpeciesSep <- cstocks_tidy %>% 
  separate(Species, into = c("Sp1", "Sp2", "Sp3", "Sp4", "Sp5"), sep = ", ")

# We summarise the complete list of species present in multispecific meadows
listSpeciesMulti <-  cstocksSpeciesSep %>% 
  filter(Meadow_type == "multispecific") %>% 
  select(Sp1:Sp5) %>% 
  gather(key = "Species", value = "value") %>% 
  group_by(value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value))

listSpeciesMulti <- listSpeciesMulti$value
listSpeciesMulti <- rbind(listSpeciesMulti, "Amphibolis spp", "Posidonia spp")

cstocksFilteredMonoMulti <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti)

cstocksFilteredMonoMulti2 <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti) %>% 
  filter(!is.na(cstocks)) %>% 
  group_by(Meadow_type, depth) %>% 
  summarise(n = n(),
            mean_stocks = mean(cstocks),
            error_stocks = std.error(cstocks),
            max_stocks = max(cstocks))

# Boxplots cstocks 
cstocksFilteredMonoMulti %>% 
  ggplot() +
  geom_boxplot(aes(x = Meadow_type, y = cstocks, fill = Meadow_type)) +
  scale_fill_manual(values = c("#31a354", "#b2df8a")) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = Meadow_type, y = max_stocks + 20, label = str_c("(", n, ")"))) +
  scale_x_discrete(labels = c("Monospecific", "Multispecific")) +
  xlab("Meadow type") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/Cstocks_meadowtype_by_depth.pdf")

# Barplots
cstocksFilteredMonoMulti2 %>% 
  ggplot(aes(x = Meadow_type, y = mean_stocks, ymin = mean_stocks-error_stocks, ymax = mean_stocks+error_stocks)) +
  geom_bar(aes(fill = Meadow_type), stat = "identity") +
  geom_errorbar(width = 0.3) +
  geom_text(aes(y = mean_stocks+error_stocks + 10, label = str_c("(", n, ")"))) +
  scale_fill_manual(values = c("#31a354", "#b2df8a")) +
  geom_errorbar(width=.3, position = position_dodge(preserve = "single", width = 0.9)) +
  xlab("") +
  scale_x_discrete(labels = c("Monospecific", "Multispecific")) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=14))
# ggsave(filename = "Figs/Cstocks_meadowtype_by_depth_barplot.pdf")


# # # # # # # # # # # # # # # # # # # #
# C stocks tropical vs. temperate seagrasses ----
# # # # # # # # # # # # # # # # # # # #

# We separate the variable Species into 5 species columns, to be able to know the different species present in multispecific meadows.
cstocks_tropical <- cstocks_tidy %>% 
  mutate(affinity = ifelse(abs(Latitude) <= 23, "tropical", ifelse(abs(Latitude) > 23 & Posi == "Posi", "Posi", "temperate"))) %>% 
  filter(!is.na(affinity)) %>% 
  select(Latitude, affinity, cstocks, depth, Posi) %>% 
  filter(!is.na(cstocks)) %>% 
  group_by(affinity, depth) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T))

# Barplot
cstocks_tropical %>% 
  ggplot(aes(x = affinity, y = mean_stock, ymin=mean_stock-std.error, ymax=mean_stock+std.error)) +
  geom_bar(aes(fill = affinity), stat = "identity") +
  scale_fill_manual(values = c("#1f78b4", "#31a354", "#b2df8a")) +
  geom_errorbar(width=.3, position = position_dodge(preserve = "single", width = 0.9)) +
  geom_text(aes(y = mean_stock+std.error + 20, label = str_c("(", n, ")"))) +
  xlab("") +
  scale_x_discrete(labels = c(expression(italic("Posidonia oceanica")), "Temperate", "Tropical")) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave("Figs/Cstocks_tropicalVStemperate_by_depthBARPLOT.pdf")


# Boxplot
cstocks_tidy %>% 
  mutate(affinity = ifelse(abs(Latitude) <= 23, "tropical", ifelse(abs(Latitude) > 23 & Posi == "Posi", "Posi", "temperate"))) %>% 
  filter(!is.na(affinity)) %>% 
  ggplot() +
  geom_boxplot(aes(x = affinity, y = cstocks, fill = affinity, colour = affinity)) +
  scale_colour_manual(values = c("#1f78b4", "#31a354", "#b2df8a")) +
  scale_fill_manual(values = c("#1f78b4", "#31a354", "#b2df8a")) +
  geom_point(data = cstocks_tropical, aes(x = affinity, y = median_stock), shape = 15, size = 2) +
  geom_text(data = cstocks_tropical, aes(x = affinity, y = max_stock + 40, label = str_c("(", n, ")"))) +
  xlab("") +
  scale_x_discrete(labels = c(expression(italic("Posidonia oceanica")), "Temperate", "Tropical")) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave("Figs/Cstocks_tropicalVStemperate_by_depthBOXPLOT.pdf")






