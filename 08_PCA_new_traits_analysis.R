# # # # # #
# Analysing seagrass trait data set provided by Hilary Kennedy
# Jordi F. Pagès
# 19-02-2021
# University of Barcelona 
# (currently working from home/CEAB)
# # # # # # 


# # # # # # # # #
# LIBRARIES ----
# # # # # # # # #

library(tidyverse)
library(tidylog)
library(ggfortify)
library(ggsci)

# # # # # # # # # # # # # # # # # # # # # # # # 
# Loading TRAITS data set and checking NAs ----
# # # # # # # # # # # # # # # # # # # # # # # #

traits <- read_csv(file = "Seagrass_traits_nov2020.csv")

SpeciesNAs <- traits %>% 
  mutate(NAs = rowSums(is.na(.)), ) %>% 
  select(Species, NAs)

colNAs <- traits %>% 
  # filter(Species != "Thalassodendron ciliatum") %>%
  select(-Species) %>% 
  is.na %>% 
  colSums
traitsNAs <- tibble(traits = names(traits)[-1], colNAs)


# # # # # # # # # # # # # # # # # # #
# PCA with all TRAITS        ----
# (we have to delete species with NAs)
# # # # # # # # # # # # # # # # # # #

# Cymodocea rotundata, Halophila ovalis, Thalassia testudinum, Zostera marina, Zostera noltii 
# are the only species with data for all traits. Since according to previous figures, there aren't big C-density differences
# among these species, let's drop 1 trait, and at least we'll be able to include P.oceanica information too. 

# So, we drop trait called L.ndf, to be able to use P.oceanica trait data too
alltraits_df  <-  traits %>% 
  select(-L.ndf) %>% 
  mutate(NAs = rowSums(is.na(.)), ) %>% 
  filter(NAs == 0) %>% 
  select(-NAs)

alltraits_PCA <- prcomp(alltraits_df[,-1], scale. = TRUE)
autoplot(alltraits_PCA, data = alltraits_df, colour = "Species")
autoplot(alltraits_PCA, 
         data = alltraits_df, 
         colour = "Species", 
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.label.hjust = 1.1)
# ggsave(filename = "Figs/Final_Figures_March2021/PCA/PCA_alltraits_biplot.pdf")

# This plot, allows us to check which variables are correlated to which, and decide which ones to drop,
# to increase the sample of species available.

# For example, according to traitsNAs:
# S.plastochrone.interv	hasn't got data for 5 species. And according to the biplot, it's super correlated to L.lifespan, among others.
# Thus, S.plastochrone can be dropped.
# BG.turnover hasn't got data for 5 species. According to biplot, it's correlated to BG.production and ABG.turnover
# Thus, BG.turnover can be dropped.
# GPP --> 3 NAs. It's a bit trickier, because it's not super positively correlated with other variables. But it is negatively correlated with others.
# Then the group ABG.turnover, BG.production has 3 NAs, and is correlated with BG.turnover, with 5 NAs... 


# # # # # # # # # # # # # # # # # # #
# PCA with all SPECIES        ----
# (we have to drop traits with NAs)
# # # # # # # # # # # # # # # # # # #

traits_full <- traitsNAs %>% 
  filter(colNAs == 0) %>% 
  select(traits) %>% 
  print()

traits_missing <- traitsNAs %>% 
  filter(colNAs > 0) %>% 
  print()

allspecies_df <- as.data.frame(traits[,which(names(traits) %in% traits_full$traits)])
rownames(allspecies_df) <- traits$Species

allspecies_PCA <- prcomp(allspecies_df, scale. = TRUE)
autoplot(allspecies_PCA, data = allspecies_df)
autoplot(allspecies_PCA, 
         data = allspecies_df, 
         label = TRUE, loadings = TRUE,
         loadings.label = TRUE, 
         loadings.label.hjust = 1.1)
# ggsave(filename = "Figs/Final_Figures_March2021/PCA/PCA_allspecies_biplot.pdf")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PCA with all SPECIES BUT Thalassodendron ciliatum       ----
# (because this species has lots of missing trait data)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

colNAs_thalasso <- traits %>% 
  filter(Species != "Thalassodendron ciliatum") %>%
  select(-Species) %>% 
  is.na %>% 
  colSums

traitsNAs_thalasso <- tibble(traits = names(traits)[-1], colNAs_thalasso)

traits_full_thalasso <- traitsNAs_thalasso %>% 
  filter(colNAs_thalasso == 0) %>% 
  select(traits) %>% 
  print(n = Inf)

traits_missing_thalasso <- traitsNAs_thalasso %>% 
  filter(colNAs_thalasso > 0) %>% 
  print()

allspecies_df_thalasso <- as.data.frame(traits[, which(names(traits) %in% traits_full_thalasso$traits)])
rownames(allspecies_df_thalasso) <- traits$Species

allspecies_df_thalassoOK <- allspecies_df_thalasso %>% 
  mutate(Species = rownames(.)) %>% 
  filter(Species != "Thalassodendron ciliatum") #%>% 
  #select(-Species)
rownames(allspecies_df_thalassoOK) <- rownames(allspecies_df_thalasso)[-13] # the position of Thalassodendron ciliatum
  
# PCA for all species except Thalassodendron ciliatum and as many traits as possible
allspecies_PCA_thalasso <- prcomp(allspecies_df_thalassoOK[,-26], scale. = TRUE)
autoplot(allspecies_PCA_thalasso, data = allspecies_df_thalassoOK)
autoplot(allspecies_PCA_thalasso, 
         data = allspecies_df_thalassoOK, 
         x = 1, y = 2,
         colour = "Species",
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.label.hjust = 1.1,
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         size = 4) +
  scale_colour_d3("category20", guide = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(-0.5,0.6)) + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave(filename = "Figs/Final_Figures_March2021/PCA/PCA_allspecies_but_Thalasso_biplot3.pdf")



