# New PLSR


# # # # # # # # #
# LIBRARIES ----
# # # # # # # # #

library(tidyverse)
library(tidylog)
library(pls)
library(caret)
library(matrixStats)
library(ggsci)


# # # # # # # # # # # # # # # # # # # # # # # # 
# Loading TRAITS data set and checking NAs ----
# # # # # # # # # # # # # # # # # # # # # # # #

# Loading traits data
traits <- read_csv(file = "Seagrass_traits_nov2020.csv")
traits2 <- read_csv(file = "Seagrass_traits_nov2020+d13soil.csv")
traits == traits2[,-(38:40)] # OK just to check that both data sets are the same. And they are. 
traits <- traits2 # We choose traits2, because it's the one with d13Csoil data.
rm(list = "traits2")

# According to what we talked with H. Kennedy, we complete surf.d13C for Z.noltii with the value from Z.marina.
traits$surf.d13C[which(traits$Species == "Zostera noltii")] <- -19.7

# We'll do it with all species that have d13C but Thalassodendron, in this way we get the most traits for most species.
traits_full <- traits %>% 
  # select(-surf_20cm.d13C, -cm20_50.d13C) %>% 
  filter(Species != "Thalassodendron ciliatum") %>% 
  # filter(Species != "Posidonia oceanica") %>%
  select(where(~!any(is.na(.)))) # This line drops any column containing a missing value

# To be used as filtering vector in the next section 
SpeciesList <- traits_full$Species


# Loading carbon stock data
source("01_DataImport&Corrections_CarbonReview.R")
rm(list = c("cstocks_tidy", "cstocks")) # Because we want to use, only, the 20 cm stocks data, not carbon densities.

cstocks_data <- cstocks_tidy20Stocks %>% 
  filter(Meadow_type == "monospecific") %>%
  mutate(Lat_Zone_Posi = ifelse(Posi == "Posi", "Posi", Lat_Zone),
         Lat_Zone = factor(Lat_Zone),
         CoreID_unique = factor(CoreID_unique),
         Geomorph = factor(recode(FIN_TYP,
                                  `0` = "Endorheic or Glaciated",
                                  `1` = "Small deltas",
                                  `2` = "Tidal systems",
                                  `3` = "Lagoons",
                                  `4` = "Fjords and fjaerds",
                                  `5` = "Large rivers",
                                  `6` = "Karst",
                                  `7` = "Arheic"))) %>% 
  select(CoreID_unique, Latitude, Species, Geomorph, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$Geomorph <- droplevels(cstocks_data$Geomorph)

# The full data set with traits + "raw" cstocks data + other contextual predictors (e.g Fintype, depth...)
cstocks_traits <- cstocks_data %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Species %in% SpeciesList) %>%
  left_join(traits_full, by = "Species") %>%
  mutate(diff_LvsSurfd13C = surf.d13C - L.d13C) %>% 
  select(Species,
         cstocks,
         depth,
         Geomorph,
         Lat_Zone,
         starts_with("L."),
         "ABG.biomass",
         "BG.biomass",
         "Ratio.AGBGbiomass",
         starts_with("R."),
         "N.leaves.shoot",
         diff_LvsSurfd13C,
         -L.d13C)



# # # # # #
# PLSR ----
# # # # # #
# set.seed(123)
training.samples <- createDataPartition(cstocks_traits$cstocks, p = 0.8, list = FALSE)
train.data  <- cstocks_traits[training.samples, ]
test.data <- cstocks_traits[-training.samples, ]

fitControl <- trainControl(method = "cv",  number = 1000)

set.seed(123)
model.initial <- train(cstocks ~ depth + Geomorph + Lat_Zone + L.width + L.length +
                         L.fibre + L.breakforce + L.lifespan + L.carbon + L.nitrogen + L.phosphorus + 
                         L.ratioNP + L.ratioCN + L.ratioCP + L.lignin + L.product.rate + L.plastochrone.interv + 
                         L.turnover + ABG.biomass + BG.biomass + Ratio.AGBGbiomass + R.diameter + 
                         R.internode.length + R.helong.rate + R.plastochrone.interv + N.leaves.shoot + diff_LvsSurfd13C,
                       data = train.data, 
                       method = "simpls",
                       scale = TRUE,
                       tuneLength = 20)

model.initial$bestTune
model.initial$finalModel

plot(model.initial)
summary(model.initial)

plot(varImp(model.initial, scale = F))
