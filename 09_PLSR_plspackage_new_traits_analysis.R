# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/152-principal-component-and-partial-least-squares-regression-essentials/#partial-least-squares-regression
# # # # # # #
# Analysing seagrass trait data set provided by Hilary Kennedy
# Jordi F. Pagès
# 23-02-2021
# University of Barcelona 
# (currently working from home/CEAB)
# # # # # # 

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

traits <- read_csv(file = "Seagrass_traits_nov2020.csv")

# We'll doing with all species but Thalassodendron, in this way we get the most traits for most species.
traits_full <- traits %>% 
  filter(Species != "Thalassodendron ciliatum") %>% 
  # filter(Species != "Posidonia oceanica") %>% 
  select(where(~!any(is.na(.)))) # This line drops any column containing a missing value

# To be used as filtering vector in the next section 
SpeciesList <- traits_full$Species


# # # # # # # # # # # # # # # # # # # # # # 
#    Loading carbon stocks data set    ----
#     to get the response variable        # 
#  that we aim at predicting with PLS     #
# # # # # # # # # # # # # # # # # # # # # #

source("1_DataImport&Corrections_CarbonReview.R")

cstocks_data <- cstocks_tidy %>% 
  filter(Meadow_type == "monospecific") %>%
  mutate(Lat_Zone_Posi = ifelse(Posi == "Posi", "Posi", Lat_Zone),
         FIN_TYP_names = factor(recode(FIN_TYP,
                                       `0` = "Endorheic or Glaciated",
                                       `1` = "Small deltas",
                                       `2` = "Tidal systems",
                                       `3` = "Lagoons",
                                       `4` = "Fjords and fjaerds",
                                       `5` = "Large rivers",
                                       `6` = "Karst",
                                       `7` = "Arheic"))) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)

# The full data set with traits + "raw" cstocks data + other contextual predictors (e.g Fintype, depth...)
cstocks_traits <- cstocks_data %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Species %in% SpeciesList) %>%
  left_join(traits_full, by = "Species") %>%
  select(Species,
         cstocks,
         depth,
         FIN_TYP_names,
         Lat_Zone,
         starts_with("L."),
         "ABG.biomass",
         "BG.biomass",
         "Ratio.AGBGbiomass",
         starts_with("R."),
         "N.leaves.shoot")

# The data set only with traits + and just one observation per species (mean, or median) + no more predictors 
# (because we can't summarise categorical variables)
# Uncomment to run PLSR with only one observation of cstocks per species (thus not repeating trait data for every cstock observation)
# cstocks_traits <- cstocks_tidy %>% 
#   filter(!is.na(cstocks)) %>% 
#   group_by(Species) %>%
#   summarise(cstocks.mean = mean(cstocks),
#             # std.error = std.error(cstocks),
#             cstocks.median = median(cstocks),
#             cstocks.n = n()) %>%
#   filter(Species %in% SpeciesList) %>%
#   left_join(traits_full, by = "Species") %>%
#   select(Species,
#          cstocks.mean,
#          starts_with("L."),
#          "ABG.biomass",
#          "BG.biomass",
#          "Ratio.AGBGbiomass",
#          starts_with("R."),
#          "N.leaves.shoot")
  


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PLS with traits_full (without Thalassodendron) NO LOOP ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# We’ll randomly split the data into training set (80% for building a predictive model) 
# and test set (20% for evaluating the model). Make sure to set seed for reproducibility.
# Split the data into training and test set
# set.seed(123)
cstocks_traits <- as.data.frame(cstocks_traits)
training.samples <- cstocks_traits$cstocks %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- cstocks_traits[as.vector(training.samples), ]
test.data <- cstocks_traits[-as.vector(training.samples), ]

# The R function train() [caret package] provides an easy workflow to compute PLS by invoking the pls package.
# It has an option named method, which can take the value pls.
# An additional argument is scale = TRUE for standardizing the variables to make them comparable.
# caret uses cross-validation to automatically identify the optimal number of principal components (ncomp) to be incorporated in the model.
# Here, we’ll test 10 different values of the tuning parameter ncomp. This is specified using the option tuneLength. 
# The optimal number of principal components is selected so that the cross-validation error (RMSE) is minimized.

# Build the model on training set
# set.seed(123)
model.initial <- train(cstocks ~ depth + FIN_TYP_names + Lat_Zone + L.width + L.length +
                         L.fibre + L.breakforce + L.lifespan + L.carbon + L.nitrogen + L.phosphorus + 
                         L.ratioNP + L.ratioCN + L.ratioCP + L.lignin + L.product.rate + L.plastochrone.interv + 
                         L.turnover + L.d13C + ABG.biomass + BG.biomass + Ratio.AGBGbiomass + R.diameter + 
                         R.internode.length + R.helong.rate + R.plastochrone.interv + N.leaves.shoot,
                       data = train.data, 
                       method = "simpls",
                       scale = TRUE,
                       trControl = trainControl("cv"),
                       tuneLength = 12)

# Plot model RMSE vs different values of components
plot(model.initial, metric = "RMSE")
plot(model.initial, metric = "Rsquared")

# Print the best tuning parameter ncomp that minimize the cross-validation error, RMSE
model.initial$bestTune

# Summarize the final model
summary(model.initial$finalModel) # With only 2 components the Rsquared of cstocks is quite high, and it flattens after 2 components. 
                                  # In contrast, the % variance explained goes up and up for the predictors til the 8-9th component. 

# Model performance metrics
r2plot <- caret::R2(model.initial %>% predict(test.data), test.data$cstocks)
# R2 = 0.67 very good!

# Make predictions an plot them in 1:1 plots to see residual variation
test.data$predictions <- model.initial %>% 
  predict(test.data)
p <- ggplot(data = test.data, mapping = aes(x = predictions, y = cstocks)) +
  geom_point(aes(colour = Species)) +
  scale_colour_d3("category20", 
                  alpha = 1,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  geom_abline() +
  coord_cartesian(xlim = c(0,10), ylim = c(0,10)) +
  ylab(bquote('Observed C density (Mg C' ~ha^-1~cm^-1* ')')) +
  xlab("Predicted C density (PLSR)") + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p + annotate(geom = "text", x = 7.5, y = 10, label = bquote(R^2 == .(round(r2plot, 2))))
# ggsave(filename = "Figs/PLSRnotLooped/cstocksVSpredicted_not_loopR2=0.52OK.pdf")
# Predictions are not super good, but overall the species effect is very obvious. With some intraspecific residual variance
# that it's escaping our model (we're not managing to predict this intraspecific variance).

# Variable Importance on the Projection
ggplot(varImp(model.initial))
# ggsave(filename = "Figs/varImportance_not_loop.pdf")
varIMP <- varImp(model.initial)

# Model results
model.coefficients <- as_tibble(model.initial$finalModel$coefficients)
model.coefficients.tibble <- tibble(predictors = rownames(model.initial$finalModel$coefficients), model.coefficients)

All_scoresPlsr <- model.initial$finalModel$scores
scoresPlsr1 <- All_scoresPlsr[1:432,1]

train.data$scoresPlsr1 <- All_scoresPlsr[1:432,1]
ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = cstocks)) +
  geom_point(aes(colour = Species)) +
  geom_smooth(method = "lm", se=T, colour = "black")
# ggsave(filename = "Figs/cstocksVSscoresPLSR1(notLoop).pdf")
 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PLS with traits_full (without Thalassodendron) LOOP!!! ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

bestTune <- NULL
bestTune.out <- NULL
r2 <- NULL
r2.out <- NULL
varImp <- NULL
imp.out <- NULL
coef1.out <- NULL
coef2.out <- NULL
cstocks_traits <- as.data.frame(cstocks_traits)
# scoresPlsr1.out <- NULL
for(i in 1:100){
  # We’ll randomly split the data into training set (80% for building a predictive model) 
  # and test set (20% for evaluating the model). Make sure to set seed for reproducibility.
  # Split the data into training and test set
  # set.seed(123)
  training.samples <- cstocks_traits$cstocks %>%
    createDataPartition(p = 0.8, list = FALSE)
  train.data  <- cstocks_traits[as.vector(training.samples), ]
  test.data <- cstocks_traits[-as.vector(training.samples), ]
  
  # The R function train() [caret package] provides an easy workflow to compute PCR and PLS by invoking the pls package.
  # It has an option named method, which can take the value pcr or pls.
  # An additional argument is scale = TRUE for standardizing the variables to make them comparable.
  # caret uses cross-validation to automatically identify the optimal number of principal components (ncomp) to be incorporated in the model.
  # Here, we’ll test 10 different values of the tuning parameter ncomp. 
  # This is specified using the option tuneLength. 
  # The optimal number of principal components is selected so that the cross-validation error (RMSE) is minimized.
  
  # Build the model on training set
  # set.seed(123)
  model.initial <- train(cstocks ~ depth + FIN_TYP_names + Lat_Zone + L.width + L.length +
                           L.fibre + L.breakforce + L.lifespan + L.carbon + L.nitrogen + L.phosphorus + 
                           L.ratioNP + L.ratioCN + L.ratioCP + L.lignin + L.product.rate + L.plastochrone.interv + 
                           L.turnover + L.d13C + ABG.biomass + BG.biomass + Ratio.AGBGbiomass + R.diameter + 
                           R.internode.length + R.helong.rate + R.plastochrone.interv + N.leaves.shoot,
                         data = train.data, 
                         method = "simpls",
                         scale = TRUE,
                         trControl = trainControl("cv"),
                         tuneLength = 12)
  
  # Print the best tuning parameter ncomp that minimizes the cross-validation error, RMSE
  bestTune <- model.initial$bestTune
  
  # Summarize the final model
  # summary(model.initial$finalModel)
  
  # Model coefficients
  model.coefficients <- as_tibble(model.initial$finalModel$coefficients)
  model.coefficients.tibble <- tibble(predictors = rownames(model.initial$finalModel$coefficients), model.coefficients)
  coef1 <- model.coefficients.tibble$`.outcome.1 comps`
  if(bestTune > 1){
    coef2 <- model.coefficients.tibble$`.outcome.2 comps`
  }
  
  # Variable Importance on the Projection
  varIMP <- varImp(model.initial)
  
  # Model performance metrics
  predictions <- model.initial %>% predict(test.data)
  r2 <- caret::R2(predictions, test.data$cstocks)
  #  Rsquare
  #  0.45
  
  if(i == 1){
    imp.out <- varIMP$importance
    coef1.out <- coef1
    if(bestTune > 1){
      coef2.out <- coef2
    }
  }
  else{
    imp.out <- cbind(imp.out, varIMP$importance)
    coef1.out <- cbind(coef1.out, coef1)
    if(bestTune > 1){
      coef2.out <- cbind(coef2.out, coef2)
    }
  }
  r2.out <- c(r2.out, r2)
  bestTune.out <- c(bestTune.out, bestTune[[1]])
  if(i<20){
    train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:432,1]
    ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = cstocks)) +
      geom_point(aes(colour = Species)) +
      geom_smooth(method = "lm", se=T, colour = "black") + 
      theme_bw() +
      theme(text = element_text(size=14),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
     file_name <- paste("Figs/PLSRLooped/SpeciesVSplsr1_i", i, ".pdf", sep = "")
    ggsave(filename = file_name)
  }
}


# Variable importance
names(imp.out) <- 1:100
# imp.out <- as.data.frame(t(imp.out))
# imp.means <- sort(colMeans(imp.out), decreasing = T)
# imp.sds <- colSds(as.matrix(imp.out))
imp.means <- rowMeans(imp.out)
imp.sds <- rowSds(as.matrix(imp.out))
imp.summarised <- as.data.frame(cbind(imp.means, imp.sds))
imp.summarised$x <- rownames(imp.summarised)


# Coefficients 1st component
names(coef1.out) <- 1:100
# coef1.out <- as.data.frame(t(coef1.out))
# coef1.means <- colMeans(coef1.out)
# coef1.sds <- colSds(as.matrix(coef1.out))
coef1.means <- rowMeans(coef1.out)
coef1.sds <- rowSds(as.matrix(coef1.out))
coef1.summarised <- as.data.frame(cbind(coef1.means, coef1.sds))
coef1.summarised$x <- rownames(imp.summarised)

# Coefficients 2nd component
# coef2.out <- as.data.frame(t(coef2.out))
# rownames(coef2.out) <- 1:length(rownames(coef2.out))
names(coef2.out) <- 1:100
# coef2.means <- colMeans(coef2.out)
# coef2.sds <- colSds(as.matrix(coef2.out))
coef2.means <- rowMeans(coef2.out)
coef2.sds <- rowSds(as.matrix(coef2.out))
coef2.summarised <- as.data.frame(cbind(coef2.means, coef2.sds))
coef2.summarised$x <- rownames(imp.summarised)

# Best tune
mean(bestTune.out) # 9.27 components
sd(bestTune.out) # 0.94 components
min(bestTune.out) # 7 components
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(bestTune.out) # 9 components is the most frequent one

bestTune.out <- as.data.frame(bestTune.out)

max(bestTune.out)

bestTune.out %>% 
  mutate(name = as.factor(bestTune.out)) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = name, y = n)) +
  geom_bar(stat = 'identity')

bestTune.out %>% 
  mutate(name = as.factor(bestTune.out)) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>%
  mutate(percent = (n/sum(n))*100) %>% 
  ggplot(aes(x = name, y = percent)) +
  geom_bar(stat = 'identity') +
  ylim(c(0,100)) + 
  ylab("%") +
  xlab('Number of components') + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave('Figs/PLSRLooped/NumberOfComponents.pdf')

# R2
mean(r2.out) # R2 = 0.4270163
sd(r2.out)
min(r2.out)
max(r2.out)

# # Scores PLSR1
# names(scoresPlsr1.out) <- 1:100
# scoresPlsr1.out <- as.data.frame(scoresPlsr1.out)
# scoresPlsr1.means <- rowMeans(scoresPlsr1.out)
# scoresPlsr1.sds <- rowSds(as.matrix(scoresPlsr1.out))
# scoresPlsr1.summarised <- as.data.frame(cbind(scoresPlsr1.means, scoresPlsr1.sds))


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# GRAPHS FOR THE PAPER HIGHLIGHTING MODEL RESULTS ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Variable importance on the projection
imp.summarised <- imp.summarised %>% 
  mutate(fill = as.numeric(imp.means < 50))

ggplot(imp.summarised, aes(x=reorder(x, imp.means, FUN = max), y = imp.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=imp.means, ymax=imp.means+imp.sds), width=.2, position=position_dodge(.9)) + 
  # scale_y_continuous(name = "Variable importance", limits = c(0,25,50,75,100)) +
  xlab(label = "") + 
  ylab(label ="Variable Importance in Projection (%)") +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/PLSRLooped/variableImportance(n=100).pdf")

# PLSR Coefficients 1st component
coef1.summarised %>% 
  # mutate(fill = ifelse(test = x %in% c("N.leaves.shoot", "BG.biomass", "L.lifespan", "L.d13C", "ABG.biomass", "L.lignin", "L.carbon"), 
  #                      yes = 1, no = 2)) %>% 
  left_join(imp.summarised, by = "x") %>% 
  ggplot(aes(x = reorder(x, coef1.means, FUN = max), y = coef1.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=coef1.means - coef1.sds, ymax = coef1.means + coef1.sds), width=.2, position=position_dodge(0.9)) + 
  ylab("PLSR coefficients 1st component") +
  xlab("") + 
  # scale_y_continuous(limits = c(-2, 2)) +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/PLSRLooped/PLSRcoefficients1stcomp(n=100).pdf")

# PLSR Coefficients 2nd component
coef2.summarised %>% 
  # mutate(fill = ifelse(test = x %in% c("N.leaves.shoot", "BG.biomass", "L.lifespan", "L.d13C", "ABG.biomass", "L.lignin", "L.carbon"), 
  #                      yes = 1, no = 2)) %>% 
  left_join(imp.summarised, by = "x") %>% 
  ggplot(aes(x = reorder(x, coef2.means, FUN = max), y = coef2.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=coef2.means - coef2.sds, ymax = coef2.means + coef2.sds), width=.2, position=position_dodge(0.9)) + 
  ylab("PLSR coefficients 2nd component") +
  xlab("") + 
  # scale_y_continuous(limits = c(-5, 5)) +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/PLSRLooped/PLSRcoefficients2ndcomp(n=100).pdf")


# Figures PLSR scores vs. RESPONSE VARIABLE (similar to Carrascal et al 2009)
# We use just the last result from the for loop as an example
train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:432,1]
ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = cstocks)) +
  geom_point(aes(colour = Species)) +
  # scale_colour_simpsons(alpha = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  geom_smooth(method = "lm", se=T, colour = "black") + 
  xlab("PLSR scores 1st component") + 
  ylab(bquote('Observed C density (Mg C' ~ha^-1~cm^-1* ')')) +
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
file_name <- paste("Figs/PLSRLooped/SpeciesVSplsr1_i", i, ".pdf", sep = "")
# ggsave(filename = file_name)


# Final model summary results
summary(model.initial)
# Data: 	X dimension: 432 33 
# Y dimension: 432 1
# Fit method: simpls
# Number of components considered: 10
# TRAINING: % variance explained
#             1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps  8 comps  9 comps  10 comps
# X           36.28    45.97    51.09    56.75    62.17    71.99    75.39    78.89    81.04     82.99
# .outcome    31.15    46.20    47.61    47.91    48.15    48.34    48.78    48.86    48.92     48.98






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# REPEATING ANALYSIS FOR TRAITS data set WITH d13 in soil and checking NAs ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Loading traits data
traits <- read_csv(file = "Seagrass_traits_nov2020.csv")
traits2 <- read_csv(file = "Seagrass_traits_nov2020+d13soil.csv")
traits == traits2[,-(38:40)] # OK just to check that both data sets are the same. And they are. 
traits <- traits2

# We'll doing with all species that have d13C but Thalassodendron, in this way we get the most traits for most species.
traits_full <- traits %>% 
  # select(-surf_20cm.d13C, -cm20_50.d13C) %>% 
  filter(Species != "Thalassodendron ciliatum") %>% 
  filter(Species != "Halodule wrightii") %>% 
  filter(Species != "Thalassia testudinum") %>% 
  filter(Species != "Zostera noltii") %>% 
  # filter(Species != "Posidonia oceanica") %>% 
  select(where(~!any(is.na(.)))) # This line drops any column containing a missing value

# To be used as filtering vector in the next section 
SpeciesList <- traits_full$Species


# Loading carbon stock data
source("1_DataImport&Corrections_CarbonReview.R")
cstocks_data <- cstocks_tidy %>% 
  filter(Meadow_type == "monospecific") %>%
  mutate(Lat_Zone_Posi = ifelse(Posi == "Posi", "Posi", Lat_Zone),
         FIN_TYP_names = factor(recode(FIN_TYP,
                                       `0` = "Endorheic or Glaciated",
                                       `1` = "Small deltas",
                                       `2` = "Tidal systems",
                                       `3` = "Lagoons",
                                       `4` = "Fjords and fjaerds",
                                       `5` = "Large rivers",
                                       `6` = "Karst",
                                       `7` = "Arheic"))) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)

# The full data set with traits + "raw" cstocks data + other contextual predictors (e.g Fintype, depth...)
cstocks_traits <- cstocks_data %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Species %in% SpeciesList) %>%
  left_join(traits_full, by = "Species") %>%
  mutate(diff_Lvs20cmd13C = surf_20cm.d13C - L.d13C) %>% 
  select(Species,
         cstocks,
         depth,
         FIN_TYP_names,
         Lat_Zone,
         starts_with("L."),
         "ABG.biomass",
         "BG.biomass",
         "Ratio.AGBGbiomass",
         starts_with("R."),
         "N.leaves.shoot",
         "surf_20cm.d13C",
         -L.d13C)


# Looped version with d13C ----
bestTune <- NULL
bestTune.out <- NULL
r2 <- NULL
r2.out <- NULL
varImp <- NULL
imp.out <- NULL
coef1.out <- NULL
coef2.out <- NULL
cstocks_traits <- as.data.frame(cstocks_traits)
# scoresPlsr1.out <- NULL
for(i in 1:100){
  # We’ll randomly split the data into training set (80% for building a predictive model) 
  # and test set (20% for evaluating the model). Make sure to set seed for reproducibility.
  # Split the data into training and test set
  # set.seed(123)
  training.samples <- cstocks_traits$cstocks %>%
    createDataPartition(p = 0.8, list = FALSE)
  train.data  <- cstocks_traits[as.vector(training.samples), ]
  test.data <- cstocks_traits[-as.vector(training.samples), ]
  
  # The R function train() [caret package] provides an easy workflow to compute PCR and PLS by invoking the pls package.
  # It has an option named method, which can take the value pcr or pls.
  # An additional argument is scale = TRUE for standardizing the variables to make them comparable.
  # caret uses cross-validation to automatically identify the optimal number of principal components (ncomp) to be incorporated in the model.
  # Here, we’ll test 10 different values of the tuning parameter ncomp. 
  # This is specified using the option tuneLength. 
  # The optimal number of principal components is selected so that the cross-validation error (RMSE) is minimized.
  
  # Build the model on training set
  # set.seed(123)
  model.initial <- train(cstocks ~ depth + FIN_TYP_names + Lat_Zone + L.width + L.length +
                           L.fibre + L.breakforce + L.lifespan + L.carbon + L.nitrogen + L.phosphorus + 
                           L.ratioNP + L.ratioCN + L.ratioCP + L.lignin + L.product.rate + L.plastochrone.interv + 
                           L.turnover + ABG.biomass + BG.biomass + Ratio.AGBGbiomass + R.diameter + 
                           R.internode.length + R.helong.rate + R.plastochrone.interv + N.leaves.shoot + surf_20cm.d13C,
                         data = train.data, 
                         method = "simpls",
                         scale = TRUE,
                         trControl = trainControl("cv"),
                         tuneLength = 12)
  
  # Print the best tuning parameter ncomp that minimizes the cross-validation error, RMSE
  bestTune <- model.initial$bestTune
  
  # Summarize the final model
  # summary(model.initial$finalModel)
  
  # Model coefficients
  model.coefficients <- as_tibble(model.initial$finalModel$coefficients)
  model.coefficients.tibble <- tibble(predictors = rownames(model.initial$finalModel$coefficients), model.coefficients)
  coef1 <- model.coefficients.tibble$`.outcome.1 comps`
  if(bestTune > 1){
    coef2 <- model.coefficients.tibble$`.outcome.2 comps`
  }
  
  # Variable Importance on the Projection
  varIMP <- varImp(model.initial)
  
  # Model performance metrics
  predictions <- model.initial %>% predict(test.data)
  r2 <- caret::R2(predictions, test.data$cstocks)
  #  Rsquare
  #  0.45
  
  if(i == 1){
    imp.out <- varIMP$importance
    coef1.out <- coef1
    if(bestTune > 1){
      coef2.out <- coef2
    }
  }
  else{
    imp.out <- cbind(imp.out, varIMP$importance)
    coef1.out <- cbind(coef1.out, coef1)
    if(bestTune > 1){
      coef2.out <- cbind(coef2.out, coef2)
    }
  }
  r2.out <- c(r2.out, r2)
  bestTune.out <- c(bestTune.out, bestTune[[1]])
  # if(i<20){
  #   train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:392,1]
  #   ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = cstocks)) +
  #     geom_point(aes(colour = Species)) +
  #     geom_smooth(method = "lm", se=T, colour = "black") + 
  #     theme_bw() +
  #     theme(text = element_text(size=14),
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank())
  #   file_name <- paste("Figs/PLSRLooped/SpeciesVSplsr1_i", i, ".pdf", sep = "")
  #   ggsave(filename = file_name)
  # }
}


# Variable importance
names(imp.out) <- 1:100
# imp.out <- as.data.frame(t(imp.out))
# imp.means <- sort(colMeans(imp.out), decreasing = T)
# imp.sds <- colSds(as.matrix(imp.out))
imp.means <- rowMeans(imp.out)
imp.sds <- rowSds(as.matrix(imp.out))
imp.summarised <- as.data.frame(cbind(imp.means, imp.sds))
imp.summarised$x <- rownames(imp.summarised)


# Coefficients 1st component
names(coef1.out) <- 1:100
# coef1.out <- as.data.frame(t(coef1.out))
# coef1.means <- colMeans(coef1.out)
# coef1.sds <- colSds(as.matrix(coef1.out))
coef1.means <- rowMeans(coef1.out)
coef1.sds <- rowSds(as.matrix(coef1.out))
coef1.summarised <- as.data.frame(cbind(coef1.means, coef1.sds))
coef1.summarised$x <- rownames(imp.summarised)

# Coefficients 2nd component
# coef2.out <- as.data.frame(t(coef2.out))
# rownames(coef2.out) <- 1:length(rownames(coef2.out))
names(coef2.out) <- 1:100
# coef2.means <- colMeans(coef2.out)
# coef2.sds <- colSds(as.matrix(coef2.out))
coef2.means <- rowMeans(coef2.out)
coef2.sds <- rowSds(as.matrix(coef2.out))
coef2.summarised <- as.data.frame(cbind(coef2.means, coef2.sds))
coef2.summarised$x <- rownames(imp.summarised)

# Best tune
mean(bestTune.out) # 9.27 components
sd(bestTune.out) # 0.94 components
min(bestTune.out) # 7 components
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(bestTune.out) # 9 components is the most frequent one

bestTune.out <- as.data.frame(bestTune.out)

max(bestTune.out)

bestTune.out %>% 
  mutate(name = as.factor(bestTune.out)) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = name, y = n)) +
  geom_bar(stat = 'identity')

bestTune.out %>% 
  mutate(name = as.factor(bestTune.out)) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>%
  mutate(percent = (n/sum(n))*100) %>% 
  ggplot(aes(x = name, y = percent)) +
  geom_bar(stat = 'identity') +
  ylim(c(0,100)) + 
  ylab("%") +
  xlab('Number of components') + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave('Figs/PLSRLooped/NumberOfComponents.pdf')

# R2
mean(r2.out) # R2 = 0.4270163
sd(r2.out)
min(r2.out)
max(r2.out)

# # Scores PLSR1
# names(scoresPlsr1.out) <- 1:100
# scoresPlsr1.out <- as.data.frame(scoresPlsr1.out)
# scoresPlsr1.means <- rowMeans(scoresPlsr1.out)
# scoresPlsr1.sds <- rowSds(as.matrix(scoresPlsr1.out))
# scoresPlsr1.summarised <- as.data.frame(cbind(scoresPlsr1.means, scoresPlsr1.sds))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# GRAPHS FOR THE PAPER HIGHLIGHTING MODEL RESULTS FOR SOIL D13C ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Variable importance on the projection
imp.summarised <- imp.summarised %>% 
  mutate(fill = as.numeric(imp.means < 50))

ggplot(imp.summarised, aes(x=reorder(x, imp.means, FUN = max), y = imp.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=imp.means, ymax=imp.means+imp.sds), width=.2, position=position_dodge(.9)) + 
  # scale_y_continuous(name = "Variable importance", limits = c(0,25,50,75,100)) +
  xlab(label = "") + 
  ylab(label ="Variable Importance in Projection (%)") +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/PLSRLooped_with_soil_d13C/Diff20cmLeaf_variableImportance(n=100).pdf")

# PLSR Coefficients 1st component
coef1.summarised %>% 
  # mutate(fill = ifelse(test = x %in% c("N.leaves.shoot", "BG.biomass", "L.lifespan", "L.d13C", "ABG.biomass", "L.lignin", "L.carbon"), 
  #                      yes = 1, no = 2)) %>% 
  left_join(imp.summarised, by = "x") %>% 
  ggplot(aes(x = reorder(x, coef1.means, FUN = max), y = coef1.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=coef1.means - coef1.sds, ymax = coef1.means + coef1.sds), width=.2, position=position_dodge(0.9)) + 
  ylab("PLSR coefficients 1st component") +
  xlab("") + 
  # scale_y_continuous(limits = c(-2, 2)) +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/PLSRLooped_with_soil_d13C/Diff20cmLeaf_PLSRcoefficients1stcomp(n=100).pdf")

# PLSR Coefficients 2nd component
coef2.summarised %>% 
  # mutate(fill = ifelse(test = x %in% c("N.leaves.shoot", "BG.biomass", "L.lifespan", "L.d13C", "ABG.biomass", "L.lignin", "L.carbon"), 
  #                      yes = 1, no = 2)) %>% 
  left_join(imp.summarised, by = "x") %>% 
  ggplot(aes(x = reorder(x, coef2.means, FUN = max), y = coef2.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=coef2.means - coef2.sds, ymax = coef2.means + coef2.sds), width=.2, position=position_dodge(0.9)) + 
  ylab("PLSR coefficients 2nd component") +
  xlab("") + 
  # scale_y_continuous(limits = c(-5, 5)) +
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/PLSRLooped_with_soil_d13C/Diff20cmLeaf_PLSRcoefficients2ndcomp(n=100).pdf")


# Figures PLSR scores vs. RESPONSE VARIABLE (similar to Carrascal et al 2009)
# We use just the last result from the for loop as an example
train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:392,1]
ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = cstocks)) +
  geom_point(aes(colour = Species)) +
  # scale_colour_simpsons(alpha = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  geom_smooth(method = "lm", se=T, colour = "black") + 
  xlab("PLSR scores 1st component") + 
  ylab(bquote('Observed C density (Mg C' ~ha^-1~cm^-1* ')')) +
  theme_bw() +
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
file_name <- paste("Figs/PLSRLooped_with_soil_d13C/Diff20cmLeaf_SpeciesVSplsr1_i", i, ".pdf", sep = "")
# ggsave(filename = file_name)


# Final model summary results
summary(model.initial)
# Data: 	X dimension: 392 33 
#         Y dimension: 392 1
# Fit method: simpls
# Number of components considered: 7
# TRAINING: % variance explained
#             1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps
# X           37.29    48.10    55.90    64.73    71.18    75.74    78.37
# .outcome    28.23    42.66    44.49    45.05    45.30    45.39    45.43



