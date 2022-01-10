# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/152-principal-component-and-partial-least-squares-regression-essentials/#partial-least-squares-regression
# # # # # # #
# Analysing seagrass trait data set provided by Hilary Kennedy
# Jordi F. Pagès
# 09-04-2021
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


# Looped version with diff_LvsSurfd13C ----
predictor.var.out <- NULL
outcome.var.out <- NULL
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
  training.samples <- createDataPartition(cstocks_traits$cstocks, p = 0.8, list = FALSE)
  train.data  <- cstocks_traits[training.samples, ]
  test.data <- cstocks_traits[-training.samples, ]
  
  # The R function train() [caret package] provides an easy workflow to compute PCR and PLS by invoking the pls package.
  # It has an option named method, which can take the value pcr or pls.
  # An additional argument is scale = TRUE for standardizing the variables to make them comparable.
  # caret uses cross-validation to automatically identify the optimal number of principal components (ncomp) to be incorporated in the model.
  # Here, we’ll test 10 different values of the tuning parameter ncomp. 
  # This is specified using the option tuneLength. 
  # The optimal number of principal components is selected so that the cross-validation error (RMSE) is minimized.
  
  # Build the model on training set
  # set.seed(123)
  fitControl <- trainControl(method = "cv",  number = 10)
  model.initial <- train(cstocks ~ Geomorph + Lat_Zone + L.width + L.length +
                           L.fibre + L.breakforce + L.lifespan + L.carbon + L.nitrogen + L.phosphorus + 
                           L.ratioNP + L.ratioCN + L.ratioCP + L.lignin + L.product.rate + L.plastochrone.interv + 
                           L.turnover + ABG.biomass + BG.biomass + Ratio.AGBGbiomass + R.diameter + 
                           R.internode.length + R.helong.rate + R.plastochrone.interv + N.leaves.shoot + diff_LvsSurfd13C,
                         data = train.data, 
                         method = "simpls",
                         scale = TRUE,
                         trControl = fitControl,
                         tuneLength = 27)
  
  # Print the best tuning parameter ncomp that minimizes the cross-validation error, RMSE
  bestTune <- model.initial$bestTune
  
  # Summarize the final model
  # summary(model.initial$finalModel)
  predictor.var <- str_extract(capture.output(summary(model.initial))[7], "[0-9]+")
  outcome.var <- str_extract(capture.output(summary(model.initial))[8], "[0-9]+")
  
  model.coefficients <- as_tibble(model.initial$finalModel$coefficients)
  model.coefficients.tibble <- tibble(predictors = rownames(model.initial$finalModel$coefficients), model.coefficients)
  coef1 <- model.coefficients.tibble$`.outcome.1 comps`
  
  # Model loadings
  model.loadings <- model.initial$finalModel$loadings
  model.loadings1.tibble <- tibble(predictors = rownames(model.initial$finalModel$loadings), model.initial$finalModel$loadings[,1])
  loadings1 <- as.numeric(model.loadings1.tibble$`model.initial$finalModel$loadings[, 1]`)
  # if(bestTune > 1){
  #   coef2 <- model.coefficients.tibble$`.outcome.2 comps`
  # }
  
  # Variable Importance on the Projection
  varIMP <- varImp(model.initial, scale = FALSE)
  
  # Model performance metrics
  predictions <- model.initial %>% predict(test.data)
  r2 <- caret::R2(predictions, test.data$cstocks)
  #  Rsquare
  #  0.45
  
  # # Make predictions an plot them in 1:1 plots to see residual variation
  # test.data$predictions <- model.initial %>% 
  #   predict(test.data)
  # p <- ggplot(data = test.data, mapping = aes(x = predictions, y = cstocks)) +
  #   geom_point(aes(colour = Species)) +
  #   scale_colour_d3("category20", 
  #                   alpha = 1,
  #                   guide = guide_legend(override.aes = list(size = 3,
  #                                                            alpha = 1))) +
  #   geom_abline() +
  #   # coord_cartesian(xlim = c(0,10), ylim = c(0,10)) +
  #   ylab(bquote('Observed C density (Mg C' ~ha^-1~cm^-1* ')')) +
  #   xlab("Predicted C density (PLSR)") + 
  #   theme_bw() +
  #   theme(text = element_text(size=14),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank())
  # p + annotate(geom = "text", x = 25, y = 160, label = bquote(R^2 == .(round(r2, 2))))
  # ggsave(filename = "Figs/PLSRnotLooped/cstocksVSpredicted_not_loopR2=0.52OK.pdf")
  # Predictions are not super good, but overall the species effect is very obvious. With some intraspecific residual variance
  # that it's escaping our model (we're not managing to predict this intraspecific variance).
  
  if(i == 1){
    imp.out <- varIMP$importance
    loadings1.out <- loadings1
    # if(bestTune > 1){
    #   coef2.out <- coef2
    # }
  }
  else{
    imp.out <- cbind(imp.out, varIMP$importance)
    loadings1.out <- cbind(loadings1.out, loadings1)
    # if(bestTune > 1){
    #   coef2.out <- cbind(coef2.out, coef2)
    # }
  }
  r2.out <- c(r2.out, r2)
  bestTune.out <- c(bestTune.out, bestTune[[1]])
  predictor.var.out <- c(predictor.var.out, predictor.var)
  outcome.var.out <- c(outcome.var.out, outcome.var)
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


# Loadings 1st component
loadings1.out <- as.data.frame(loadings1.out)
rownames(loadings1.out) <- rownames(imp.summarised)
names(loadings1.out) <- 1:100
loadings1.means <- rowMeans(loadings1.out)
loadings1.sds <- rowSds(as.matrix(loadings1.out))
loadings1.summarised <- as.data.frame(cbind(loadings1.means, loadings1.sds))
loadings1.summarised$x <- rownames(imp.summarised)


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
# 
# # Variable importance on the projection
# imp.summarised <- imp.summarised %>% 
#   mutate(fill = as.numeric(imp.means < 1))
# 
# ggplot(imp.summarised, aes(x=reorder(x, imp.means, FUN = max), y = imp.means, fill = fill)) + 
#   geom_bar(stat="identity", colour = "black", show.legend = F) +
#   geom_errorbar(aes(ymin=imp.means, ymax=imp.means+imp.sds), width=.2, position=position_dodge(.9)) + 
#   # scale_y_continuous(name = "Variable importance", limits = c(0,25,50,75,100)) +
#   xlab(label = "") + 
#   ylab(label ="Variable Importance in Projection (%)") +
#   coord_flip() + 
#   theme_bw() +
#   theme(text = element_text(size=14),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())
# # ggsave("Figs/Final_Figures_March2021/PLSR/Leaf_variableImportance(n=100).pdf")

# Another plot for Variable importance
imp_plot <- ggplot(data = pivot_longer(as_tibble(t(imp.out)), cols = 1:32), aes(y = value, x = reorder(name, value, FUN = median))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab(label = "") + 
  ylab(label ="Variable Importance in Projection") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))
# ggsave("Figs/Final_Figures_March2021/PLSR/VariableImportance(n=100)_WITH_POSI.pdf")

imp.summarised <- pivot_longer(as_tibble(t(imp.out)), cols = 1:32) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(fill = as.numeric(mean < 1))
  

# PLSR Loadings 1st component
load_plot <- loadings1.summarised %>% 
  # mutate(fill = ifelse(test = x %in% c("N.leaves.shoot", "BG.biomass", "L.lifespan", "L.d13C", "ABG.biomass", "L.lignin", "L.carbon"), 
  #                      yes = 1, no = 2)) %>% 
  left_join(imp.summarised, by = c("x" = "name")) %>% 
  filter(fill == 0) %>% 
  ggplot(aes(x = reorder(x, loadings1.means, FUN = max), y = loadings1.means, fill = fill)) + 
  geom_bar(stat="identity", colour = "black", show.legend = F) +
  geom_errorbar(aes(ymin=loadings1.means - loadings1.sds, ymax = loadings1.means + loadings1.sds), width=.2, position=position_dodge(0.9)) + 
  ylab("PLSR loadings 1st component") +
  xlab("") + 
  # scale_y_continuous(limits = c(-2, 2)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))
# ggsave("Figs/Final_Figures_March2021/PLSR/Leaf_PLSRloadings1stcomp(n=100).pdf")


# Another plot for Loadings importance
# ggplot(data = pivot_longer(as_tibble(t(loadings1.out)), cols = 1:33), aes(y = value, x = reorder(name, value, FUN = median))) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_hline(yintercept = 0, lty = 2) +
#   xlab(label = "") + 
#   ylab(label ="Variable Importance in Projection") +
#   coord_flip() +
#   theme_bw() +
#   theme(text = element_text(size=14),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())
# ggsave("Figs/Final_Figures_March2021/PLSR/VariableImportance(n=100)_WITH_POSI.pdf")

# Figures PLSR scores vs. RESPONSE VARIABLE (similar to Carrascal et al 2009)
# We use just the last result from the for loop as an example
train.data$scoresPlsr1 <- model.initial$finalModel$scores[1:370,1]
plsr_line_plot <- ggplot(data = train.data, mapping = aes(x = scoresPlsr1, y = cstocks)) +
  geom_point(aes(colour = Species)) +
  # scale_colour_simpsons(alpha = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  geom_smooth(method = "lm", se=T, colour = "black") + 
  xlab("PLSR scores 1st component") + 
  ylab(bquote(~C[20]~ 'stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# file_name <- paste("Figs/Final_Figures_March2021/PLSR/Leaf_SpeciesVSplsr1_iOK", i, ".pdf", sep = "")
# ggsave(filename = file_name)

library(cowplot)
panel1 <- plot_grid(imp_plot, 
          load_plot + geom_hline(yintercept = 0),
          nrow = 1, labels = "AUTO", rel_widths = c(1,1), align = "h")

panel2 <- plot_grid(NULL,
          plsr_line_plot + theme(legend.position = "none"),
          NULL,
          ncol = 3, rel_widths = c(1,3,1), labels = c("", "C", ""))
plot_grid(panel1,
          panel2,
          nrow = 2, rel_heights = c(3,2), align = "hv")
# ggsave(filename = "Figs/Final_Figures_March2021/PanelPLSR_sizedNODEPTH_NOPOSI.pdf", width = 185, height = 185, units = "mm")


# Final model summary results
summary(model.initial)
# WITH POSIDONIA OCEANICA
# Data: 	X dimension: 432 32 
# Y dimension: 432 1
# Fit method: simpls
# Number of components considered: 13
# TRAINING: % variance explained
#           1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps  8 comps  9 comps  10 comps  11 comps  12 comps  13 comps
# X           39.14    48.58    55.37    64.46    70.91    76.22    78.61    81.27    83.21     85.84     88.03     89.60     92.17
# .outcome    27.26    44.24    46.73    47.61    48.21    48.47    48.74    48.89    49.01     49.03     49.06     49.08     49.09

# WITHOUT POSIDONIA OCEANICA
# Data: 	X dimension: 370 32 
# Y dimension: 370 1
# Fit method: simpls
# Number of components considered: 15
# TRAINING: % variance explained
#           1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps  8 comps  9 comps  10 comps  11 comps  12 comps  13 comps
# X          36.112    48.46    57.79    64.72    70.22    72.70    76.67    79.85    81.08     84.36     87.58     88.88     
# .outcome    7.858    17.11    22.02    23.83    25.13    26.47    26.99    27.28    27.50     27.57     27.62     27.72


