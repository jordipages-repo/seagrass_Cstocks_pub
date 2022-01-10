# # # # # #
# Analysing the Forqurean's data set provided by Hilary Kennedy
# Jordi F. Pagès
# 15-10-2019
# University of Barcelona
# # # # # # 


# # # # # # # # #
# LIBRARIES ----
# # # # # # # # #

library(tidyverse)
library(tidylog)
library(stringr)
library(forcats)
library(RColorBrewer)
library(metafor)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loading the main data set and checking for errors ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

forqurean <- read_csv(file = "Data_forqurean.csv")
glimpse(forqurean)

source("01_DataImport&Corrections_CarbonReview.R")
rm(list = c("cstocks_tidy", "cstocks_tidy20Stocks"))

# We will only keep cores with more than 3 depths
forqurean_enough <- forqurean %>% 
  left_join(cstocks, key = "CoreID") %>% # This left_join is to add info about the study source from where we took core info
  group_by(CoreID) %>% 
  summarise(n = n(),
            Species = first(Species),
            Source = first(Source)) %>% 
  filter(n>3)

# And join with original data set to recover all columns (because by summarising we lost some columns), and we also join 
# with cstocks to get info about the study source
forqureanOK <- forqurean_enough %>% 
  left_join(forqurean, key = "CoreID") %>% 
  mutate(Carbon_density = 0.01*Organic_carbon*Dry_bulk_density)

 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Checking if all cores have lower %C as depth increases ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# We'll first nest data to then use purrr
by_coreID <- forqureanOK %>% 
  group_by(CoreID) %>% 
  nest()

# Model fitting function (just a linear model) for organic carbon
core_model_organic_carbon <- function(df){
  lm(Organic_carbon ~ Depth_centre_slice, data = df)
}

# Model fitting function (just a linear model) for dry bulk density
core_model_bulkd <- function(df){
  lm(Dry_bulk_density ~ Depth_centre_slice, data = df)
}

# Model fitting function (just a linear model) for carbon density
core_model_carbon_density <- function(df){
  lm(Carbon_density ~ Depth_centre_slice, data = df)
}

# We now apply the function (linear model) to each data frame in each row of the nested data frame
by_coreID <- by_coreID %>% 
  mutate(model_organic_carbon = map(data, core_model_organic_carbon),
         model_bulkd = map(data, core_model_bulkd),
         model_carbon_density = map(data, core_model_carbon_density))

# Now we want the coefficients. We use broom:tidy
by_coreID <- by_coreID %>% 
  mutate(coefs_organic_carbon = map(model_organic_carbon, broom::tidy),
         coefs_bulkd = map(model_bulkd, broom::tidy),
         coefs_carbon_density = map(model_carbon_density, broom::tidy)) %>% 
  unnest(c(coefs_organic_carbon, coefs_bulkd, coefs_carbon_density), .drop = T) %>% 
  rename(term_organic_carbon = term,
         estimate_organic_carbon = estimate,
         std.error_organic_carbon = std.error,
         statistic_organic_carbon = statistic,
         p.value_organic_carbon = p.value,
         term_bulkd = term1,
         estimate_bulkd = estimate1,
         std.error_bulkd = std.error1,
         statistic_bulkd = statistic1,
         p.value_bulkd = p.value1,
         term_carbon_density = term2,
         estimate_carbon_density = estimate2,
         std.error_carbon_density = std.error2,
         statistic_carbon_density = statistic2,
         p.value_carbon_density = p.value2)

# We only want the Depth_centre_slice coeffs
all_coefs <- by_coreID %>% 
  select(-data, -model_organic_carbon, -model_bulkd, -model_carbon_density) %>% 
  filter(term_organic_carbon == "Depth_centre_slice" & term_bulkd == "Depth_centre_slice" & term_carbon_density == "Depth_centre_slice") %>% 
  print(n = Inf)



# # # # # # # # # # # # # # # # # # # # # # # #
# Meta-analysis to check overall           ----    
# significance of slopes of % organic carbon  #
# vs depth and dry bulk density vs depth      #
# # # # # # # # # # # # # # # # # # # # # # # #

# we need the slopes, which can directly be used as effect sizes (Koricheva book)
# and also the variance, so we calculate it from the std.errors
all_coefs <- all_coefs %>% 
  left_join(forqurean_enough, key = "CoreID") %>% 
  mutate(variance_organic_carbon = (std.error_organic_carbon^2)*n,
         variance_bulkd = (std.error_bulkd^2)*n,
         variance_carbon_density = (std.error_carbon_density^2)*n,
         Species = factor(Species),
         Source = ifelse(is.na(Source), "not found", Source))


# # # # # # # # # # # # # # # # # # # 
# ORGANIC CARBON META-ANALYSIS ----
# # # # # # # # # # # # # # # # # # #

# Fixed Effects Meta-analysis for organic carbon
organic_carbon_meta <- rma(yi = estimate_organic_carbon, 
                           vi = variance_organic_carbon, 
                           method = "FE", 
                           data = all_coefs)
summary(organic_carbon_meta)
# Fixed-Effects Model (k = 254)
# 
# logLik    deviance         AIC         BIC        AICc 
# 686.0280    173.6426  -1370.0560  -1366.5186  -1370.0401   
# 
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  0.69
# 
# Test for Heterogeneity:
#   Q(df = 253) = 173.6426, p-val = 1.0000
# 
# Model Results:
#   
# estimate    se     zval    pval    ci.lb   ci.ub 
# -0.0007  0.0005  -1.5919  0.1114  -0.0017  0.0002    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test for heterogeneity validates this fixed effects model, because p-value is high.

# However, we have reasons to think there might be study level heterogeneity. LEt's check it.
# Random Effects Meta-analysis for organic carbon
organic_carbon_meta2 <- rma(yi = estimate_organic_carbon, 
                           vi = variance_organic_carbon, 
                           method = "REML", 
                           test = "knha",
                           data = all_coefs)
summary(organic_carbon_meta2)
# Random-Effects Model (k = 254; tau^2 estimator: REML)
# 
# logLik    deviance         AIC         BIC        AICc 
# 685.5885  -1371.1769  -1367.1769  -1360.1101  -1367.1289   
# 
# tau^2 (estimated amount of total heterogeneity): 0.0000 (SE = 0.0000)
# tau (square root of estimated tau^2 value):      0.0027
# I^2 (total heterogeneity / total variability):   11.44%
# H^2 (total variability / sampling variability):  1.13
# 
# Test for Heterogeneity:
#   Q(df = 253) = 173.6426, p-val = 1.0000
# 
# Model Results:
#   
#   estimate      se     tval    pval    ci.lb    ci.ub 
# -0.0017  0.0005  -3.4322  0.0007  -0.0027  -0.0007  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Adding species as a moderator
organic_carbon_meta3 <- rma(yi = estimate_organic_carbon, 
                            vi = variance_organic_carbon, 
                            mods = ~Species,
                            method = "REML", 
                            test = "knha",
                            data = all_coefs)
summary(organic_carbon_meta3)
# Species is not a significant moderator.

# Adding study as a moderator, to see if we need to include it as a random effect to deal with NON INDEPENDENCE OF STUDIES WITHIN PAPERS
organic_carbon_meta4 <- rma.mv(yi = estimate_organic_carbon, 
                               V = variance_organic_carbon,
                               mods = ~Source,
                               method = "REML", 
                               data = all_coefs)
summary(organic_carbon_meta4)
# Source is significant P = 0.0503

# Let's add it as a random effect.
# Adding a random effect (Source), to account for all those studies coming from the same paper. We first need to take NA's out.
organic_carbon_meta5 <- rma.mv(yi = estimate_organic_carbon,
                        V = variance_organic_carbon,
                        random = ~ 1|Source,
                        method = "REML",
                        data = all_coefs)
summary(organic_carbon_meta5)
# Multivariate Meta-Analysis Model (k = 254; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc 
# 687.7859  -1375.5717  -1371.5717  -1364.5049  -1371.5237   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed  factor 
# sigma^2    0.0000  0.0020     17     no  Source 
# 
# Test for Heterogeneity:
#   Q(df = 253) = 173.6426, p-val = 1.0000
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub 
# -0.0017  0.0009  -1.9844  0.0472  -0.0034  -0.0000  * 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Now let's include moderator = Species
organic_carbon_meta6 <- rma.mv(yi = estimate_organic_carbon,
                               V = variance_organic_carbon,
                               mods = ~Species,
                               random = ~ 1|Source,
                               method = "REML",
                               data = all_coefs)
summary(organic_carbon_meta6)
# Species not significant moderator

# Thus, the best model, is the one with Source as random effects, but without species as moderator.
# That's organic_carbon_meta5

# Funnel plots
funnel(organic_carbon_meta5)
funnel(trimfill(organic_carbon_meta))


# # # # # # # # # # # # # # # # # # 
# BULK DENSITY META-ANALYSIS ----
# # # # # # # # # # # # # # # # # #

# Fixed effects meta-analysis for BULK DENSITY
bulkd_meta <- rma(yi = estimate_bulkd, 
                  vi = variance_bulkd,
                  method = "FE",
                  data = all_coefs)
summary(bulkd_meta)
# Fixed-Effects Model (k = 254)
# 
# logLik    deviance         AIC         BIC        AICc 
# 864.2459    209.1145  -1726.4918  -1722.9545  -1726.4760   
# 
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  0.83
# 
# Test for Heterogeneity:
#   Q(df = 253) = 209.1145, p-val = 0.9797
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub 
# 0.0026  0.0003  7.8872  <.0001  0.0020  0.0033  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test for heterogeneity validates this fixed effects model, because p-value is high.
# However, we have reasons to think there might be study level heterogeneity. Let's check it.
# Random Effects Meta-analysis for organic carbon
bulkd_meta2 <- rma(yi = estimate_bulkd, 
                   vi = variance_bulkd, 
                   method = "REML", 
                   test = "knha",
                   data = all_coefs)
summary(bulkd_meta2)
# Random-Effects Model (k = 254; tau^2 estimator: REML)
# 
# logLik    deviance         AIC         BIC        AICc 
# 870.7110  -1741.4221  -1737.4221  -1730.3553  -1737.3741   
# 
# tau^2 (estimated amount of total heterogeneity): 0.0000 (SE = 0.0000)
# tau (square root of estimated tau^2 value):      0.0030
# I^2 (total heterogeneity / total variability):   24.09%
# H^2 (total variability / sampling variability):  1.32
# 
# Test for Heterogeneity:
#   Q(df = 253) = 209.1145, p-val = 0.9797
# 
# Model Results:
#   
#   estimate      se    tval    pval   ci.lb   ci.ub 
# 0.0033  0.0003  9.7829  <.0001  0.0027  0.0040  *** 

# There is some heterogeneity, 24%

# Adding species as a moderator
bulkd_meta3 <- rma(yi = estimate_bulkd,
                   vi = variance_bulkd, 
                   mods = ~Species,
                   method = "REML", 
                   test = "knha",
                   data = all_coefs)
summary(bulkd_meta3)
# Species appears to be significant here.

# Adding study as a moderator, to see if we need to include it as a random effect to deal with NON INDEPENDENCE OF STUDIES WITHIN PAPERS
bulkd_meta4 <- rma.mv(yi = estimate_bulkd, 
                      V = variance_bulkd,
                      mods = ~Source,
                      method = "REML", 
                      data = all_coefs)
summary(bulkd_meta4)
# Source is significant P = 0.0001

# Let's add it as a random effect.
# Adding a random effect (Source), to account for all those studies coming from the same paper. We first need to take NA's out.
bulkd_meta5 <- rma.mv(yi = estimate_bulkd,
                      V = variance_bulkd,
                      random = ~ 1|Source,
                      method = "REML",
                      data = all_coefs)
summary(bulkd_meta5)
# Multivariate Meta-Analysis Model (k = 254; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc 
# 866.0125  -1732.0250  -1728.0250  -1720.9582  -1727.9770   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed  factor 
# sigma^2    0.0000  0.0026     17     no  Source 
# 
# Test for Heterogeneity:
#   Q(df = 253) = 209.1145, p-val = 0.9797
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub 
# 0.0037  0.0008  4.4554  <.0001  0.0021  0.0053  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Now let's include moderator = Species
bulkd_meta6 <- rma.mv(yi = estimate_bulkd,
                      V = variance_bulkd,
                      mods = ~Species,
                      random = ~ 1|Source,
                      method = "REML",
                      data = all_coefs)
summary(bulkd_meta6)
# Species has to be included.

# bulkd_meta6 is the best-selected model.

# Funnel plots
funnel(bulkd_meta5)
funnel(trimfill(bulkd_meta))


# # # # # # # # # # # # # # # # # # 
# CARBON DENSITY META-ANALYSIS ----
# # # # # # # # # # # # # # # # # #

# Fixed effects meta-analysis for CARBON DENSITY
carbon_density_meta <- rma(yi = estimate_carbon_density, 
                           vi = variance_carbon_density, 
                           method = "FE",
                           data = all_coefs)
summary(carbon_density_meta)
# Fixed-Effects Model (k = 254)
# 
# logLik    deviance         AIC         BIC        AICc 
# 1868.2117    131.5190  -3734.4234  -3730.8860  -3734.4075   
# 
# I^2 (total heterogeneity / total variability):   0.00%
# H^2 (total variability / sampling variability):  0.52
# 
# Test for Heterogeneity:
#   Q(df = 253) = 131.5190, p-val = 1.0000
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub 
# 0.0000  0.0000  0.5161  0.6058  -0.0000  0.0000      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Test for heterogeneity validates this fixed effects model, because p-value is high.

# However, we have reasons to think there might be study level heterogeneity. LEt's check it.
# Random Effects Meta-analysis for organic carbon
carbon_density_meta2 <- rma(yi = estimate_carbon_density, 
                            vi = variance_carbon_density, 
                            method = "REML", 
                            test = "knha",
                            data = all_coefs)
summary(carbon_density_meta2)
# Random-Effects Model (k = 254; tau^2 estimator: REML)
# 
# logLik    deviance         AIC         BIC        AICc 
# 1860.1405  -3720.2811  -3716.2811  -3709.2143  -3716.2331   
# 
# tau^2 (estimated amount of total heterogeneity): 0.0000 (SE = 0.0000)
# tau (square root of estimated tau^2 value):      0.0000
# I^2 (total heterogeneity / total variability):   5.53%
# H^2 (total variability / sampling variability):  1.06
# 
# Test for Heterogeneity:
#   Q(df = 253) = 131.5190, p-val = 1.0000
# 
# Model Results:
#   
#   estimate      se     tval    pval    ci.lb   ci.ub 
# -0.0000  0.0000  -0.3498  0.7267  -0.0000  0.0000    
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Adding species as a moderator
carbon_density_meta3 <- rma(yi = estimate_carbon_density, 
                            vi = variance_carbon_density, 
                            mods = ~Species,
                            method = "REML", 
                            test = "knha",
                            data = all_coefs)
summary(carbon_density_meta3)
# Species is a significant moderator.

# Adding study as a moderator, to see if we need to include it as a random effect to deal with NON INDEPENDENCE OF STUDIES WITHIN PAPERS
carbon_density_meta4 <- rma.mv(yi = estimate_carbon_density, 
                               V = variance_carbon_density,
                               mods = ~Source,
                               method = "REML", 
                               data = all_coefs)
summary(carbon_density_meta4)
# Source is significant P = 0.0274

# Let's add it as a random effect.
# Adding a random effect (Source), to account for all those studies coming from the same paper. We first need to take NA's out.
carbon_density_meta5 <- rma.mv(yi = estimate_carbon_density,
                               V = variance_carbon_density,
                               random = ~ 1|Source,
                               method = "REML",
                               data = all_coefs)
summary(carbon_density_meta5)
# Multivariate Meta-Analysis Model (k = 254; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc 
# 1861.6324  -3723.2649  -3719.2649  -3712.1981  -3719.2169   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed  factor 
# sigma^2    0.0000  0.0000     17     no  Source 
# 
# Test for Heterogeneity:
#   Q(df = 253) = 131.5190, p-val = 1.0000
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub 
# -0.0000  0.0000  -0.3695  0.7117  -0.0000  0.0000    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Now let's include moderator = Species
carbon_density_meta6 <- rma.mv(yi = estimate_carbon_density,
                               V = variance_carbon_density,
                               mods = ~Species,
                               random = ~ 1|Source,
                               method = "REML",
                               data = all_coefs)
summary(carbon_density_meta6)
# Species not significant moderator

# Thus, the best model, is the one with Source as random effects, but without species as moderator.
# That's carbon_density_meta5

# Funnel plots
funnel(carbon_density_meta5)
funnel(trimfill(carbon_density_meta))




# Data to be plotted
data_organic_carbon <- tibble("Organic carbon", 
                              organic_carbon_meta5$beta, 
                              organic_carbon_meta5$ci.lb, 
                              organic_carbon_meta5$ci.ub, 
                              "*")
names(data_organic_carbon) <- c("variable", 
                                "estimate", 
                                "lower", 
                                "upper", 
                                "label")
data_bulk <- tibble("Dry bulk density", 
                    bulkd_meta5$beta,
                    bulkd_meta5$ci.lb,
                    bulkd_meta5$ci.ub,
                    "***")
names(data_bulk) <- c("variable", 
                      "estimate",
                      "lower", 
                      "upper", 
                      "label")
data_carbon_density <- tibble("Carbon density",
                              carbon_density_meta5$beta,
                              carbon_density_meta5$ci.lb, 
                              carbon_density_meta5$ci.ub,
                              " ")
names(data_carbon_density) <- c("variable", 
                                "estimate", 
                                "lower", 
                                "upper", 
                                "label")
data_meta <- rbind(data_organic_carbon, data_bulk, data_carbon_density)

# Nice forest plot for both carbon and bulk density
ggplot(data_meta, aes(x = variable, y = estimate, ymax = upper, ymin = lower)) +
  geom_point(aes(size = 1.2)) +
  geom_errorbar(aes(width = 0.1)) +
  geom_text(aes(label = label), size = 8, position = position_nudge(x = 0.1)) +
  coord_flip() +
  geom_hline(aes(yintercept = 0), lty = 2, size = 1) + 
  xlab("") +
  ylab("Change along core depth") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/Final_Figures_March2021/Meta-analysis.pdf")





