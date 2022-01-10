# # # # # #
# Analysing the data set provided by Hilary Kennedy
# Jordi F. Pagès
# 30-07-2019
# CEAB
# # # # # # 

library(nlme)
library(car)
library(modelr)
library(multcomp)
# library(multcompView)
# library(emmeans)

# # # # # # # # # # # # # # # # # # # # # 
# Loading the data clean and corrected  ----
# # # # # # # # # # # # # # # # # # # # # 

source("1_DataImport&Corrections_CarbonReview.R")
source("mcheck_function.R")

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
  # filter(Species != "Posidonia oceanica") %>% # To check the effect of Species without posidonia oceanica
  # select(Latitude, Species, Source, depth, cstocks, CoreID_unique) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  
# filter(Species != "Posidonia oceanica" | cstocks != 1.9)

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Modelling SPECIES, DEPTH, FIN_TYP_names (geomorphology) and LAT_ZONE (affinity) - OPTION 1. ALL CORES (the GOOD one) ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# USING A BETTER APPROACH AT THE END OF THIS SCRIPT! (MODELLING HETEROGENEITY AT THE FIN_TYP_names level too)

# Removing species with less than 10 samples
SpeciesMoreThan10Samples <- cstocks_data %>% 
  group_by(Species) %>% 
  summarise(n = n()) %>% 
  filter(n>=10) %>% 
  select(Species)
SpeciesMoreThan10Samples <- as.data.frame(SpeciesMoreThan10Samples)
SpeciesMoreThan10Samples <- SpeciesMoreThan10Samples$Species

cstocks_data <- cstocks_data %>% 
  filter(Species %in% SpeciesMoreThan10Samples)

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$Species <- droplevels(cstocks_data$Species) 


# Full linear model
m0 <- lm(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_data)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We do need weights!

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
anova(m0.w.lme, m0.gls.w) # Source is needed

m0.w.lme1 <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 random = ~1|CoreID_unique,
                 weights = varIdent(form = ~1|Species), 
                 control = lmeControl(maxIter = 100, msMaxIter = 100),
                 data = cstocks_data)
anova(m0.w.lme1, m0.gls.w) # CoreID_unique is needed
mcheck(m0.w.lme1)
# We know that not including this is a bit like committing pseudo-replication, but this variable messes up model validation
# therefore, we are not going to include it. Instead, we will run similar models with subsets of this data set to 
# ensure that pseudo-replication is not an issue.

# m0.w.lme2 <- lme(log(cstocks) ~ Species + depth, 
#                 random = ~1|Source/CoreID_unique,
#                 weights = varIdent(form = ~1|Species), 
#                 control = lmeControl(maxIter = 100, msMaxIter = 100),
#                 data = cstocks_data)
# anova(m0.w.lme2, m0.gls.w) # Both random effects are needed.


# # 2 random effects... but how should we specify them exactly?
# # Let's start without weights
# m0.w.lmer1 <- lme4::lmer(log(cstocks) ~ Species + depth + (1|Source),
#                         data = cstocks_data)
# m0.w.lmer2 <- lme4::lmer(log(cstocks) ~ Species + depth + (1|Source) + (1|CoreID_unique),
#                          data = cstocks_data)
# m0.w.lmer3 <- lme4::lmer(log(cstocks) ~ Species + depth + (1|Source/CoreID_unique),
#                          data = cstocks_data)
# anova(m0.w.lmer1, m0.w.lmer2) # 2random effects better than 1
# anova(m0.w.lmer1, m0.w.lmer3) # 2 random effects better than 1
# anova(m0.w.lmer2, m0.w.lmer3) # No significant difference between crossed random effects or nested random effects. 
# 
# So the nlme version with nested random effects is the best.


# As said, although including coreID makes sense, it messes up model validation, and 
# thus we are not including it.

# Model selection 
m1.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -Species)
anova(m1.w.lme, m2.w.lme) # Species CANNOT be dropped.

m4.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m4.w.lme) # depth can be dropped... 

m5.w.lme <- update(m4.w.lme, ~. -FIN_TYP_names)
anova(m4.w.lme, m5.w.lme) # FIN_TYP_names can be dropped... 

m6.w.lme <- update(m5.w.lme, ~. -Lat_Zone)
anova(m5.w.lme, m6.w.lme) # Lat_Zone can be dropped (but p-value = 0.0647)... 

mfinal <- lme(log(cstocks) ~ Species, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species),
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
Anova(mfinal)
summary(mfinal)
mcheck(mfinal) # Residuals in qqplot are good! (if we do not include coreID)


# For another analysis with latitude later
cstocks_data_clean_Species <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes")

# Checking for extreme values. deleting residuals bigger 1.5
cstocks_data_clean <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred") %>% 
  filter(abs(logSpeciesRes)<1.5) # 7 extreme values deleted (1%)

cstocks_data_clean$Species <- droplevels(cstocks_data_clean$Species) 

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_data_clean)
m0.gls.w <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data_clean)
anova(m0.gls, m0.gls.w) # We do need weights!

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_clean)
anova(m0.gls.w, m0.w.lme) # The random effect Source is needed

# Model selection 
m1.w.lme <- update(m0.w.lme, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -Species)
anova(m1.w.lme, m2.w.lme) # Species CANNOT be dropped.

m4.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m4.w.lme) # depth CAN be dropped.

m5.w.lme <- update(m4.w.lme, ~. -FIN_TYP_names)
anova(m4.w.lme, m5.w.lme) # FIN_TYP_names CAN be dropped.

m6.w.lme <- update(m5.w.lme, ~. -Lat_Zone)
anova(m5.w.lme, m6.w.lme) # Lat_Zone CAN be dropped.


mfinal <- lme(log(cstocks) ~ Species, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data_clean)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
#           Chisq Df Pr(>Chisq)    
# Species 218.94  16  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(mfinal)

cstocks_data_clean <- cstocks_data_clean %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred")

# Model validation
mcheck(mfinal) #  OK FIT

# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(mfinal))

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(cstocks_data_clean$Species))

# Multiple comparisons with library(multcompar)
multcompar <- glht(mfinal, linfct=mcp(Species="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar))
# Significant multiple comparisons
summary_multcompar %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf) %>% 
  write_csv(file = "Figs/table_multcompar_species.csv")
# A tibble: 19 x 7
# term    contrast                                          null.value estimate std.error statistic adj.p.value
# 1 Species Posidonia oceanica - Amphibolis antarctica             0    1.69      0.316      5.36 0.0000220   
# 2 Species Posidonia oceanica - Cymodocea rotundata               0    1.46      0.296      4.92 0.0000603   
# 3 Species Halophila stipulacea - Cymodocea serrulata             0   -0.420     0.120     -3.48 0.0311      
# 4 Species Posidonia oceanica - Cymodocea serrulata               0    1.23      0.265      4.65 0.000307    
# 5 Species Halophila stipulacea - Enhalus acoroides               0   -0.417     0.106     -3.94 0.00616     
# 6 Species Posidonia oceanica - Enhalus acoroides                 0    1.24      0.264      4.68 0.000290    
# 7 Species Posidonia oceanica - Halodule uninervis                0    1.46      0.271      5.37 0.00000782  
# 8 Species Posidonia australis - Halophila ovalis                 0    0.706     0.179      3.95 0.00573     
# 9 Species Posidonia oceanica - Halophila ovalis                  0    1.89      0.313      6.04 0.0000000436
# 10 Species Posidonia oceanica - Halophila stipulacea              0    1.65      0.272      6.07 0.0000000305
# 11 Species Thalassia hemprichii - Halophila stipulacea            0    0.404     0.101      3.99 0.00520     
# 12 Species Posidonia oceanica - Posidonia australis               0    1.18      0.279      4.25 0.00199     
# 13 Species Posidonia sinuosa - Posidonia oceanica                 0   -1.44      0.292     -4.92 0.0000783   
# 14 Species Thalassia hemprichii - Posidonia oceanica              0   -1.25      0.263     -4.75 0.000117    
# 15 Species Thalassia testudinum - Posidonia oceanica              0   -1.01      0.115     -8.80 0           
# 16 Species Thalassodendron ciliatum - Posidonia oceanica          0   -1.78      0.306     -5.81 0.000000410 
# 17 Species Zostera marina - Posidonia oceanica                    0   -1.20      0.249     -4.83 0.000175    
# 18 Species Zostera noltii - Posidonia oceanica                    0   -1.49      0.126    -11.8  0           
# 19 Species Zostera noltii - Thalassia testudinum                  0   -0.478     0.130     -3.67 0.0171 

# I have checked it, and these multiple comparisons do not change if we include CoreID_unique as a random factor, or not.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THERE IS A STRONG EFFECT OF SPECIES IDENTITY ON CARBON STOCKS.#
# MAINLY BECAUSE OF POSIDONIA, BUT NOT ONLY BECAUSE OF IT.      #
# TO CHECK THE EFFECT OF SPECIES WITHOUT POSIDONIA OCEANICA,    #
# UNCOMMENT FILTER AT BEGINNING OF MODELLING.                   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# So, results:
#     Random effect source is needed (Effect of CoreID_unique is not very important.)
#     Depth does not influence cstocks
#     Species identity INFLUENCES cstocks
#     Geomorphology does NOT influence cstocks
#     Affinity (temperate, polar or tropical) is on the limit of significance... 
# 



# # # # # # # # # # # # # # # # # # # # #
# Modelling SPECIES EFFECT, depths together OPTION 2. ONLY CORES WITH BOTH DEPTHS ----
# # # # # # # # # # # # # # # # # # # # #

# We now filter those cores that have information for both 0-20cm and 20-50cm depths.
# These cores are the ones to which really make sense to apply the random effect coreID. So let's see.
coreID_both_depths <- cstocks_data %>% 
  group_by(CoreID_unique) %>% 
  summarise(n = n()) %>% 
  filter(n>1) %>% 
  select(CoreID_unique)
coreID_both_depths <- as.array(coreID_both_depths$CoreID_unique)
cstocks_data_both_depths <- cstocks_data %>% 
  filter(CoreID_unique %in% coreID_both_depths)

# Full linear model
m0 <- lm(log(cstocks) ~ Species+depth, data = cstocks_data_both_depths)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Species + depth, data = cstocks_data_both_depths)
m0.gls.w <- gls(log(cstocks) ~ Species + depth, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data_both_depths)
anova(m0.gls, m0.gls.w) # We do need weights!


# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_both_depths)
anova(m0.w.lme, m0.gls.w) # Source is needed

m0.w.lme1 <- lme(log(cstocks) ~ Species + depth, 
                 random = ~1|CoreID_unique,
                 weights = varIdent(form = ~1|Species), 
                 control = lmeControl(maxIter = 100, msMaxIter = 100),
                 data = cstocks_data_both_depths)
anova(m0.w.lme1, m0.gls.w) # CoreID_unique is needed
# Again, CoreID messes up validation plots. Therefore, we're not going to use this random effect.
# In this case, it doesn't make that much of sense, because at least in this sense the data set is balanced: 
# we have as much 0-20 as 20-50 cores.

# m0.w.lme2 <- lme(log(cstocks) ~ Species + depth, 
#                  random = ~1|Source/CoreID_unique,
#                  weights = varIdent(form = ~1|Species), 
#                  control = lmeControl(maxIter = 100, msMaxIter = 100),
#                  data = cstocks_data_both_depths)
# anova(m0.w.lme2, m0.gls.w) # Both random effects are needed.

# Model selection 
m1.w.lme <- lme(log(cstocks) ~ Species + depth, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_both_depths, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -Species)
anova(m1.w.lme, m2.w.lme) # Species CANNOT be dropped.

m4.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m4.w.lme) # depth can be dropped... 

mfinal <- lme(log(cstocks) ~ Species, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species),
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data_both_depths)
car::Anova(mfinal) # for some reason it doesn't work with this data set...
anova(mfinal)
summary(mfinal)
mcheck(mfinal) # Residuals in qqplot are good! (if we do not include coreID)

# So, results:
#     Random effect source is needed
#     Depth does not influence cstocks
#     Species identity INFLUENCES cstocks

#####

# # # # # # # # # # # # # # # # # # # # #
# Modelling SPECIES EFFECT, depths together OPTION 3. ONLY 1 SAMPLE PER COREID ----
# # # # # # # # # # # # # # # # # # # # #

# We now take only 1 random sample for those coreIDs that have 2 depths.
coreID_both_depths <- cstocks_data %>% 
  group_by(CoreID_unique) %>% 
  summarise(n = n()) %>% 
  filter(n>1) %>% 
  select(CoreID_unique)
coreID_both_depths <- as.array(coreID_both_depths$CoreID_unique)

# If we take a look at the cores that are NOT in the data set coreID_both_depths...:
cstocks_data_subsample1 <- cstocks_data %>% 
  filter(!CoreID_unique %in% coreID_both_depths)
# We see that all of them are from 0-20 cm depths. 
# Therefore, now, we only want the 20-50cm data for the cores %in% coreID_both_depths
cstocks_data_subsample2 <- cstocks_data %>% 
  filter(CoreID_unique %in% coreID_both_depths) %>% 
  filter(depth == "20-50 cm")

# Now we want both subsamples in a single data set
cstocks_subsampleOK <- bind_rows(cstocks_data_subsample1, cstocks_data_subsample2)

# Just to check that we only have one observation per coreID
cstocks_subsampleOK %>% 
  group_by(CoreID_unique) %>% 
  summarise(n = n()) %>% 
  filter(n>1)
# Good! No cores left with more than 1 sample.

# Modelling per se
# Full linear model
m0 <- lm(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_subsampleOK)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_subsampleOK)
m0.gls.w <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_subsampleOK)
anova(m0.gls, m0.gls.w) # We do need weights!


# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_subsampleOK)
anova(m0.w.lme, m0.gls.w) # Source is needed

# Model selection 
m1.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_subsampleOK, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -Species)
anova(m1.w.lme, m2.w.lme) # Species CANNOT be dropped.

m4.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m4.w.lme) # depth can be dropped... 

m5.w.lme <- update(m4.w.lme, ~. -FIN_TYP_names)
anova(m4.w.lme, m5.w.lme) # FIN_TYP_names can be dropped... 

m6.w.lme <- update(m5.w.lme, ~. -Lat_Zone)
anova(m5.w.lme, m6.w.lme) # Lat_Zone can be dropped... 

mfinal <- lme(log(cstocks) ~ Species, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species),
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_subsampleOK)
Anova(mfinal)
mcheck(mfinal) # Good!
mcheck2(mfinal) # Residuals in qqplot are good! (if we do not include coreID)


# Checking for extreme values. deleting residuals bigger 1.5
cstocks_subsampleOK_clean <- cstocks_subsampleOK %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred") %>% 
  filter(abs(logSpeciesRes)<1.5) # 5 extreme values deleted (1%)

cstocks_subsampleOK_clean$Species <- droplevels(cstocks_subsampleOK_clean$Species) 

mfinal <- lme(log(cstocks) ~ Species, 
             random = ~1|Source,
             weights = varIdent(form = ~1|Species),
             control = lmeControl(maxIter = 100, msMaxIter = 100),
             data = cstocks_subsampleOK_clean)
Anova(mfinal)
mcheck(mfinal) # Good!
mcheck2(mfinal) # Residuals in qqplot are good! (if we do not include coreID)

# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(mfinal))

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(cstocks_subsampleOK_clean$Species))

# Multiple comparisons with library(multcompar)
multcompar <- glht(mfinal, linfct=mcp(Species="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar, level = 0.06))
# Significant multiple comparisons
summary_multcompar %>%
  filter(adj.p.value<0.051) %>%
  print(n = Inf)
# A tibble: 14 x 7
# term    contrast                                      null.value estimate std.error statistic   adj.p.value
#  1 Species Posidonia oceanica - Amphibolis antarctica             0     1.57     0.349      4.51 0.000772     
#  2 Species Posidonia oceanica - Cymodocea rotundata               0     1.31     0.338      3.89 0.00743      
#  3 Species Posidonia oceanica - Cymodocea serrulata               0     1.02     0.297      3.42 0.0398       
#  4 Species Posidonia oceanica - Enhalus acoroides                 0     1.03     0.300      3.45 0.0368       
#  5 Species Posidonia oceanica - Halodule uninervis                0     1.25     0.303      4.13 0.00342      
#  6 Species Posidonia oceanica - Halophila ovalis                  0     1.57     0.364      4.33 0.00155      
#  7 Species Posidonia oceanica - Halophila stipulacea              0     1.45     0.306      4.72 0.000175     
#  8 Species Posidonia oceanica - Posidonia australis               0     1.03     0.307      3.35 0.0509       
#  9 Species Posidonia sinuosa - Posidonia oceanica                 0    -1.29     0.329     -3.93 0.00639      
# 10 Species Thalassia hemprichii - Posidonia oceanica              0    -1.07     0.296     -3.62 0.0207       
# 11 Species Thalassia testudinum - Posidonia oceanica              0    -1.08     0.169     -6.41 0.00000000877
# 12 Species Thalassodendron ciliatum - Posidonia oceanica          0    -1.38     0.351     -3.93 0.00662      
# 13 Species Zostera marina - Posidonia oceanica                    0    -1.14     0.285     -3.99 0.00538      
# 14 Species Zostera noltii - Posidonia oceanica                    0    -1.48     0.146    -10.1  0            


# So, results:
#     Random effect source is needed
#     Depth does not influence cstocks
#     Species identity INFLUENCES cstocks


####
# ALL MODELS DISCARD DEPTH AS SIGNIFICANT PREDICTOR OF CSTOCKS, AND 
# POINT TO SPECIES IDENTITY AS A VERY SIGNIFICANT PREDICTOR.
####

#####



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# BARPLOTS - C stocks by SPECIES by depth with sample size coloured ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# AFER MANY TRIALS THE FINAL MODEL TO PRESENT IN THE PAPER IS 
# mfinal from Modelling [...] - OPTION 1. ALL CORES (the GOOD one)

# We want to get the mean. 
means <- cstocks_data_clean %>% 
  filter(!is.na(cstocks)) %>% 
  # filter(Meadow_type == "monospecific") %>% 
  summarise(n = n(),
            mean = mean(cstocks),
            median = median(cstocks),
            max = max(cstocks))

cstocks_to_order <- cstocks_data_clean %>% 
  # filter(Meadow_type == "monospecific") %>%
  filter(!is.na(cstocks)) %>%
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev(),
         Pred = exp(logSpeciesPred)) %>% 
  group_by(Species) %>% 
  summarise(n = n(),
            mean_cstocks = mean(cstocks),
            error_cstocks = std.error(cstocks),
            mean_pred = mean(Pred),
            error_pred = std.error(Pred)) %>% 
  ungroup() %>% 
  arrange(mean_cstocks) %>% 
  mutate(order = row_number(),
         nscale = ifelse(n<=10, "n≤10", "n>10"))

cstocks_to_order <- left_join(cstocks_to_order, letters_multcompar, by = 'Species')

ggplot(cstocks_to_order, aes(x = order, y = mean_cstocks)) +
  geom_bar(aes(x = order, y = mean_cstocks, fill = nscale), stat = "identity") +
  scale_fill_manual(name = "Sample size", 
                    labels = c("n > 10", bquote("n" <= 10)),
                    values = c("#00BFC4", "#F8766D")) +
  geom_errorbar(aes(ymin = mean_cstocks-error_cstocks, ymax = mean_cstocks+error_cstocks), width = 0.3) +
  geom_hline(data = means, aes(yintercept = mean), lty = 2, colour = "darkgrey") +
  geom_text(aes(y = ifelse(n>1, mean_cstocks+error_cstocks + 0.3, mean_cstocks + 0.3), label = str_c("(", n, ")")), size = 3) +
  geom_text(aes(y = ifelse(n>1, mean_cstocks+error_cstocks + 1, mean_cstocks + 1), label = letters)) +
  # geom_point(aes(y = mean_pred)) +
  # geom_errorbar(aes(ymin = mean_pred - error_pred, ymax = mean_pred + error_pred), width = 0, lwd = 3, colour = "blue") +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  # facet_wrap(~depth, scales = "free") + 
  coord_flip(ylim = c(0,5)) +
  scale_x_continuous(breaks = cstocks_to_order$order,
                     labels = cstocks_to_order$Species,
                     expand = c(0.008,0.008)) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.justification = c(1,0),
        legend.position = c(0.97, 0.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14)) +
  theme(legend.position = "none")
# ggsave(filename = "Figs/Cstocks_by_Species_BARPLOT_TUKEY_ALLCORES_OPTION1.pdf")




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# BOXPLOTS - C stocks by SPECIES by depth with sample size coloured ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# To count the sample size for each species and depth
listSpeciesN <- cstocks_data_clean %>%
  filter(!is.na(cstocks)) %>% 
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev(),
         Pred = exp(logSpeciesPred)) %>% 
  group_by(Species) %>% 
  summarise(n = n(),
            median_stocks = median(cstocks),
            max_stocks = max(cstocks),
            mean_pred = mean(Pred)) %>% 
  ungroup() %>% 
  arrange(median_stocks) %>% 
  mutate(order = row_number(),
         nscale = ifelse(n<=10, "n≤10", "n>10")) %>% 
  left_join(letters_multcompar, by = 'Species')

# We join the above data set with the data set 
cstocksN <- cstocks_data_clean %>% 
  filter(!is.na(cstocks)) %>% 
  left_join(listSpeciesN, by = c("Species"))

ggplot(cstocksN, aes(x = as.factor(order), y = cstocks)) +
  # geom_boxplot(aes(fill = nscale, colour = nscale), fatten = 1) +
  geom_violin(aes(fill = nscale, colour = nscale)) +
  scale_fill_manual(name = "Sample size", 
                    labels = c("n > 10", bquote("n" <= 10)),
                    values = c("#00BFC4", "#F8766D")) +
  scale_colour_manual(name = "Sample size", 
                      labels = c("n > 10", bquote("n" <= 10)),
                      values = c("#00BFC4", "#F8766D")) +
  geom_point(aes(y = median_stocks), shape = 15, size = 1.5) +
  geom_hline(data = means, aes(yintercept = median), lty = 2, colour = "darkgrey") +
  geom_text(data = listSpeciesN, aes(x = as.factor(order), y = max_stocks + 1, label = str_c("(", n, ")")), size = 3) +
  geom_text(data = listSpeciesN, aes(y = ifelse(n>1, max_stocks + 2, max_stocks + 2), label = letters)) +
  # geom_point(aes(y = mean_pred)) +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip(ylim = c(0, 15)) +
  scale_x_discrete(breaks = listSpeciesN$order,
                   labels = listSpeciesN$Species,
                   expand = c(0.025,0.025)) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.justification = c(1,0),
        legend.position = c(0.97, 0.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14)) +
  theme(legend.position = "none")
# ggsave(filename = "Figs/Cstocks_by_Species_VIOLINPLOT_TUKEY_ALLCORES_OPTION1.pdf")



# # # # # # # # # # # # # # # # # # # # #
# Modelling DEPTH EFFECT, species together ----
# # # # # # # # # # # # # # # # # # # # #

# cstocks_data <- cstocks_tidy %>%
#   filter(Meadow_type == "monospecific") %>%
#   select(Species, depth, Source, cstocks) %>%
#   filter_all(all_vars(!is.na(.)))

# Full linear model
m0 <- lm(log(cstocks) ~ depth, data = cstocks_data)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ depth, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ depth,
                weights = varIdent(form = ~1|Species),
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We do need weights!

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ depth,
                random = ~1|Source,
                weights = varIdent(form = ~1|Species),
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data)
anova(m0.gls.w, m0.w.lme) # The random effect Source is needed

mfinal <- lme(log(cstocks) ~ depth,
              random = ~1|Source,
              weights = varIdent(form = ~1|Species),
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
Anova(mfinal)
summary(mfinal)

# Model validation
mcheck(mfinal) # there are some very high residuals.
# min(resid(mfinal))

# Checking for extreme values. deleting residuals bigger 1.5
cstocks_data_clean <- cstocks_data %>%
  add_residuals(model = mfinal, var = "logDepthRes") %>%
  add_predictions(model = mfinal, var = "logDepthPred") %>%
  filter(abs(logDepthRes)<1.5) # 7 extreme values deleted (1%)


cstocks_data_clean$Species <- droplevels(cstocks_data_clean$Species)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ depth, data = cstocks_data_clean)
m0.gls.w <- gls(log(cstocks) ~ depth,
                weights = varIdent(form = ~1|Species),
                data = cstocks_data_clean)
anova(m0.gls, m0.gls.w) # We do need weights!

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ depth,
                random = ~1|Source,
                weights = varIdent(form = ~1|Species),
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_clean)
anova(m0.gls.w, m0.w.lme) # The random effect Source is needed

mfinal <- lme(log(cstocks) ~ depth,
              random = ~1|Source,
              weights = varIdent(form = ~1|Species),
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data_clean)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
#
# Response: log(cstocks)
#        Chisq Df Pr(>Chisq)
# depth 0.7989  1     0.3714
summary(mfinal)

# Model validation
mcheck(mfinal) #  OK FIT

cstocks_data_clean <- cstocks_data_clean %>%
  add_residuals(model = mfinal, var = "logDepthRes") %>%
  add_predictions(model = mfinal, var = "logDepthPred")


# Multiple comparisons with library(multcompar)
multcompar <- glht(mfinal, linfct=mcp(depth="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar))
# Significant multiple comparisons
summary_multcompar %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DEPTH DOES NOT INFLUENCE CARBON STOCKS. NO DOUBT ABOUT IT.  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# BARPLOTS - C stocks by DEPTH with sample size coloured  ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cstocks_data_clean %>% 
#   mutate(Pred = exp(logDepthPred)) %>% 
#   group_by(depth) %>% 
#   summarise(n = n(),
#             mean_cstocks = mean(cstocks),
#             error_cstocks = std.error(cstocks),
#             mean_pred = mean(Pred),
#             error_pred = std.error(Pred)) %>% 
#   left_join(letters_multcompar, by = "depth") %>% 
#   ggplot(aes(x = rev(depth), y = mean_cstocks, fill = depth)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("#189ad3", "#005073")) +
#   scale_x_discrete(labels = rev(c("0-20 cm", "20-50 cm"))) +
#   geom_errorbar(aes(ymin = mean_cstocks-error_cstocks, ymax = mean_cstocks+error_cstocks), width = 0.3) +
#   geom_point(aes(y = mean_pred)) +
#   # geom_errorbar(aes(ymin = mean_pred-error_pred, ymax = mean_pred+error_pred), width = 0.3) +
#   geom_text(aes(y = mean_cstocks+error_cstocks + 0.5, label = str_c("(", n, ")")), size = 3) +
#   geom_text(aes(y = mean_cstocks+error_cstocks + 1, label = letters)) +
#   xlab("") +
#   ylab(bquote('Carbon density (Mg C' ~ha^-1~cm^-1* ')')) +
#   coord_flip() +
#   theme_bw() +
#   theme(legend.position = 'none',
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 14))
# # ggsave(filename = "Figs/Cstocks_by_depth_barplot_Tukey_fitted.pdf")  
# 
# 
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # BOXPLOTS - C stocks by DEPTH with sample size coloured  ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# To count the sample size for each species and depth
listSpeciesN <- cstocks_data_clean %>%
  filter(!is.na(cstocks)) %>%
  mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev(),
         Pred = exp(logDepthPred)) %>%
  group_by(depth) %>%
  summarise(n = n(),
            median_stocks = median(cstocks),
            max_stocks = max(cstocks),
            mean_pred = mean(Pred)) %>%
  left_join(letters_multcompar, by = 'depth')

# We join the above data set with the data set
cstocksN <- cstocks_data_clean %>%
  filter(!is.na(cstocks)) %>%
  left_join(listSpeciesN, by = c("depth"))

# We want to get the mean.
means <- cstocks_data_clean %>%
  filter(!is.na(cstocks)) %>%
  # filter(Meadow_type == "monospecific") %>%
  summarise(n = n(),
            mean = mean(cstocks),
            median = median(cstocks),
            max = max(cstocks))

cstocksN %>%
  mutate(depth = factor(depth,
                        levels = c("20-50 cm", "0-20 cm"))) %>%
  ggplot(aes(x = depth, y = cstocks)) +
  geom_boxplot(aes(fill = depth, colour = depth)) +
  scale_fill_manual(values = rev(c("#189ad3", "#005073"))) +
  scale_colour_manual(values = rev(c("#189ad3", "#005073"))) +
  geom_hline(data = means, aes(yintercept = median), lty = 2, colour = "darkgrey") +
  stat_summary(fun = median, geom = "point", shape = 15, size = 2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 1), size = 3) +
  geom_text(data = listSpeciesN, aes(y = ifelse(n>1, max_stocks + 2, max_stocks + 2), label = letters)) +
  # geom_point(aes(y = mean_pred)) +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/Cstocks_by_Depth_BOXPLOT_no_fitted.pdf")

#####


# # # # # # # # # # # # # # # # # # # # #
# Modelling MEADOW TYPE EFFECT ----
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

cstocks_data <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti) %>% 
  select(Meadow_type, depth, Source, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  

# Full linear model
m0 <- lm(log(cstocks) ~ Meadow_type*depth, data = cstocks_data)
step(m0) # interaction cannot be dropped.

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Meadow_type*depth, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Meadow_type*depth, 
                weights = varIdent(form = ~1|depth), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We don't need weights!

# Do we need random effects?
m0.lme <- lme(log(cstocks) ~ Meadow_type*depth, 
                random = ~1|Source,
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data)
anova(m0.gls, m0.lme) # The random effect Source is needed

# Model selection
m1.lme <- update(m0.lme, method = "ML")
m2.lme <- update(m1.lme, ~. -Meadow_type:depth)
anova(m1.lme, m2.lme) # The interaction can be dropped.

m3.lme <- update(m2.lme, ~. -Meadow_type)
anova(m2.lme, m3.lme) # Meadow type can't be dropped.

m4.lme <- update(m2.lme, ~. -depth)
anova(m2.lme, m4.lme) # Depth can't be dropped.

mfinal <- lme(log(cstocks) ~ Meadow_type + depth, 
              random = ~1|Source,
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
Anova(mfinal)
summary(mfinal)

# Model validation
mcheck(mfinal) # there are some very high residuals.
# min(resid(mfinal))

# Checking for extreme values. deleting residuals bigger 1.5
cstocks_data_clean <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logMeadowRes") %>% 
  add_predictions(model = mfinal, var = "logMeadowPred") %>% 
  filter(abs(logMeadowRes)<1.5) # 7 extreme values deleted (1%)


# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Meadow_type*depth, data = cstocks_data_clean)
m0.gls.w <- gls(log(cstocks) ~ Meadow_type*depth, 
                weights = varIdent(form = ~1|depth), 
                data = cstocks_data_clean)
anova(m0.gls, m0.gls.w) # We don't need weights!

# Do we need random effects?
m0.lme <- lme(log(cstocks) ~ Meadow_type*depth, 
                random = ~1|Source,
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_clean)
anova(m0.gls, m0.lme) # The random effect Source is needed

# Model selection
m1.lme <- update(m0.lme, method = "ML")
m2.lme <- update(m1.lme, ~. -Meadow_type:depth)
anova(m1.lme, m2.lme) # The interaction can be dropped.

m3.lme <- update(m2.lme, ~. -Meadow_type)
anova(m2.lme, m3.lme) # Meadow type can't be dropped.

m4.lme <- update(m2.lme, ~. -depth)
anova(m2.lme, m4.lme) # Depth can't be dropped.

mfinal <- lme(log(cstocks) ~ Meadow_type + depth, 
              random = ~1|Source,
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data_clean)

Anova(mfinal) # Meadow_type is not significant
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
# Chisq Df Pr(>Chisq)  
# Meadow_type 3.4857  1    0.06190 .
# depth       2.8507  1    0.09133 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(mfinal)


# Model validation
mcheck(mfinal) #  Underdispersion... tricky...

# Multiple comparisons with library(multcompar)
multcompar <- glht(mfinal, linfct=mcp(Meadow_type="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar))
# Significant multiple comparisons
summary_multcompar %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THERE DOESN'T SEEM TO BE A STRONG MEADOW_TYPE EFFECT, ALTHOUGH THE DATA SET HAS SOME UNDERDISPERSION  #
# WHICH MAKES THE RESULTS INFERENCE A BIT MORE TRICKY. BUT IF MEADOW_TYPE WAS REALLY IMPORTANT          #
# IT WOULD HAVE BEEN SIGNNIFICANT, AS IT HAPPENS FOR DEPTH.                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # #
# Plotting MEADOW TYPE EFFECT ----
# # # # # # # # # # # # # # # # # # # # #
cstocksFilteredMonoMulti2 <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti) %>% 
  filter(!is.na(cstocks)) %>% 
  group_by(Meadow_type) %>% 
  summarise(n = n(),
            mean_stocks = mean(cstocks),
            error_stocks = std.error(cstocks),
            max_stocks = max(cstocks)) %>% 
  left_join(letters_multcompar, by = 'Meadow_type')

cstocks_data_clean %>% 
  ggplot(aes(x = Meadow_type, y = cstocks)) +
  geom_boxplot(aes(fill = Meadow_type, colour = Meadow_type)) +
  stat_summary(fun=median, geom="point", shape = 15, size = 2) +
  scale_fill_manual(values = c("#31a354", "#b2df8a")) +
  scale_colour_manual(values = c("#31a354", "#b2df8a")) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = Meadow_type, y = max_stocks + 0.5, label = str_c("(", n, ")"))) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = Meadow_type, y = max_stocks + 1, label = letters)) +
  scale_x_discrete(labels = c("Monospecific", "Multispecific")) +
  xlab("Meadow type") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1~cm^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/Cstocks_meadowtype_by_depth_modelled.pdf")


# # # # # # # # # # # # # # # # # # # # #
# # Modelling AFFINITY EFFECT, with depth ----
# # # # # # # # # # # # # # # # # # # # #
# 
# # We create the affinity variable based on latitude
# cstocks_data <- cstocks_tidy %>% 
#   mutate(affinity = ifelse(abs(Latitude) <= 23, "tropical", ifelse(abs(Latitude) > 23 & Posi == "Posi", "Posi", "temperate")),
#          affinity = factor(affinity)) %>% 
#   filter(!is.na(affinity)) %>% 
#   select(affinity, depth, Posi, Source, cstocks, CoreID_unique) %>% 
#   filter_all(all_vars(!is.na(.))) 
# 
# # Full linear model
# m0 <- lm(log(cstocks) ~ affinity*depth, data = cstocks_data)
# step(m0) # interaction cannot be dropped.
# mcheck(m0)
# 
# # Do we need weights?
# m0.gls <- gls(log(cstocks) ~ affinity*depth, data = cstocks_data)
# m0.gls.w <- gls(log(cstocks) ~ affinity*depth, 
#                 weights = varIdent(form = ~1|affinity), 
#                 data = cstocks_data)
# anova(m0.gls, m0.gls.w) # We need weights for affinity!
# 
# # Do we need random effects?
# m0.w.lme <- lme(log(cstocks) ~ affinity*depth, 
#                 random = ~1|Source,
#                 weights = varIdent(form = ~1|affinity), 
#                 control = lmeControl(maxIter = 100, msMaxIter = 100),
#                 data = cstocks_data)
# anova(m0.gls.w, m0.w.lme) # The random effect Source is needed
# 
# # m0.w.lme2 <- lme(log(cstocks) ~ affinity*depth, 
# #                  random = ~1|Source/CoreID_unique,
# #                  weights = varIdent(form = ~1|affinity), 
# #                  control = lmeControl(maxIter = 100, msMaxIter = 100),
# #                  data = cstocks_data)
# # anova(m0.gls.w, m0.w.lme2) # The random effect CoreID_unique is needed
# # But CoreID_unique always messes the residuals. We won't include this random effect, in the end. 
# 
# # Model selection 
# m1.w.lme <- update(m0.w.lme, method = "ML")
# m2.w.lme <- update(m1.w.lme, ~. -affinity:depth)
# anova(m1.w.lme, m2.w.lme) # the interaction CAN be dropped.
# 
# m3.w.lme <- update(m2.w.lme, ~. -depth)
# anova(m2.w.lme, m3.w.lme) # depth CAN be dropped.
# 
# m4.w.lme <- update(m2.w.lme, ~. -affinity)
# anova(m2.w.lme, m4.w.lme) # affinity CANNOT be dropped.
# 
# mfinal <- lme(log(cstocks) ~ affinity + depth, 
#               random = ~1|Source,
#               weights = varIdent(form = ~1|affinity), 
#               control = lmeControl(maxIter = 100, msMaxIter = 100),
#               data = cstocks_data)
# Anova(mfinal)
# summary(mfinal)
# 
# # Model validation
# mcheck(mfinal) # not to bad
# # min(resid(mfinal))
# 
# # Checking for extreme values. deleting residuals bigger 1.5
# cstocks_data_clean <- cstocks_data %>% 
#   add_residuals(model = mfinal, var = "logAffinityRes") %>% 
#   add_predictions(model = mfinal, var = "logAffinityPred") %>% 
#   filter(abs(logAffinityRes)<1.5) # 9 extreme values deleted (1%)
# 
# 
# # Do we need weights?
# m0.gls <- gls(log(cstocks) ~ affinity*depth, data = cstocks_data_clean)
# m0.gls.w <- gls(log(cstocks) ~ affinity*depth, 
#                 weights = varIdent(form = ~1|affinity), 
#                 data = cstocks_data_clean)
# anova(m0.gls, m0.gls.w) # We need weights for affinity!
# 
# # Do we need random effects?
# m0.w.lme <- lme(log(cstocks) ~ affinity*depth, 
#                 random = ~1|Source,
#                 weights = varIdent(form = ~1|affinity), 
#                 control = lmeControl(maxIter = 100, msMaxIter = 100),
#                 data = cstocks_data_clean)
# anova(m0.gls.w, m0.w.lme) # The random effect Source is needed
# 
# # Model selection 
# m1.w.lme <- update(m0.w.lme, method = "ML")
# m2.w.lme <- update(m1.w.lme, ~. -affinity:depth)
# anova(m1.w.lme, m2.w.lme) # the interaction CAN be dropped.
# 
# m3.w.lme <- update(m2.w.lme, ~. -depth)
# anova(m2.w.lme, m3.w.lme) # depth CAN be dropped.
# 
# m4.w.lme <- update(m3.w.lme, ~. -affinity)
# anova(m2.w.lme, m4.w.lme) # affinity CANNOT be dropped.
# 
# mfinal <- lme(log(cstocks) ~ affinity, 
#               random = ~1|Source,
#               weights = varIdent(form = ~1|affinity), 
#               control = lmeControl(maxIter = 100, msMaxIter = 100),
#               data = cstocks_data_clean)
# Anova(mfinal)
# # Analysis of Deviance Table (Type II tests)
# # Response: log(cstocks)
# #            Chisq Df Pr(>Chisq)
# # affinity 122.73  2  < 2.2e-16 ***
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# summary(mfinal)
# 
# # Model validation
# mcheck(mfinal) # there are some very high residuals.
# # min(resid(mfinal))
# 
# cstocks_data_clean <- cstocks_data_clean %>% 
#   add_residuals(model = mfinal, var = "logAffinityRes") %>% 
#   add_predictions(model = mfinal, var = "logAffinityPred") 
# 
# 
# # MULTIPLE COMPARISONS! 
# # Initially not working because of unused levels in the factor Species (a legacy of the subset)
# length(coef(mfinal))
# 
# # If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
# length(levels(cstocks_data_clean$affinity))
# 
# # Multiple comparisons
# multcompar <- glht(mfinal, linfct=mcp(affinity="Tukey"))
# summary_multcompar <- broom::tidy(summary(multcompar))
# 
# letters_multcompar <- broom::tidy(cld(multcompar))
# 
# # Significant multiple comparisons
# summary_multcompar %>% 
#   filter(adj.p.value<0.05) %>% 
#   print(n = Inf)
# # lhs                    rhs estimate std.error statistic p.value
# # temperate - Posi         0  -1.19      0.107    -11.1     0    
# # tropical - Posi          0  -1.21      0.137     -8.80    0    
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # THERE REALLY ISN'T AN AFFINITY EFFECT, BECAUSE TEMPERATE-TROPICAL COMPARISON IS NON-SIGNIFICANT #
# # WE'RE SEEING A SPECIES EFFECT HERE!                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# 
# # # # # # # # # # # # # # # # # # # # #
# # Plotting AFFINITY EFFECT, with depth ----
# # # # # # # # # # # # # # # # # # # # #
# 
# cstocks_data_clean %>%
#   mutate(Pred = exp(logAffinityPred)) %>% 
#   group_by(affinity) %>% 
#   summarise(n = n(),
#             mean_stock = mean(cstocks),
#             median_stock = median(cstocks),
#             std.error = std.error(cstocks), 
#             max_stock = max(cstocks),
#             mean_pred = mean(Pred)) %>% 
#   left_join(letters_multcompar, by = "affinity") %>% 
#   ggplot(aes(x = affinity, y = mean_stock, ymin=mean_stock-std.error, ymax=mean_stock+std.error)) +
#   geom_bar(aes(fill = affinity), stat = "identity") +
#   scale_fill_manual(values = c("#1f78b4", "#31a354", "#b2df8a")) +
#   geom_errorbar(width=.3, position = position_dodge(preserve = "single", width = 0.9)) +
#   geom_point(aes(y = mean_pred)) +
#   geom_text(aes(y = mean_stock+std.error + 0.5, label = str_c("(", n, ")"))) +
#   geom_text(aes(y = mean_stock+std.error + 1, label = letters)) +
#   xlab("") +
#   scale_x_discrete(labels = c("Mediterranean", "Temperate", "Tropical")) +
#   ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
#   coord_flip() +
#   theme_bw() +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 14))
# # ggsave("Figs/Cstocks_tropicalVStemperate_by_depthBARPLOT_TUKEY.pdf")
# 
# 
# # # # # # # # # # # # # # # # # # # # #
# # BoxPlotting AFFINITY EFFECT, with depth ----
# # # # # # # # # # # # # # # # # # # # #
# 
# cstocks_data_box <- cstocks_data_clean %>% 
#   group_by(affinity) %>% 
#   summarise(n = n(),
#             mean_stock = mean(cstocks, na.rm = T),
#             median_stock = median(cstocks, na.rm = T),
#             std.error = std.error(cstocks, na.rm = T), 
#             max_stock = max(cstocks, na.rm = T)) %>% 
#   left_join(letters_multcompar, by = "affinity")
# 
# cstocks_data_clean %>% 
#   ggplot() +
#   geom_boxplot(aes(x = affinity, y = cstocks, fill = affinity)) +
#   scale_fill_manual(values = c("#1f78b4", "#31a354", "#b2df8a")) +
#   # geom_point(data = cstocks_data_box, aes(x = affinity, y = median_stock), shape = 15, size = 2) +
#   geom_text(data = cstocks_data_box, aes(x = affinity, y = max_stock + 2, label = str_c("(", n, ")"))) +
#   geom_text(data = cstocks_data_box, aes(x = affinity, y = max_stock + 4, label = letters)) +
#   xlab("") +
#   scale_x_discrete(labels = c("Mediterranean", "Temperate", "Tropical")) +
#   ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
#   coord_flip() +
#   theme_bw() +
#   theme(legend.position = "none",
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 14))
# # ggsave("Figs/Cstocks_tropicalVStemperate_by_depthBOXPLOT_TUKEY.pdf")
#####


# # # # # # # # # # # # # # # # # # # #
# Modelling LATITUDINAL EFFECT (ABSOLUTE LATITUDE), with depth no SPECIES ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data <- cstocks_tidy %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  filter(Meadow_type == "monospecific") %>% 
  select(Source, depth, Latitude, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  

# Full linear model
m0 <- lm(log(cstocks) ~ Latitude*depth, data = cstocks_data)
step(m0) # Interaction not dropped

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Latitude*depth, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Latitude*depth, 
              weights = varIdent(form = ~1|depth),
              data = cstocks_data)
anova(m0.gls, m0.gls.w) # We don't need weights!

# Do we need random effects?
m0.lme <- lme(log(cstocks) ~ Latitude*depth, 
                random = ~1|Source,
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data)
anova(m0.gls, m0.lme) # The random effect Source is needed

# Model selection 
m1.lme <- update(m0.lme, method = "ML")
m2.lme <- update(m1.lme, ~. -Latitude:depth)
anova(m1.lme, m2.lme) # Interaction CANNOT be dropped.

mfinal <- lme(log(cstocks) ~ Latitude*depth, 
              random = ~1|Source,
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
Anova(mfinal)

# Model validation
mcheck(mfinal) # quite ok
# min(resid(mfinal))

# Checking for extreme values. deleting residuals bigger 1.5
cstocks_data_clean <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred") %>% 
  filter(abs(logSpeciesRes)<1.5) # 10 extreme values deleted (1%)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Latitude*depth, data = cstocks_data_clean)
m0.gls.w <- gls(log(cstocks) ~ Latitude*depth, 
                weights = varIdent(form = ~1|depth), 
                data = cstocks_data_clean)
anova(m0.gls, m0.gls.w) # We DON'T need weights!

# Do we need random effects?
m0.lme <- lme(log(cstocks) ~ Latitude*depth, 
                random = ~1|Source,
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_clean)
anova(m0.gls, m0.lme) # The random effect Source is needed

# Model selection 
m1.lme <- update(m0.lme, method = "ML")
m2.lme <- update(m1.lme, ~. -Latitude:depth)
anova(m1.lme, m2.lme) # Interaction CANNOT be dropped.

mfinal <- lme(log(cstocks) ~ Latitude*depth, 
              random = ~1|Source,
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data_clean)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
#                 Chisq Df Pr(>Chisq)  
# Latitude       2.5897  1     0.1076  
# depth          0.0275  1     0.8684  
# Latitude:depth 3.5128  1     0.0609 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(mfinal)

# Model validation
mcheck(mfinal) #  GOOD FIT

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# POTENTIAL INTERACTION BETWEEN DEPTH AND LATITUDE ON CARBON STOCKS, WITH EFFECTS OF DEPTH BEING MODULATED #
# BY LATITUDE? MARGINALLY SIGNIFICANT... BUT WE SHOULD BEAR IN MIND WE HAVEN'T CONTROLLED BY SPECIES,      #
# SO WE SHOULD CHECK THE RESIDUALS OF THE SPECIES*DEPTH MODEL. SEE BELOW.                                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # #
# Modelling LATITUDINAL EFFECT USING THE RESIDUALS FROM THE MODEL WITH SPECIES AND DEPTH ----
# # # # # # # # # # # # # # # # # # # #


# # # # # #  # 
# IMPORTANT! #
# 1st you have to RUN SPECIES MODELS TO GET THE RESIDUALS! #
# # # # #  # #

cstocks_data <- cstocks_data_clean_Species %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  # filter(Meadow_type == "monospecific") %>% not needed, because this filter was already done in the species analysis
  filter_all(all_vars(!is.na(.)))  

# Full linear model
m0 <- lm(logSpeciesRes ~ Latitude*depth, data = cstocks_data)
step(m0) 

mcheck(m0) # +-OK!

# Let's get rid of biggest residuals
cstocks_data <- cstocks_data_clean_Species %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  filter(abs(logSpeciesRes)<1.5) %>% 
  # filter(Meadow_type == "monospecific") %>% not needed, because this filter was already done in the species analysis
  filter_all(all_vars(!is.na(.)))  

# Full linear model
m0 <- lm(logSpeciesRes ~ Latitude*depth, data = cstocks_data)
step(m0) # Latitude not dropped.

m1 <- lm(logSpeciesRes~Latitude+depth, data = cstocks_data)
Anova(m1)
# Anova Table (Type II tests)
# Response: logSpeciesRes
#             Sum Sq  Df F value  Pr(>F)  
# Latitude    0.949   1  3.6833 0.05539 .
# depth       0.667   1  2.5857 0.10831  
# Residuals 168.322 653  
summary(m1)
mcheck(m1) # OK!

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# LATITUDINAL EFFECT ON CARBON STOCKS IS ON THE LIMIT OF SIGNIFICANCE,  #
# ONCE YOU HAVE TAKEN SPECIES AND DEPTH INTO ACCOUNT!                   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # #
# Modelling Lat_Zone_Posi EFFECT, with depth ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data <- cstocks_tidy %>% 
  filter(Meadow_type == "monospecific") %>%
  mutate(Lat_Zone_Posi = factor(ifelse(Posi == "Posi", "Posi", Lat_Zone)),
         FIN_TYP_names = factor(recode(FIN_TYP,
                                       `0` = "Endorheic or Glaciated",
                                       `1` = "Small deltas",
                                       `2` = "Tidal systems",
                                       `3` = "Lagoons",
                                       `4` = "Fjords and fjaerds",
                                       `5` = "Large rivers",
                                       `6` = "Karst",
                                       `7` = "Arheic"))) %>% 
  # filter(Species != "Posidonia oceanica") %>% # To check the effect of Species without posidonia oceanica
  # select(Latitude, Species, Source, depth, cstocks, CoreID_unique) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  
# filter(Species != "Posidonia oceanica" | cstocks != 1.9)

# Full linear model
m0 <- lm(log(cstocks) ~ Lat_Zone_Posi*depth, data = cstocks_data)
step(m0) # interaction cannot be dropped.
mcheck(m0)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Lat_Zone_Posi+depth, data = cstocks_data) # The model with interaction is singular. So we drop interaction.
m0.gls.w <- gls(log(cstocks) ~ Lat_Zone_Posi+depth, 
                weights = varIdent(form = ~1|Lat_Zone_Posi), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We need weights for Lat_Zone_Posi!

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Lat_Zone_Posi + depth, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Lat_Zone_Posi), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data)
anova(m0.gls.w, m0.w.lme) # The random effect Source is needed

# Model selection 
m1.w.lme <- update(m0.w.lme, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m2.w.lme) # depth CAN be dropped.

m3.w.lme <- update(m2.w.lme, ~. -Lat_Zone_Posi)
anova(m2.w.lme, m3.w.lme) # Lat_Zone_Posi CANNOT be dropped.

mfinal <- lme(log(cstocks) ~ Lat_Zone_Posi, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Lat_Zone_Posi), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
Anova(mfinal)
summary(mfinal)

# Model validation
mcheck(mfinal) # not to bad
# min(resid(mfinal))

# Checking for extreme values. deleting residuals bigger 1.5
cstocks_data_clean <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logLat_Zone_PosiRes") %>% 
  add_predictions(model = mfinal, var = "logLat_Zone_PosiPred") %>% 
  filter(abs(logLat_Zone_PosiRes)<1.5) # 9 extreme values deleted (1%)

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Lat_Zone_Posi + depth, data = cstocks_data_clean)
m0.gls.w <- gls(log(cstocks) ~ Lat_Zone_Posi + depth, 
                weights = varIdent(form = ~1|Lat_Zone_Posi), 
                data = cstocks_data_clean)
anova(m0.gls, m0.gls.w) # We need weights for Lat_Zone_Posi!

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Lat_Zone_Posi + depth, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Lat_Zone_Posi), 
                control = lmeControl(maxIter = 100, msMaxIter = 100),
                data = cstocks_data_clean)
anova(m0.gls.w, m0.w.lme) # The random effect Source is needed

# Model selection 
m1.w.lme <- update(m0.w.lme, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m2.w.lme) # depth CAN be dropped.

m3.w.lme <- update(m2.w.lme, ~. -Lat_Zone_Posi)
anova(m2.w.lme, m3.w.lme) # Lat_Zone_Posi CANNOT be dropped.

mfinal <- lme(log(cstocks) ~ Lat_Zone_Posi, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Lat_Zone_Posi), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data_clean)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
#            Chisq Df Pr(>Chisq)
# affinity 127.36   3  < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(mfinal)

# Model validation
mcheck(mfinal) # there are some very high residuals.
# min(resid(mfinal))

cstocks_data_clean <- cstocks_data_clean %>% 
  add_residuals(model = mfinal, var = "logLat_Zone_PosiRes") %>% 
  add_predictions(model = mfinal, var = "logLat_Zone_PosiPred") 


# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(mfinal))

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(cstocks_data_clean$Lat_Zone_Posi))

# Multiple comparisons
multcompar <- glht(mfinal, linfct=mcp(Lat_Zone_Posi="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar))

# Significant multiple comparisons
summary_multcompar %>% 
  filter(adj.p.value<0.05) %>% 
  print(n = Inf)
# A tibble: 4 x 7
# term          contrast             null.value estimate std.error statistic adj.p.value
# <chr>         <chr>                     <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 Lat_Zone_Posi Posi - Polar                  0    2.30     0.662       3.47     0.00228
# 2 Lat_Zone_Posi Temperate - Posi              0   -1.31     0.120     -10.8      0      
# 3 Lat_Zone_Posi Tropical - Posi               0   -1.09     0.112      -9.69     0      
# 4 Lat_Zone_Posi Tropical - Temperate          0    0.217    0.0850      2.55     0.0416 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THERE APPEARS TO BE A SUBTLE LAT_ZONE_POSI EFFECT. APART FROM POSI TO THE REST,                 #
# THERE'S A JUST SIGNIFICANT DIFFERENCE BETWEEN TEMPERATE AND TROPICAL                            #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # #
# BARPlotting LAT_ZONE_POSI EFFECT    ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_clean %>%
  mutate(Pred = exp(logLat_Zone_PosiPred)) %>% 
  group_by(Lat_Zone_Posi) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks),
            median_stock = median(cstocks),
            std.error = std.error(cstocks), 
            max_stock = max(cstocks),
            mean_pred = mean(Pred)) %>% 
  left_join(letters_multcompar, by = "Lat_Zone_Posi") %>% 
  ggplot(aes(x = reorder(Lat_Zone_Posi, mean_stock, FUN = mean), y = mean_stock, ymin=mean_stock-std.error, ymax=mean_stock+std.error)) +
  geom_bar(aes(fill = Lat_Zone_Posi), stat = "identity") +
  scale_fill_manual(values = c("#1f78b4", "#31a354", "#b2df8a", "#f9d382")) +
  geom_errorbar(width=.3, position = position_dodge(preserve = "single", width = 0.9)) +
  # geom_point(aes(y = mean_pred)) +
  geom_text(aes(y = mean_stock+std.error + 0.5, label = str_c("(", n, ")"))) +
  geom_text(aes(y = mean_stock+std.error + 1, label = letters)) +
  xlab("") +
  scale_x_discrete(labels = c("Polar","Tropical", "Temperate", expression(italic("Posidonia oceanica")))) +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave("Figs/Cstocks_LAT_ZONE_POSI_BARPLOT_TUKEY.pdf")


# # # # # # # # # # # # # # # # # # # #
# BOXPlotting LAT_ZONE_POSI EFFECT, with depth ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data_clean %>% 
  group_by(Lat_Zone_Posi) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar, by = "Lat_Zone_Posi")

cstocks_data_clean %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Lat_Zone_Posi, cstocks, FUN = median), y = cstocks, fill = Lat_Zone_Posi, colour = Lat_Zone_Posi)) +
  scale_fill_manual(values = c("#1f78b4", "#31a354", "#b2df8a", "#f9d382")) +
  scale_colour_manual(values = c("#1f78b4", "#31a354", "#b2df8a", "#f9d382")) +
  geom_point(data = cstocks_data_box, aes(x = Lat_Zone_Posi, y = median_stock), shape = 15, size = 2) +
  geom_text(data = cstocks_data_box, aes(x = Lat_Zone_Posi, y = max_stock + 2, label = str_c("(", n, ")"))) +
  geom_text(data = cstocks_data_box, aes(x = Lat_Zone_Posi, y = max_stock + 4, label = letters)) +
  xlab("") +
  scale_x_discrete(labels = c("Polar","Tropical", "Temperate", expression(italic("Posidonia oceanica")))) +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave("Figs/Cstocks_LAT_ZONE_POSI_BOXPLOT_TUKEY.pdf")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Modelling FIN_TYP_names plus covariates ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Removing species with less than 10 samples
SpeciesMoreThan10Samples <- cstocks_data %>% 
  group_by(Species) %>% 
  summarise(n = n()) %>% 
  filter(n>=10) %>% 
  select(Species)
SpeciesMoreThan10Samples <- as.data.frame(SpeciesMoreThan10Samples)
SpeciesMoreThan10Samples <- SpeciesMoreThan10Samples$Species

cstocks_data <- cstocks_data %>% 
  filter(Species %in% SpeciesMoreThan10Samples)

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$Species <- droplevels(cstocks_data$Species) 

# # Option to remove cores with FIN_TYP_names too far away from FIN_TYP (high nearest_distance)
# cstocks_data <- cstocks_data %>%
  # filter(Nearest_distance < mean(Nearest_distance))
  # filter(Nearest_distance < quantile(Nearest_distance, probs = 0.95))

# Full linear model
m0 <- lm(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_data)

# Do we need SPECIES weights?
m0.gls <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We do need SPECIES weights!

# Do we need FIN_TYP_names weights?
m0.gls.w2 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|FIN_TYP_names), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w2) # We do need FIN_TYP_names weights!

# Do we need both type of weights?
m0.gls.w3 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 data = cstocks_data)
AIC(m0.gls, m0.gls.w, m0.gls.w2, m0.gls.w3)
anova(m0.gls.w, m0.gls.w3)
# YEP! BOTH WEIGHTS ARE NEEDED!

# Do we need also weights for depth?
m0.gls.w4 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names*depth), 
                 data = cstocks_data)
anova(m0.gls.w3, m0.gls.w4) # NOPE

# Do we need weights for Lat_Zone?
m0.gls.w5 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names*Lat_Zone), 
                 data = cstocks_data)
anova(m0.gls.w3, m0.gls.w5) # Nope!

# Therefore, we do need weights for SPECIES AND FIN_TYP_names

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data)
anova(m0.w.lme, m0.gls.w3) # Source is needed

# As said, although including coreID makes sense, it messes up model validation, and 
# thus we are not including it.

# Model selection 
m1.w.lme <- update(m0.w.lme, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -Species)
anova(m1.w.lme, m2.w.lme) # Species CANNOT be dropped.

m4.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m4.w.lme) # depth can be dropped... 

m5.w.lme <- update(m4.w.lme, ~. -FIN_TYP_names)
anova(m4.w.lme, m5.w.lme) # FIN_TYP_names CANNOT be dropped... 

m6.w.lme <- update(m4.w.lme, ~. -Lat_Zone)
anova(m4.w.lme, m6.w.lme) # Lat_Zone can be dropped

mfinal <- lme(log(cstocks) ~ Species + FIN_TYP_names, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species*FIN_TYP_names),
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data)
Anova(mfinal)
summary(mfinal)
mcheck(mfinal) # Residuals in qqplot are ok! (if we do not include coreID)


# For another analysis with latitude later
cstocks_data_clean_Species <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes")

# Checking for extreme values. deleting residuals bigger 1.5
cstocks_data_clean <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred") %>% 
  filter(abs(logSpeciesRes)<2) # 7 extreme values deleted (1%)

cstocks_data_clean$Species <- droplevels(cstocks_data_clean$Species) 

# Do we need SPECIES weights?
m0.gls <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, data = cstocks_data_clean)
m0.gls.w <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data_clean)
anova(m0.gls, m0.gls.w) # We do need SPECIES weights!

# Do we need FIN_TYP_names weights?
m0.gls.w2 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|FIN_TYP_names), 
                 data = cstocks_data_clean)
anova(m0.gls, m0.gls.w2) # We do need FIN_TYP_names weights!

# Do we need both type of weights?
m0.gls.w3 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 data = cstocks_data_clean)
AIC(m0.gls, m0.gls.w, m0.gls.w2, m0.gls.w3)
anova(m0.gls.w, m0.gls.w3)
# YEP! BOTH WEIGTHS ARE NEEDED!

# Do we need also weights for depth?
m0.gls.w4 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names*depth), 
                 data = cstocks_data_clean)
anova(m0.gls.w3, m0.gls.w4) # NOPE

# Do we need weights for Lat_Zone?
m0.gls.w5 <- gls(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names*Lat_Zone), 
                 data = cstocks_data_clean)
anova(m0.gls.w3, m0.gls.w5) # Nope!

# Therefore, we do need weights for SPECIES AND FIN_TYP_names

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data_clean)
anova(m0.w.lme, m0.gls.w3) # Source is needed


# Model selection 
m1.w.lme <- update(m0.w.lme, method = "ML")
m2.w.lme <- update(m1.w.lme, ~. -Species)
anova(m1.w.lme, m2.w.lme) # Species CANNOT be dropped.

m4.w.lme <- update(m1.w.lme, ~. -depth)
anova(m1.w.lme, m4.w.lme) # depth CAN be dropped.

m5.w.lme <- update(m4.w.lme, ~. -FIN_TYP_names)
anova(m4.w.lme, m5.w.lme) # FIN_TYP_names CANNOT be dropped.

m6.w.lme <- update(m4.w.lme, ~. -Lat_Zone)
anova(m4.w.lme, m6.w.lme) # Lat_Zone CAN be dropped.

mfinal <- lme(log(cstocks) ~ Species + FIN_TYP_names, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species*FIN_TYP_names),
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data_clean)

Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
#                 Chisq Df Pr(>Chisq)    
# Species       194.591 16  < 2.2e-16 ***
# FIN_TYP_names  58.793  6  7.914e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(mfinal)

cstocks_data_clean <- cstocks_data_clean %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred")

# Model validation
mcheck(mfinal) #  OK FIT

# Checking if we solved heterogeneity
boxplot(resid(mfinal) ~ Species, cstocks_data_clean)
boxplot(resid(mfinal, type = "pearson") ~ Species, cstocks_data_clean) # FANTASTIC! Check difference between both type of residuals.

boxplot(resid(mfinal) ~ FIN_TYP_names, cstocks_data_clean)
boxplot(resid(mfinal, type = "pearson") ~ FIN_TYP_names, cstocks_data_clean) # Good. Check difference between both type of residuals.


# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(mfinal))

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(cstocks_data_clean$Species))

# Multiple comparisons with library(multcompar) FOR SPECIES!
multcompar_species <- glht(mfinal, linfct=mcp(Species="Tukey"))
summary_multcompar_species <- broom::tidy(summary(multcompar_species))
letters_multcompar_species <- broom::tidy(cld(multcompar_species))
# Significant multiple comparisons
summary_multcompar_species %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)
# # A tibble: 15 x 7
# term    contrast                                     null.value estimate std.error statistic adj.p.value
# <chr>   <chr>                                             <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 Species Posidonia oceanica - Amphibolis antarctica            0    1.29     0.341       3.78  0.0106    
# 2 Species Zostera marina - Amphibolis antarctica                0    0.748    0.225       3.33  0.0495    
# 3 Species Halophila ovalis - Cymodocea serrulata                0   -0.820    0.201      -4.09  0.00336   
# 4 Species Halophila stipulacea - Cymodocea serrulata            0   -0.410    0.103      -3.98  0.00462   
# 5 Species Halophila ovalis - Enhalus acoroides                  0   -0.865    0.205      -4.22  0.00184   
# 6 Species Halophila stipulacea - Enhalus acoroides              0   -0.456    0.0930     -4.90  0.0000511 
# 7 Species Thalassodendron ciliatum - Enhalus acoroides          0   -0.235    0.0605     -3.88  0.00773   
# 8 Species Posidonia oceanica - Halophila ovalis                 0    1.56     0.341       4.57  0.000421  
# 9 Species Thalassia hemprichii - Halophila ovalis               0    0.861    0.201       4.29  0.00125   
# 10 Species Zostera marina - Halophila ovalis                    0    1.02     0.235       4.34  0.00105   
# 11 Species Posidonia oceanica - Halophila stipulacea            0    1.15     0.301       3.82  0.00913   
# 12 Species Thalassia hemprichii - Halophila stipulacea          0    0.451    0.0851      5.30  0.00000490
# 13 Species Zostera marina - Halophila stipulacea                0    0.609    0.159       3.83  0.00891   
# 14 Species Posidonia sinuosa - Posidonia oceanica               0   -1.15     0.316      -3.65  0.0171    
# 15 Species Zostera noltii - Posidonia oceanica                  0   -1.21     0.115     -10.5   0  

# # Multiple comparisons with library(multcompar) FOR depth!
# multcompar_depth <- glht(mfinal, linfct=mcp(depth="Tukey"))
# summary_multcompar_depth <- broom::tidy(summary(multcompar_depth))
# letters_multcompar_depth <- broom::tidy(cld(multcompar_depth))
# # Significant multiple comparisons
# summary_multcompar_depth %>%
#   filter(adj.p.value<0.05) %>%
#   print(n = Inf)
# # NOOOOO SIGNIFICANT DIFFERENCES FOR DEPTH!

# Multiple comparisons with library(multcompar) FOR FIN_TYP_NAMES!
multcompar_FIN_TYP_names <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
summary_multcompar_FIN_TYP_names <- broom::tidy(summary(multcompar_FIN_TYP_names))
letters_multcompar_FIN_TYP_names <- broom::tidy(cld(multcompar_FIN_TYP_names))
# Significant multiple comparisons
summary_multcompar_FIN_TYP_names %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)
# SIGNIFICANT DIFFERENCES FOR DEPTH!
# term          contrast                               null.value estimate std.error statistic  adj.p.value
# <chr>         <chr>                                       <dbl>    <dbl>     <dbl>     <dbl>        <dbl>
# 1 FIN_TYP_names Endorheic or Glaciated - Arheic                 0   -1.32     0.438      -3.01 0.0316     
# 2 FIN_TYP_names Fjords and fjaerds - Arheic                     0   -1.49     0.336      -4.44 0.000162   
# 3 FIN_TYP_names Lagoons - Endorheic or Glaciated                0    1.46     0.424       3.44 0.00740    
# 4 FIN_TYP_names Tidal systems - Endorheic or Glaciated          0    1.50     0.426       3.53 0.00557    
# 5 FIN_TYP_names Lagoons - Fjords and fjaerds                    0    1.64     0.323       5.07 0.00000356 
# 6 FIN_TYP_names Small deltas - Fjords and fjaerds               0    1.29     0.318       4.05 0.000800   
# 7 FIN_TYP_names Tidal systems - Fjords and fjaerds              0    1.68     0.312       5.38 0.000000811
# 8 FIN_TYP_names Small deltas - Lagoons                          0   -0.347    0.0824     -4.22 0.000318   
# 9 FIN_TYP_names Tidal systems - Small deltas                    0    0.389    0.107       3.64 0.00369      


# # # # # # # # # # # # # # # # # # # #
# BARPlotting FIN_TYP_names ----
# # # # # # # # # # # # # # # # # # # #
library(ggsci)
cstocks_data_clean %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(mean_cstock = mean(cstocks),
            std.error = std.error(cstocks),
            n = n()) %>%
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names") %>% 
  ggplot(aes(x = reorder(FIN_TYP_names, mean_cstock, FUN = mean), y = mean_cstock, ymin = mean_cstock - std.error, ymax = mean_cstock + std.error)) +
  geom_bar(aes(fill = FIN_TYP_names), stat = "identity") +
  scale_fill_d3() +
  geom_errorbar(width=.3, position = position_dodge(preserve = "single", width = 0.9)) +
  geom_text(aes(y = mean_cstock+std.error + 0.5, label = str_c("(", n, ")"))) +
  geom_text(aes(y = mean_cstock+std.error + 1, label = letters)) +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave("Figs/cstocks_by_geomorph_barplot.pdf")



# # # # # # # # # # # # # # # # # # # #
# BOXPlotting FIN_TYP_names ----
# # # # # # # # # # # # # # # # # # # #
library(ggsci)

a_fun <- function(x){
  return(data.frame(y = max(x), label = "a"))
}

cstocks_data_box <- cstocks_data_clean %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names")


cstocks_data_clean %>%
# cstocks_tidy %>%
  filter(!is.na(FIN_TYP_names)) %>% 
  filter(!is.na(cstocks)) %>% 
  # filter(Nearest_distance < mean(Nearest_distance)) %>%
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names") %>% 
  ggplot(aes(x = reorder(FIN_TYP_names, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = FIN_TYP_names, colour = FIN_TYP_names)) +
  scale_colour_d3() +
  scale_fill_d3() +
  stat_summary(fun = median, geom = "point", shape = 15, size = 2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 1)) +
  # stat_summary(fun.data = a_fun, geom = "text", position = position_nudge(y = 2)) +
  geom_text(data = cstocks_data_box, aes(x = FIN_TYP_names, y = max_stock + 2, label = letters)) +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(file = "Figs/cstocks_by_geomorph_allFINtyp.pdf")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# BOXPLOTS - C stocks by SPECIES by depth and FIN_TYP_NAMES with sample size coloured ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

cstocks_data_box <- cstocks_data_clean %>% 
  group_by(Species) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_species, by = "Species")


cstocks_data_clean %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  left_join(letters_multcompar_species, by = "Species") %>% 
  ggplot(aes(x = reorder(Species, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = Species, colour = Species)) +
  scale_colour_d3("category20") +
  scale_fill_d3("category20") +
  stat_summary(fun = median, geom = "point", shape = 15, size = 2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 1), size = 3) +
  geom_text(data = cstocks_data_box, aes(x = Species, y = max_stock + 2, label = letters)) +
  # geom_hline(aes(yintercept = median(cstocks_data_clean$cstocks)), lty = 2, colour = "darkgrey") +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"),)
# ggsave(file = "Figs/Cstocks_by_Species_BOXPLOT_ALL_FINTYP.pdf")
