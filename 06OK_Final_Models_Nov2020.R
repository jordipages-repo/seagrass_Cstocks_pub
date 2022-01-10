# # # # # #
# Analysing CSTOCKS with different predictors
# Jordi F. Pagès
# November 2020
# Working from home
# # # # # # 

library(nlme)
library(car)
library(modelr)
library(multcomp)
library(ggsci)
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

# Therefore, we do need weights for SPECIES and FIN_TYP_names

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data)
anova(m0.w.lme, m0.gls.w3) # Source is needed

# As said in script 6_Final_Models_according_to_Outline.R, 
# although including coreID makes sense, it messes up model validation, and thus we are not including it. 
# It doesn't change final results anyway.

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
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
#                 Chisq Df Pr(>Chisq)    
# Species       171.839 16  < 2.2e-16 ***
# FIN_TYP_names  56.098  6  2.781e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

mcheck(mfinal) # Residuals in qqplot are ok! (if we do not include coreID), just a couple of residuals > abs(2), 
                # which could be removed, in my opinion.

# For another analysis with latitude later
cstocks_data_clean_Species <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes")

# Checking for extreme values. deleting residuals bigger abs(2)
cstocks_data_clean <- cstocks_data %>% 
  add_residuals(model = mfinal, var = "logSpeciesRes") %>% 
  add_predictions(model = mfinal, var = "logSpeciesPred") %>% 
  filter(abs(logSpeciesRes)<2) # 2 extreme values deleted (<1%)

cstocks_data_clean$Species <- droplevels(cstocks_data_clean$Species) 

# Rechecking model selection with the edited data set (without extreme values)
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

# Therefore, we do need weights for SPECIES and FIN_TYP_names

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


# Multiple comparisons with library(multcompar) FOR depth!
# !! Before doing this, you have to temporarily add depth to mfinal
multcompar_depth <- glht(mfinal, linfct=mcp(depth="Tukey"))
summary_multcompar_depth <- broom::tidy(summary(multcompar_depth))
letters_multcompar_depth <- broom::tidy(cld(multcompar_depth))
# Significant multiple comparisons
summary_multcompar_depth %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)



# # # # # # # # # # # # # # # # # # # #
# BOXPlotting FIN_TYP_names ----
# # # # # # # # # # # # # # # # # # # #

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



# # # # # # # # # # # # # # # # # # # #
# BOXPLOTS - C stocks by SPECIES  ----
# # # # # # # # # # # # # # # # # # # #

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



# # # # # # # # # # # # # # # # # # # #
# # BOXPLOTS - C stocks by DEPTH   ----
# # # # # # # # # # # # # # # # # # # # 

cstocks_data_box <- cstocks_data_clean %>% 
  group_by(depth) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_depth, by = "depth")

cstocks_data_clean %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  left_join(letters_multcompar_depth, by = "depth") %>% 
  mutate(depth = factor(depth, levels = c("20-50 cm", "0-20 cm"))) %>%
  ggplot(aes(x = depth, y = cstocks)) +
  geom_boxplot(aes(fill = depth, colour = depth)) +
  scale_fill_manual(values = rev(c("#189ad3", "#005073"))) +
  scale_colour_manual(values = rev(c("#189ad3", "#005073"))) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 1), size = 3) +
  geom_text(data = cstocks_data_box, aes(x = depth, y = max_stock + 2, label = letters)) +
  # geom_hline(aes(yintercept = median(cstocks_data_clean$cstocks)), lty = 2, colour = "darkgrey") +
  xlab("") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/Cstocks_by_Depth_BOXPLOT_final.pdf")



# # # # # # # # # # # # # # # # # # # # #
# Modelling MEADOW TYPE EFFECT ----
# # # # # # # # # # # # # # # # # # # # #

# Meadow type cannot be assessed with the previous model, because there, we're assessing the effect of Species, and
# by definition, the effect of Species must be tested on MONOSPECIFIC meadows.

# Hence, we model Meadow_type below:

# We separate the variable Species into 5 species columns, to be able to know the different species present in multispecific meadows.
cstocksSpeciesSep <- cstocks_tidy %>% 
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

# Now we select those cores either from multispecific meadows, or from monospecific meadows, 
# made up of species that also can occurr in multispecific meadows (e.g. not posidonia oceanica)
cstocks_data <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti) %>% 
  select(CoreID_unique, Meadow_type, Latitude, FIN_TYP_names, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  # select(Meadow_type, depth, Source, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  

# To remove typologies that do not have multispecific meadows (e.g. polar zones)
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)


# In these models where we include meadow_type, I won't be including Lat_Zone_Posi, because they're collinear.


# Full linear model
m0 <- lm(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, data = cstocks_data)
step(m0) # Apparently, Meadow_type, FIN_TYP_names will be important.

# Do we need weights?
m0.gls <- gls(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
                weights = varIdent(form = ~1|FIN_TYP_names), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # Best with FIN_TYP_names weights.

m0.gls.w2 <- gls(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
                weights = varIdent(form = ~1|Meadow_type), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w2) # Best WITHOUT Meadow_type weights.

# Do we need random effects?
m0.lme <- lme(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
              random = ~1|Source,
              weights = varIdent(form = ~1|FIN_TYP_names), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
anova(m0.gls.w, m0.lme) # The random effect Source is needed


# Model selection
m1.lme <- update(m0.lme, method = "ML")
m2.lme <- update(m1.lme, ~. -Meadow_type)
anova(m1.lme, m2.lme) # We can't drop Meadow_type

m3.lme <- update(m1.lme, ~. -depth)
anova(m1.lme, m3.lme) #  We can't drop depth

m4.lme <- update(m1.lme, ~. -FIN_TYP_names)
anova(m1.lme, m4.lme) # FIN_TYP can't be dropped.

mfinal <- lme(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
              random = ~1|Source,
              weights = varIdent(form = ~1|FIN_TYP_names), 
              data = cstocks_data)
Anova(mfinal)
summary(mfinal)

# Model validation
mcheck(mfinal) # there are some weird residuals...  underdispersion? But not too bad
# min(resid(mfinal))


# Multiple comparisons with library(multcompar)
# For Meadow_type 
# multcompar <- glht(mfinal, linfct=mcp(Meadow_type="Tukey"))
# summary_multcompar <- broom::tidy(summary(multcompar))
# letters_multcompar <- broom::tidy(cld(multcompar))
# # Significant multiple comparisons
# summary_multcompar %>%
#   filter(adj.p.value<0.05) %>%
#   print(n = Inf)

# For FIN_typ_names 
# multcompar <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
# summary_multcompar <- broom::tidy(summary(multcompar))
# letters_multcompar <- broom::tidy(cld(multcompar))
# # Significant multiple comparisons
# summary_multcompar %>%
#   filter(adj.p.value<0.05) %>%
#   print(n = Inf)

# For depth 
multcompar <- glht(mfinal, linfct=mcp(depth="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar))
# Significant multiple comparisons
summary_multcompar %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)



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

cstocks_data %>% 
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
# ggsave("Figs/Cstocks_meadowtype_modelled_withFINandDepth.pdf")



# # # # # # # # # # # # # # # # # # # #
# Plotting FIN_TYPE_names EFFECT ----
# # # # # # # # # # # # # # # # # # # #
cstocksFilteredMonoMulti2 <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti) %>% 
  filter(!is.na(cstocks)) %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(n = n(),
            mean_stocks = mean(cstocks),
            error_stocks = std.error(cstocks),
            max_stocks = max(cstocks)) %>% 
  left_join(letters_multcompar, by = 'FIN_TYP_names')

cstocks_data %>% 
  ggplot(aes(x = reorder(FIN_TYP_names, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = FIN_TYP_names, colour = FIN_TYP_names)) +
  stat_summary(fun=median, geom="point", shape = 15, size = 2) +
  scale_colour_d3() +
  scale_fill_d3() +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = FIN_TYP_names, y = max_stocks + 0.5, label = str_c("(", n, ")"))) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = FIN_TYP_names, y = max_stocks + 1, label = letters)) +
  xlab("Meadow type") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1~cm^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/Cstocks_FIN_modelled_withMEadowTypeandDepth.pdf")


# # # # # # # # # # # # # # # # # # # # #
# Plotting Depth EFFECT ----
# # # # # # # # # # # # # # # # # # # # #
cstocksFilteredMonoMulti2 <- cstocksSpeciesSep %>% 
  filter(Sp1 %in% listSpeciesMulti) %>% 
  filter(!is.na(cstocks)) %>% 
  group_by(depth) %>% 
  summarise(n = n(),
            mean_stocks = mean(cstocks),
            error_stocks = std.error(cstocks),
            max_stocks = max(cstocks)) %>% 
  left_join(letters_multcompar, by = 'depth')

cstocks_data %>% 
  ggplot(aes(x = depth, y = cstocks)) +
  geom_boxplot(aes(fill = depth, colour = depth)) +
  stat_summary(fun=median, geom="point", shape = 15, size = 2) +
  scale_fill_manual(values = (c("#189ad3", "#005073"))) +
  scale_colour_manual(values = (c("#189ad3", "#005073"))) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = depth, y = max_stocks + 0.5, label = str_c("(", n, ")"))) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = depth, y = max_stocks + 1, label = letters)) +
  # scale_x_discrete(labels = c("Monospecific", "Multispecific")) +
  xlab("Depth") +
  ylab(bquote('Carbon density (Mg C' ~ha^-1~cm^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/Cstocks_depth_modelled_withFINandMeadowtype.pdf")
