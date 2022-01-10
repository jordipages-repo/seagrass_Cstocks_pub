# # # # # #
# Analysing CSTOCKS with stepwise linear models
# Jordi F. Pagès
# March 2021
# Working from CEAB and home
# # # # # # 

# # # # # # # # # # # # # # # # # # # 
# 0. Loading packages and data  ----
# # # # # # # # # # # # # # # # # # # 

library(nlme)
library(car)
library(modelr)
library(multcomp)
library(ggsci)

source("01_DataImport&Corrections_CarbonReview.R")
source("mcheck_function.R")
rm(list = c("cstocks_tidy", "cstocks")) # Because we want to use, only, the 20 cm stocks data, not carbon densities.

cstocks_data <- cstocks_tidy20Stocks %>% 
  filter(Meadow_type == "monospecific") %>%
  mutate(Lat_Zone_Posi = ifelse(Posi == "Posi", "Posi", Lat_Zone),
         Lat_Zone = factor(Lat_Zone),
         # CoreID_unique = factor(CoreID_unique),
         FIN_TYP_names = factor(recode(FIN_TYP,
                                       `0` = "Endorheic or Glaciated",
                                       `1` = "Small deltas",
                                       `2` = "Tidal systems",
                                       `3` = "Lagoons",
                                       `4` = "Fjords and fjaerds",
                                       `5` = "Large rivers",
                                       `6` = "Karst",
                                       `7` = "Arheic"))) %>% 
  # filter(Species != "Posidonia oceanica") %>% # To check the effect of Species without Posidonia oceanica
  # select(Latitude, Species, Source, depth, cstocks, CoreID_unique) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  
# filter(Species != "Posidonia oceanica" | cstocks != 1.9)

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Modelling the response 20cm-carbon stocks. ONLY CORES WITH 0-20 and 20-50 data, to check depth effects ----
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

# We now filter the data set to get only those cores for which we have both 0-20 and 20-50 data
coreID_both_depths <- cstocks_data %>% 
  group_by(CoreID_unique) %>% 
  summarise(n = n()) %>% 
  filter(n>1) %>% 
  select(CoreID_unique)
coreID_both_depths <- as.array(coreID_both_depths$CoreID_unique)
cstocks_data_both_depths <- cstocks_data %>% 
  filter(CoreID_unique %in% coreID_both_depths)

# To avoid singularity error below
cstocks_data_both_depths <- cstocks_data_both_depths %>% 
  filter(FIN_TYP_names != "Karst")

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data_both_depths$FIN_TYP_names <- droplevels(cstocks_data_both_depths$FIN_TYP_names)
cstocks_data_both_depths$Species <- droplevels(cstocks_data_both_depths$Species)  
cstocks_data_both_depths$CoreID_unique <- droplevels(cstocks_data_both_depths$CoreID_unique)

# Full linear model
m0 <- lm(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, data = cstocks_data_both_depths)
summary(m0)

# Do we need SPECIES weights?
m0.gls <- gls(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, data = cstocks_data_both_depths)
m0.gls.w <- gls(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data_both_depths)
anova(m0.gls, m0.gls.w) # We do need SPECIES weights!

# Do we need FIN_TYP_names weights?
m0.gls.w2 <- gls(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|FIN_TYP_names), 
                 data = cstocks_data_both_depths)
anova(m0.gls, m0.gls.w2) # We do need FIN_TYP_names weights!

# Do we need both type of weights?
m0.gls.w3 <- gls(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 data = cstocks_data_both_depths)
AIC(m0.gls, m0.gls.w, m0.gls.w2, m0.gls.w3)
anova(m0.gls.w, m0.gls.w3)
# No need to add both. Species is enough.

# Do we need also weights for depth?
m0.gls.w4 <- gls(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*depth), 
                 data = cstocks_data_both_depths)
anova(m0.gls.w, m0.gls.w4) # NOPE

# Do we need weights for Lat_Zone?
m0.gls.w5 <- gls(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*Lat_Zone), 
                 data = cstocks_data_both_depths)
anova(m0.gls.w, m0.gls.w5) # Nope!

# Therefore, we need weights for SPECIES only

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data_both_depths)
anova(m0.w.lme, m0.gls.w) # Source is needed as random intercept

m0.w.lme1 <- lme(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|CoreID_unique,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data_both_depths)
anova(m0.w.lme1, m0.gls.w) # CoreID_unique is needed as random intercept

m0.w.lme2 <- lme(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source/CoreID_unique,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data_both_depths)
anova(m0.w.lme, m0.w.lme2) # Nested model is better than any single random effect.

# m0.w.lme3 <- lme(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone,
#                  random = ~1 + CoreID_unique|Source,
#                  weights = varIdent(form = ~1|Species),
#                  control = lmeControl(maxIter = 1000, msMaxIter = 1000),
#                  data = cstocks_data_both_depths)
# anova(m0.w.lme, m0.w.lme3) # Nested model is better than random intercept and slope.

AIC(m0.w.lme, m0.w.lme1, m0.w.lme2)
#           df      AIC
# m0.w.lme  51 662.2859
# m0.w.lme1 51 611.6578
# m0.w.lme2 52 596.2468  This is the best model. The one with nested random effects.
anova(m0.w.lme1, m0.w.lme2) # And log-likelihood ratio tests confirms this.

# Model selection 
m1.w.lme.full <- update(m0.w.lme2, method = "ML")
m1.w.lme.sub <- update(m1.w.lme.full, .~. -Species:depth)
anova(m1.w.lme.full, m1.w.lme.sub) # Species:depth interaction CANNOT be dropped.

m1.w.lme.sub2 <- update(m1.w.lme.full, ~. -FIN_TYP_names)
anova(m1.w.lme.full, m1.w.lme.sub2) # FIN_TYP_names CAN be dropped... 

m2.w.lme.full <- m1.w.lme.sub2
m2.w.lme.sub <- update(m2.w.lme.full, ~. -Lat_Zone)
anova(m2.w.lme.full, m2.w.lme.sub) # Lat_Zone CAN be dropped... 

mfinal <- lme(log(cstocks) ~ Species*depth, 
              random = ~1|Source/CoreID_unique,
              weights = varIdent(form = ~1|Species),
              method = "REML",
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data_both_depths)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# 
# Response: log(cstocks)
#                Chisq Df Pr(>Chisq)    
# Species       74.202 14  3.314e-10 ***
# depth          8.360  1   0.003836 ** 
# Species:depth 32.534 14   0.003362 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(mfinal)

# Checking if we solved heterogeneity
boxplot(resid(mfinal) ~ Species, cstocks_data_both_depths)
boxplot(resid(mfinal, type = "pearson") ~ Species, cstocks_data_both_depths) # FANTASTIC! Check difference between both type of residuals.

# Validation
mcheck(mfinal) # Residuals not GOOD! OBVIOUS OVERDISPERSION! But this can happen with random effects. 
# See https://stats.stackexchange.com/questions/351352/adding-an-observation-level-random-term-messes-up-residuals-vs-fitted-plot-why

# Fitted versus response (MUCH BETTER WITH this nested model than with the random intercept!!!)
plot(fitted(mfinal) ~ log(cstocks_data_both_depths$cstocks), col="darkgrey", 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")


# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(mfinal))

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(cstocks_data_both_depths$Species))

# Multiple comparisons with library(multcompar) FOR SPECIES!
multcompar_species <- glht(mfinal, linfct=mcp(Species="Tukey"))
summary_multcompar_species <- broom::tidy(summary(multcompar_species))
letters_multcompar_species <- broom::tidy(cld(multcompar_species))
# Significant multiple comparisons
summary_multcompar_species %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf) #%>% 
  # write_csv(file = "Figs/Final_Figures_March2021/Only_cores_both_depths/multcomp_Species.csv")
# # A tibble: 15 x 7
#   term    contrast                                        null.value estimate std.error statistic adj.p.value
# 1 Species Posidonia oceanica - Cymodocea rotundata                0     1.67     0.427      3.92  0.00612   
# 2 Species Posidonia oceanica - Cymodocea serrulata                0     1.94     0.415      4.67  0.000184  
# 3 Species Posidonia oceanica - Enhalus acoroides                  0     1.72     0.349      4.93  0.0000385 
# 4 Species Posidonia oceanica - Halodule uninervis                 0     1.76     0.360      4.88  0.000105  
# 5 Species Posidonia oceanica - Halodule wrightii                  0     1.84     0.509      3.62  0.0203    
# 6 Species Posidonia australis - Halophila ovalis                  0     1.09     0.234      4.66  0.000343  
# 7 Species Posidonia oceanica - Halophila ovalis                   0     2.22     0.408      5.43  0.00000745
# 8 Species Posidonia oceanica - Halophila stipulacea               0     1.99     0.366      5.44  0.0000190 
# 9 Species Thalassodendron ciliatum - Posidonia australis          0    -1.39     0.413     -3.36  0.0468    
# 10 Species Posidonia sinuosa - Posidonia oceanica                  0    -1.30     0.382     -3.41  0.0406    
# 11 Species Thalassia hemprichii - Posidonia oceanica               0    -1.62     0.353     -4.60  0.000284  
# 12 Species Thalassodendron ciliatum - Posidonia oceanica           0    -2.51     0.474     -5.30  0.0000286 
# 13 Species Zostera marina - Posidonia oceanica                     0    -1.70     0.318     -5.36  0.00000550 


# Multiple comparisons with library(multcompar) FOR FIN_TYP_NAMES... we have to include it for a second.
multcompar_FIN_TYP_names <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
summary_multcompar_FIN_TYP_names <- broom::tidy(summary(multcompar_FIN_TYP_names))
letters_multcompar_FIN_TYP_names <- broom::tidy(cld(multcompar_FIN_TYP_names))
# Significant multiple comparisons
summary_multcompar_FIN_TYP_names %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)
# NO SIGNIFICANT MULTIPLE COMPARISONS.


# Multiple comparisons with library(multcompar) FOR depth!
multcompar_depth <- glht(mfinal, linfct=mcp(depth="Tukey"))
summary_multcompar_depth <- broom::tidy(summary(multcompar_depth))
letters_multcompar_depth <- broom::tidy(cld(multcompar_depth))
# Significant multiple comparisons
summary_multcompar_depth %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)
# NO SIGNIFICANT MULTIPLE COMPARISONS.


# Multiple comparisons with depth*Species with interaction
# tmp <- expand.grid(Species = unique(cstocks_data_both_depths$Species),
#                    depth = unique(cstocks_data_both_depths$depth))
# X <- model.matrix(~ depth * Species, data = tmp)
# glht(mfinal, linfct = X)
# Tukey <- contrMat(table(cstocks_data_both_depths$Species), "Tukey")
# K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
# rownames(K1) <- paste(levels(cstocks_data_both_depths$depth)[1], rownames(K1), sep = ":")
# K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
# rownames(K2) <- paste(levels(cstocks_data_both_depths$depth)[2], rownames(K2), sep = ":")
# K <- rbind(K1, K2)
# colnames(K) <- c(colnames(Tukey), colnames(Tukey))
# summary_multcompar_interaction <- broom::tidy(summary(glht(mfinal, linfct = K %*% X)))
# summary_multcompar_interaction %>%
#   filter(adj.p.value<0.05) %>%
#   print(n = Inf)


# The GOOD way of doing Multiple comparisons with depth*Species with interaction
# CAUTION! This piece of code takes ~10 min (or more) to run
cstocks_data_both_depths$Interaction <- interaction(cstocks_data_both_depths$Species, cstocks_data_both_depths$depth, drop=T)
mfinalI <- lme(log(cstocks) ~ Interaction - 1,
               random = ~ 1|Source/CoreID_unique,
               weights = varIdent(form = ~1|Species),
               method = "REML",
               control = lmeControl(maxIter = 1000, msMaxIter = 1000),
               data = cstocks_data_both_depths)
summary_multcompar_interaction <- broom::tidy(summary(glht(mfinalI, linfct=mcp(Interaction = "Tukey"))))
summary_multcompar_interaction %>%
  # filter(adj.p.value<0.05) %>%
  # print(n = Inf) #%>% 
  # write_csv(file = "Figs/Final_Figures_March2021/Only_cores_both_depths/multcomp_interaction_terms2.csv")



# # # # # # # # # # # # # # # # # # # #
# 1.1 BOXPlotting FIN_TYP_names ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data_both_depths %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names")


cstocks_data_both_depths %>%
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
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 20)) +
  # stat_summary(fun.data = a_fun, geom = "text", position = position_nudge(y = 2)) +
  geom_text(data = cstocks_data_box, aes(x = FIN_TYP_names, y = max_stock + 40, label = letters)) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(file = "Figs/Final_Figures_March2021/Only_cores_both_depths/cstocks_vs_FINtype.pdf")



# # # # # # # # # # # # # # # # # # # #
# 1.2 BOXPLOTS - C stocks by SPECIES ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data_both_depths %>% 
  group_by(Species) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_species, by = "Species")

cstocks_data_both_depths %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  left_join(letters_multcompar_species, by = "Species") %>% 
  ggplot(aes(x = reorder(Species, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = Species, colour = Species)) +
  scale_colour_d3("category20") +
  scale_fill_d3("category20") +
  stat_summary(fun = median, geom = "point", shape = 15, size = 2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 20), size = 3) +
  geom_text(data = cstocks_data_box, aes(x = Species, y = max_stock + 40, label = letters)) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  # facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"))
# ggsave(file = "Figs/Final_Figures_March2021/Only_cores_both_depths/cstocks_by_Species.pdf")



# # # # # # # # # # # # # # # # # # #
# 1.3.1 BOXPLOTS - C stocks by DEPTH ----
# # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data_both_depths %>% 
  group_by(depth) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_depth, by = "depth")

cstocks_data_both_depths %>%
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
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 20), size = 3) +
  geom_text(data = cstocks_data_box, aes(x = depth, y = max_stock + 40, label = letters)) +
  # geom_hline(aes(yintercept = median(cstocks_data_clean$cstocks)), lty = 2, colour = "darkgrey") +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/Final_Figures_March2021/Only_cores_both_depths/cstocks_by_Depth.pdf")


# # # # # # # # # # # # # # # # # # # # # # # # #
# 1.3.2 BOXPLOTS - C stocks by SPECIES*DEPTH ----
# # # # # # # # # # # # # # # # # # # # # # # # # 

cstocks_data_both_depths %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  mutate(depth = factor(depth, levels = c("20-50 cm", "0-20 cm"))) %>%
  ggplot(aes(x = reorder(Species, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = depth, colour = depth)) +
  scale_fill_manual(values = rev(c("#189ad3", "#005073"))) +
  scale_colour_manual(values = rev(c("#189ad3", "#005073"))) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1, aes(group = depth), position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = n_fun, geom = "text", size = 3, aes(y = cstocks + 10, group = depth), position = position_dodge(width = 0.75)) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.85,0.15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"))
# ggsave(filename = "Figs/Final_Figures_March2021/Only_cores_both_depths/cstocks_by_SpeciesXdepth.pdf")


# # # # # # # # # # # # # # # # # # # # # # # # #
# 1.3.3 BOXPLOTS - FACET C stocks by SPECIES*DEPTH ----
# # # # # # # # # # # # # # # # # # # # # # # # # 

cstocks_data_both_depths %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  mutate(depth = factor(depth, levels = rev(c("20-50 cm", "0-20 cm")))) %>%
  ggplot(aes(x = reorder(Species, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = depth, colour = depth)) +
  scale_fill_manual(values = c("#189ad3", "#005073")) +
  scale_colour_manual(values = c("#189ad3", "#005073")) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1, aes(group = depth), position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = n_fun, geom = "text", size = 3, aes(y = cstocks + 20, group = depth), position = position_dodge(width = 0.75)) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  coord_flip() +
  facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = c(0.85,0.15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"))
# ggsave(filename = "Figs/Final_Figures_March2021/Only_cores_both_depths/cstocks_by_SpeciesXdepthFACET.pdf")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1.4. Trying GLMM with Gamma distribution on the "both depths" data set ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# The BEST DISTRIBUTION FOR OVERDISPERSED DATA is the Gamma distribution.

# but let's check the mean vs. variance of the response variable
mean(cstocks_data_both_depths$cstocks)
var(cstocks_data_both_depths$cstocks)
hist(log(cstocks_data_both_depths$cstocks))
shapiro.test(log(cstocks_data_both_depths$cstocks))

# Let's try a gamma glmm
detach("package:nlme", unload=TRUE)
library(lme4)

m1 <- glmer(cstocks ~ Species*depth + FIN_TYP_names + Lat_Zone + (1|Source/CoreID_unique), 
            family = Gamma(link = "inverse"),
            data = cstocks_data_both_depths)

# Let's first try to solve the lme4 convergence warnings I get
# Some of the text below taken from https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
# In glmer(), large differences in the scales of parameters often lead to problems 
# (especially with calculating standard deviations of fixed effects)
# So, if I had continuous variables I would rescale them. However, all of my variables are categorical.
# If the fit is singular or near-singular, there might be a higher chance of a false positive 
# a higher chance that the model has actually misconverged (because the optimization problem is difficult on the boundary);
# and a reasonable argument that the random effects model should be simplified.
# The definition of singularity is that some of the constrained parameters of the random effects 
# theta parameters are on the boundary (equal to zero, or very very close to zero, say <10−6):
# tt <- getME(m1,"theta")
# ll <- getME(m1,"lower")
# min(tt[ll==0])
# # 0.02483161 NOT A PROBLEM IN THIS CASE
# # Try restarting from previous fit
# ss <- getME(m1,c("theta","fixef")) # it gets the fit from previous run.
# m2 <- update(m1, start = ss, control = glmerControl(optCtrl=list(maxfun=2e4))) # and use it as a starting point for next run.
# # Got rid of one of the warnings.
# # Let's try a different optimiser
# m3 <- update(m1,start = ss, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

# Anyway, let's check which link function is best
m2 <- update(m1, family = Gamma(link = "identity"))
anova(m1, m2) # Identity link is better.


# Checking random structure of the model
m1 <- glmer(cstocks ~ Species*depth + FIN_TYP_names + Lat_Zone + (1|Source/CoreID_unique), 
            family = Gamma(link = "identity"),
            data = cstocks_data_both_depths)
m2 <- glmer(cstocks ~ Species*depth + FIN_TYP_names + Lat_Zone + (1|Source), 
            family = Gamma(link = "identity"),
            data = cstocks_data_both_depths)
m3 <- glmer(cstocks ~ Species*depth + FIN_TYP_names + Lat_Zone + (1|CoreID_unique), 
            family = Gamma(link = "identity"),
            data = cstocks_data_both_depths)
m4 <- glmer(cstocks ~ Species*depth + FIN_TYP_names + Lat_Zone + (1|Source) + (1|CoreID_unique), 
            family = Gamma(link = "identity"),
            data = cstocks_data_both_depths)
AIC(m1, m2, m3, m4)
#    df  AIC
# m1 38 2547.417
# m2 37 2645.030
# m3 37 2544.312 # Best model
# m4 38 2547.417
anova(m1, m3) # And the difference with is not significant. We should use the simplest, m3, which also has highest AIC.
# No differences between the nested model and the one with 2 random effects (not nested), same model in fact 
# see an explanation of why m1 == m4 in Robert Long's answer to this question:
# https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified

# OK, so there is no doubt that the best random structure for our model is the one with just CoreID_unique.
# Now let's select the fixed effects by dropping predictors one by one.
m1.full <- glmer(cstocks ~ Species*depth + FIN_TYP_names + Lat_Zone + (1|CoreID_unique), 
            family = Gamma(link = "identity"),
            data = cstocks_data_both_depths)
m1.sub <- update(m1.full, .~. -FIN_TYP_names)
anova(m1.full, m1.sub) # We can drop FIN_TYP_names

m2.sub <- update(m1.sub, .~. -Lat_Zone)
anova(m1.sub, m2.sub) # We can drop Lat_Zone

m3.sub <- update(m2.sub, .~. -Species:depth)
anova(m3.sub, m2.sub) # MOdels are different. Stop dropping. 

# Best model
mfinal <- glmer(cstocks ~ Species*depth + (1|CoreID_unique), 
                family = Gamma(link = "identity"),
                data = cstocks_data_both_depths)
ss <- getME(mfinal,c("theta","fixef")) # it gets the fit from previous run.
mfinalOK <- update(mfinal, start = ss, control = glmerControl(optCtrl=list(maxfun=2e4))) # and use it as a starting point for next run.
Anova(mfinalOK)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: cstocks
#                 Chisq Df   Pr(>Chisq)    
# Species       249.2268 14  < 2.2e-16 ***
# depth           3.7275  1   0.053524 .  
# Species:depth  31.4815 14   0.004744 ** 
mcheck(mfinalOK)
boxplot(resid(mfinalOK) ~ Species, cstocks_data_both_depths)


#### I get convergence warnings every time I run these models... so I think I prefer sticking to the lme models.



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 2. Modelling the response 20cm-carbon stocks. ALL CORES ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Clear R memory!!!
# Then run section 0. 
# And then start here.

# Removing species with less than 10 samples
SpeciesMoreThan10Samples <- cstocks_data %>% 
  group_by(Species) %>% 
  summarise(n = n()) %>% 
  filter(n>=10) %>% 
  select(Species)
SpeciesMoreThan10Samples <- as.data.frame(SpeciesMoreThan10Samples)
SpeciesMoreThan10Samples <- SpeciesMoreThan10Samples$Species
# SpeciesMoreThan10Samples <- SpeciesMoreThan10Samples[-10]
cstocks_data <- cstocks_data %>% 
  filter(Species %in% SpeciesMoreThan10Samples)

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$Species <- droplevels(cstocks_data$Species) 

# WE HAVE ALREADY ESTABLISHED THE MARGINAL EFFECT OF DEPTH IN DETERMINING CSTOCKS. 
# SINCE WE ARE NOW NOT INTERESTED IN DEPTH EFFECTS PER SE, BUT WE WANT TO CONTROL ITS 
# POTENTIAL MARGINAL INFLUENCE, WE WILL ADD IT AS A RANDOM (NOT FIXED) EFFECT. 
# ALSO, IT'S DIFFICULT TO ADD IT AS FIXED, BECAUSE NOT ALL SPECIES HAVE ALL DEPTHS. 
# E.G. ZOSTERA NOLTII DOESN'T HAVE 20-50CM DATA.

# Full linear model
m0 <- lm(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, data = cstocks_data)
summary(m0)
Anova(m0)

# Do we need SPECIES weights?
m0.gls <- gls(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We do need SPECIES weights!

# Do we need FIN_TYP_names weights?
m0.gls.w2 <- gls(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|FIN_TYP_names), 
                 data = cstocks_data)
anova(m0.gls, m0.gls.w2) # We do need FIN_TYP_names weights!

# Do we need both type of weights?
m0.gls.w3 <- gls(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 data = cstocks_data)
AIC(m0.gls, m0.gls.w, m0.gls.w2, m0.gls.w3)
anova(m0.gls.w, m0.gls.w3)
# We need both weigths. Best model is m0.gls.w3

# Do we need weights for Lat_Zone?
m0.gls.w5 <- gls(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                 weights = varIdent(form = ~1|Species*Lat_Zone), 
                 data = cstocks_data)
anova(m0.gls.w, m0.gls.w5) # Nope!

# Therefore, we need weights for SPECIES and FIN_TYP_names 

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data)
anova(m0.w.lme, m0.gls.w3) # Source is needed as random intercept

m0.w.lme1 <- lme(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                 random = ~1|CoreID_unique,
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                 data = cstocks_data)
anova(m0.w.lme1, m0.gls.w3) # CoreID_unique is needed as random intercept

m0.w.lme2 <- lme(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                 random = ~1|Source/CoreID_unique,
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                 data = cstocks_data)
anova(m0.w.lme, m0.w.lme2) # Nested model is better than any single random effect.

m0.w.lme3 <- lme(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone, 
                 random = ~1|Source/CoreID_unique/depth,
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                 data = cstocks_data)
anova(m0.w.lme, m0.w.lme3) # Nested + nested model is better than a single random effect.

AIC(m0.w.lme, m0.w.lme1, m0.w.lme2, m0.w.lme3)
#           df      AIC
# m0.w.lme  67 1156.137
# m0.w.lme1 67 1225.750
# m0.w.lme2 68 1115.632 # this appears to be the best model
# m0.w.lme3 69 1117.632
anova(m0.w.lme2, m0.w.lme3) # And log-likelihood ratio tests confirms the model with Source/CoreID_unique as the best.

# Trying some other options for random structure with lme4
detach("package:nlme", unload=TRUE)
library(lme4)
m1 <- lmer(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + (1|Source), data = cstocks_data)
m2 <- lmer(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + (1|CoreID_unique), data = cstocks_data)
m2.bis <- lmer(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + (1|depth), data = cstocks_data)
m3 <- lmer(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + (1|Source/CoreID_unique), data = cstocks_data)
m4 <- lmer(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + (1|CoreID_unique) + (1|depth), data = cstocks_data)
m5 <- lmer(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + (1|Source/CoreID_unique) + (1|depth), data = cstocks_data)
AIC(m1, m2, m2.bis, m3, m4, m5)
anova(m3, m5)
Anova(m3)
Anova(m5)
# Again the best selected random structure is the one with Source/CoreID_unique 
# (depth can be dropped, but 0.07... in any case, same results with Anova(model))

# Model selection 
m1.w.lme.full <- update(m0.w.lme2, method = "ML")
m1.w.lme.sub <- update(m1.w.lme.full, .~. -Lat_Zone)
anova(m1.w.lme.full, m1.w.lme.sub) # Lat_Zone CAN be dropped.

m2.w.lme.full <- m1.w.lme.sub
m2.w.lme.sub <- update(m2.w.lme.full, ~. -FIN_TYP_names)
anova(m2.w.lme.full, m2.w.lme.sub) # FIN_TYP_names CANNOT be dropped... 

mfinal <- lme(log(cstocks) ~ Species + FIN_TYP_names, 
              random = ~1|Source/CoreID_unique,
              weights = varIdent(form = ~1|Species*FIN_TYP_names),
              method = "REML",
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# 
# Response: log(cstocks)
#                Chisq Df Pr(>Chisq)    
# Species       105.636 16  3.012e-15 ***
# FIN_TYP_names  36.424  6  2.280e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(mfinal)

mfinal_lme4 <- lme4::lmer(log(cstocks) ~ Species + FIN_TYP_names + (1|Source/CoreID_unique), data = cstocks_data)
icc_specs(mfinal_lme4) %>%
  mutate_if(is.numeric, round, 2)
plot_variance(mfinal_lme4)

# Checking if we solved heterogeneity
boxplot(resid(mfinal) ~ Species, cstocks_data)
boxplot(resid(mfinal, type = "pearson") ~ Species, cstocks_data) # FANTASTIC! Check difference between both type of residuals.

boxplot(resid(mfinal) ~ FIN_TYP_names, cstocks_data)
boxplot(resid(mfinal, type = "pearson") ~ FIN_TYP_names, cstocks_data) # FANTASTIC! Check difference between both type of residuals.

# Validation
mcheck(mfinal) # Residuals not GOOD! OBVIOUS OVERDISPERSION! But this can happen with random effects. 
# See https://stats.stackexchange.com/questions/351352/adding-an-observation-level-random-term-messes-up-residuals-vs-fitted-plot-why

# Fitted versus response (MUCH BETTER WITH this nested model than with the random intercept!!!)
plot(fitted(mfinal) ~ log(cstocks_data$cstocks), col="darkgrey", 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")


# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(mfinal))

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(cstocks_data$Species))

# Multiple comparisons with library(multcompar) FOR SPECIES!
multcompar_species <- glht(mfinal, linfct=mcp(Species="Tukey"))
summary_multcompar_species <- broom::tidy(summary(multcompar_species))
letters_multcompar_species <- broom::tidy(cld(multcompar_species, level = 0.06))
# Significant multiple comparisons
summary_multcompar_species %>%
  filter(adj.p.value<0.06) %>%
  print(n = Inf)# %>%
  # write_csv(file = "Figs/Final_Figures_March2021/All_cores/multcomp_Species.csv")
# A tibble: 4 x 7
# term    contrast                                      null.value estimate std.error statistic adj.p.value
# <chr>   <chr>                                              <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 Species Posidonia oceanica - Halophila ovalis                  0     1.22     0.327      3.71    1.55e- 2
# 2 Species Posidonia oceanica - Halophila stipulacea              0     1.09     0.320      3.41    4.45e- 2
# 3 Species Thalassodendron ciliatum - Posidonia oceanica          0    -1.09     0.326     -3.34    5.44e- 2
# 4 Species Zostera noltii - Posidonia oceanica                    0    -1.26     0.145     -8.70    3.33e-16


# Multiple comparisons with library(multcompar) FOR FIN_TYP_NAMES
multcompar_FIN_TYP_names <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
summary_multcompar_FIN_TYP_names <- broom::tidy(summary(multcompar_FIN_TYP_names))
letters_multcompar_FIN_TYP_names <- broom::tidy(cld(multcompar_FIN_TYP_names, level = 0.05))
# Significant multiple comparisons
summary_multcompar_FIN_TYP_names %>%
  filter(adj.p.value<0.06) %>%
  print(n = Inf)
# term          contrast                               null.value estimate std.error statistic adj.p.value
# <chr>         <chr>                                       <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 FIN_TYP_names Fjords and fjaerds - Arheic                     0    -1.39     0.348     -3.98  0.00104   
# 2 FIN_TYP_names Lagoons - Endorheic or Glaciated                0     1.30     0.440      2.95  0.0389    
# 3 FIN_TYP_names Tidal systems - Endorheic or Glaciated          0     1.43     0.437      3.26  0.0144    
# 4 FIN_TYP_names Lagoons - Fjords and fjaerds                    0     1.45     0.329      4.41  0.000166  
# 5 FIN_TYP_names Small deltas - Fjords and fjaerds               0     1.26     0.320      3.94  0.00139   
# 6 FIN_TYP_names Tidal systems - Fjords and fjaerds              0     1.58     0.312      5.06  0.00000405



# # # # # # # # # # # # # # # # # # # #
# 2.1 BOXPLOTS - C stocks by SPECIES ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data %>% 
  group_by(Species) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_species, by = "Species") %>% 
  mutate(genus = str_split(string = Species, pattern = " ", simplify = TRUE)[,1],
         abb_genus = str_c(str_sub(genus, 1, 1), "."),
         epitetsp = str_split(string = Species, pattern = " ", simplify = TRUE)[,2],
         abb_species = str_c(abb_genus, epitetsp, sep = " "))

plot_sp <- cstocks_data %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  left_join(letters_multcompar_species, by = "Species") %>% 
  mutate(genus = str_split(string = Species, pattern = " ", simplify = TRUE)[,1],
         abb_genus = str_c(str_sub(genus, 1, 1), "."),
         epitetsp = str_split(string = Species, pattern = " ", simplify = TRUE)[,2],
         abb_species = str_c(abb_genus, epitetsp, sep = " ")) %>% 
  ggplot(aes(x = reorder(abb_species, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(aes(fill = abb_species, colour = abb_species)) +
  scale_colour_d3("category20") +
  scale_fill_d3("category20") +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 5), size = 2) +
  geom_text(data = cstocks_data_box, aes(x = abb_species, y = max_stock + 8, label = letters), size = 4) +
  xlab("Seagrass species") +
  ylab(bquote(~C[20]~ 'stock (Mg C' ~ha^-1* ')')) +
  # coord_flip() +
  # facet_wrap(~depth) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(face = "italic", 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))
# ggsave(file = "Figs/Final_Figures_March2021/All_cores/cstocks_by_Species.pdf")



# # # # # # # # # # # # # # # # # # # #
# 2.2 BOXPlotting FIN_TYP_names ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(n = n(),
            mean_stock = mean(cstocks, na.rm = T),
            median_stock = median(cstocks, na.rm = T),
            std.error = std.error(cstocks, na.rm = T), 
            max_stock = max(cstocks, na.rm = T)) %>% 
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names")


plot_geomorph <- cstocks_data %>%
  # cstocks_tidy %>%
  filter(!is.na(FIN_TYP_names)) %>% 
  filter(!is.na(cstocks)) %>% 
  # filter(Nearest_distance < mean(Nearest_distance)) %>%
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names") %>% 
  ggplot(aes(x = reorder(FIN_TYP_names, cstocks, FUN = median), y = cstocks)) +
  geom_boxplot(fill = "#B2B2B2", colour = "#B2B2B2") +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 5), size = 2) +
  # stat_summary(fun.data = a_fun, geom = "text", position = position_nudge(y = 2)) +
  geom_text(data = cstocks_data_box, aes(x = FIN_TYP_names, y = max_stock + 8, label = letters), size = 4) +
  xlab("Coastal morphology") +
  ylab(bquote(~C[20]~ 'stock (Mg C' ~ha^-1* ')')) +
  # coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))
# ggsave(file = "Figs/Final_Figures_March2021/All_cores/cstocks_vs_FINtype.pdf")

library(cowplot)
plot_grid(plot_sp, 
          plot_geomorph + ylab(NULL),
          ncol = 2, nrow = 1, labels = "AUTO", align = "hv", rel_widths = c(4,3))
# ggsave(filename = "Figs/Final_Figures_March2021/All_cores/panelModelSpeciesSIZED_NOPOSI.pdf", width = 185, height = 130, units = "mm")


# R squared, to report in the paper
MuMIn::r.squaredGLMM(mfinal)
#       R2m       R2c
# 0.3234583 0.9298394 Sí posi oceanica
# 0.2714486 0.8972954 No posi oceanica

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 3. Modelling the response 20cm-carbon stocks. Focusing on the effect of MEADOW TYPE (mono vs multispecific) ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

rm(list = ls())

library(nlme)
library(car)
library(modelr)
library(multcomp)
library(ggsci)

source("01_DataImport&Corrections_CarbonReview.R")
source("mcheck_function.R")
rm(list = c("cstocks_tidy", "cstocks")) # Because we want to use, only, the 20 cm stocks data, not carbon densities.

# Meadow type cannot be assessed with the previous model, because there, we're assessing the effect of Species, and
# by definition, the effect of Species must be tested on MONOSPECIFIC meadows.
# We separate the variable Species into 5 species columns, to be able to know the different species present in multispecific meadows.
cstocksSpeciesSep <- cstocks_tidy20Stocks %>% 
  mutate(Lat_Zone_Posi = ifelse(Posi == "Posi", "Posi", Lat_Zone),
         Lat_Zone = factor(Lat_Zone),
         CoreID_unique = factor(CoreID_unique),
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
listSpeciesMulti <- c(listSpeciesMulti, "Amphibolis spp", "Posidonia spp")

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

m1.lme <- lme(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
              random = ~1|CoreID_unique,
              weights = varIdent(form = ~1|FIN_TYP_names), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)
anova(m0.gls.w, m1.lme) # The random effect CoreID_unique is needed

m2.lme <- lme(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
              random = ~1|Source/CoreID_unique,
              weights = varIdent(form = ~1|FIN_TYP_names), 
              control = lmeControl(maxIter = 100, msMaxIter = 100),
              data = cstocks_data)

AIC(m0.lme, m1.lme, m2.lme)
anova(m1.lme, m2.lme) # Clearly, the model with CoreID_unique nested in Source is better.

# Model selection
m2.lme <- update(m2.lme, method = "ML")
m3.lme <- update(m2.lme, ~. -Meadow_type)
anova(m2.lme, m3.lme) # We can't drop Meadow_type

m4.lme <- update(m2.lme, ~. -depth)
anova(m2.lme, m4.lme) #  We can't drop depth

m5.lme <- update(m2.lme, ~. -FIN_TYP_names)
anova(m2.lme, m5.lme) # FIN_TYP can't be dropped.

mfinal <- lme(log(cstocks) ~ Meadow_type + depth + FIN_TYP_names, 
              random = ~1|Source/CoreID_unique,
              weights = varIdent(form = ~1|FIN_TYP_names), 
              data = cstocks_data)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# 
# Response: log(cstocks)
#                 Chisq Df Pr(>Chisq)    
# Meadow_type    5.3115  1  0.0211851 *  
# depth         14.8509  1  0.0001164 ***
# FIN_TYP_names 31.6374  4  2.269e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(mfinal)


# BUT really, we're not interested in depth and Fin_typ_names. Here we want to focus on meadow-type. Thus, we move them to the random part.
# To take them into account, but not checking its significance.
detach("package:nlme", unload=TRUE)
library(lme4)
mfinal2 <- lmer(log(cstocks) ~ Meadow_type + (1|Source/CoreID_unique) + (1|depth) + (1|FIN_TYP_names), data = cstocks_data)
Anova(mfinal2)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: log(cstocks)
# Chisq Df Pr(>Chisq)  
# Meadow_type 4.4887  1    0.03412 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Model validation
mcheck(mfinal2) # OK


# Multiple comparisons with library(multcompar)
# For Meadow_type 
multcompar_meadowType <- glht(mfinal2, linfct=mcp(Meadow_type="Tukey"))
summary_multcompar_meadowType <- broom::tidy(summary(multcompar_meadowType))
letters_multcompar_meadowType <- broom::tidy(cld(multcompar_meadowType))
# Significant multiple comparisons
summary_multcompar_meadowType %>%
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
  left_join(letters_multcompar_meadowType, by = 'Meadow_type')

cstocks_data %>% 
  ggplot(aes(x = Meadow_type, y = cstocks)) +
  geom_boxplot(aes(fill = Meadow_type, colour = Meadow_type)) +
  stat_summary(fun=median, geom="point", shape = 15, size = 2) +
  scale_fill_manual(values = c("#31a354", "#b2df8a")) +
  scale_colour_manual(values = c("#31a354", "#b2df8a")) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = Meadow_type, y = max_stocks + 5, label = str_c("(", n, ")"))) +
  geom_text(data = cstocksFilteredMonoMulti2, aes(x = Meadow_type, y = max_stocks + 10, label = letters)) +
  scale_x_discrete(labels = c("Monospecific", "Multispecific")) +
  xlab("Meadow type") +
  ylab(bquote(~C[20]~ 'stock (Mg C' ~ha^-1* ')')) +
  # coord_flip() +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("Figs/Final_Figures_March2021/All_cores/meadow_type.pdf", width = 88, height = 82, units = "mm")

# R squared, to report in the paper
MuMIn::r.squaredGLMM(mfinal2)
#       R2m       R2c
# 0.007449527 0.7865999
