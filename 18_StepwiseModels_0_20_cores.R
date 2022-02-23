# # # # # #
# Analysing 0-20 cm CSTOCKS with stepwise linear models
# Jordi F. Pagès
# February 2022
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
  filter(depth == "0-20 cm") %>% 
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
  filter(Species != "Posidonia oceanica") %>% # To check the effect of Species without Posidonia oceanica
  # select(Latitude, Species, Source, depth, cstocks, CoreID_unique) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Nearest_distance, Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks) %>% 
  filter_all(all_vars(!is.na(.)))  
# filter(Species != "Posidonia oceanica" | cstocks != 1.9)

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1. Modelling the response 20cm-carbon stocks. ONLY CORES WITH 0-20 cm ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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



# Fixed effects model selection 
m1.w.lme.full <- update(m0.w.lme, method = "ML")
m1.w.lme.sub <- update(m1.w.lme.full, .~. -Lat_Zone)
anova(m1.w.lme.full, m1.w.lme.sub) # Lat_Zone CAN be dropped.

m2.w.lme.full <- m1.w.lme.sub
m2.w.lme.sub <- update(m2.w.lme.full, ~. -FIN_TYP_names)
anova(m2.w.lme.full, m2.w.lme.sub) # FIN_TYP_names CANNOT be dropped... 

m3.w.lme.sub <- update(m2.w.lme.full, .~. -Species)
anova(m2.w.lme.full, m3.w.lme.sub) # Species CANNOT be dropped... 

mfinal <- lme(log(cstocks) ~ Species + FIN_TYP_names, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species*FIN_TYP_names),
              method = "REML",
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data)
Anova(mfinal)
# WITH POSIDONIA!
# Analysis of Deviance Table (Type II tests)
# 
# Response: log(cstocks)
#                Chisq Df Pr(>Chisq)    
# Species       116.011 16  < 2.2e-16 ***
# FIN_TYP_names  56.911  6  1.904e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# WITHOUT POSIDONIA!
# Analysis of Deviance Table (Type II tests)
# 
# Response: log(cstocks)
#                Chisq Df Pr(>Chisq)    
# Species       53.443 15  3.254e-06 ***
# FIN_TYP_names 58.996  6  7.197e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(mfinal)

mfinal_lme4 <- lme4::lmer(log(cstocks) ~ Species + FIN_TYP_names + (1|Source), data = cstocks_data)
specr::icc_specs(mfinal_lme4) %>%
  mutate_if(is.numeric, round, 2)
specr::plot_variance(mfinal_lme4)

# Checking if we solved heterogeneity
boxplot(resid(mfinal) ~ Species, cstocks_data)
boxplot(resid(mfinal, type = "pearson") ~ Species, cstocks_data) # FANTASTIC! Check difference between both type of residuals.

boxplot(resid(mfinal) ~ FIN_TYP_names, cstocks_data)
boxplot(resid(mfinal, type = "pearson") ~ FIN_TYP_names, cstocks_data) # FANTASTIC! Check difference between both type of residuals.

# Validation
mcheck(mfinal) # Residuals quite GOOD! A tiny bit of OVERDISPERSION! But this can happen with random effects. 
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
letters_multcompar_species <- broom::tidy(cld(multcompar_species, level = 0.01))
# Significant multiple comparisons
summary_multcompar_species %>%
  filter(adj.p.value<0.01) %>%
  print(n = Inf)# %>%
# write_csv(file = "Figs/Final_Figures_March2021/All_cores/multcomp_Species.csv")

# With POSIDONIA!
# A tibble: 3 × 7
#   term    contrast                                   null.value estimate std.error statistic adj.p.value
# <chr>   <chr>                                           <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 Species Posidonia oceanica - Amphibolis antarctica          0     1.39     0.338      4.10    3.39e- 3
# 2 Species Posidonia oceanica - Halophila stipulacea           0     1.23     0.313      3.92    6.20e- 3
# 3 Species Zostera noltii - Posidonia oceanica                 0    -1.20     0.152     -7.87    5.80e-14


# Without POSIDONIA!
# A tibble: 2 x 7
# term    contrast                                      null.value estimate std.error statistic adj.p.value
# <chr>   <chr>                                              <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 Species Thalassodendron ciliatum - Enhalus acoroides          0   -0.251    0.0633     -3.98     0.00463
# 2 Species Thalassia hemprichii - Halophila stipulacea           0    0.416    0.110       3.79     0.00967


# Multiple comparisons with library(multcompar) FOR FIN_TYP_NAMES
multcompar_FIN_TYP_names <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
summary_multcompar_FIN_TYP_names <- broom::tidy(summary(multcompar_FIN_TYP_names))
letters_multcompar_FIN_TYP_names <- broom::tidy(cld(multcompar_FIN_TYP_names, level = 0.01))
# Significant multiple comparisons
summary_multcompar_FIN_TYP_names %>%
  filter(adj.p.value<0.01) %>%
  print(n = Inf)

# WITH POSIDONIA!
# A tibble: 8 × 7
#   term          contrast                               null.value estimate std.error statistic adj.p.value
# <chr>         <chr>                                       <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 FIN_TYP_names Fjords and fjaerds - Arheic                     0   -1.52      0.336     -4.52 0.0000803  
# 2 FIN_TYP_names Lagoons - Endorheic or Glaciated                0    1.73      0.487      3.56 0.00492    
# 3 FIN_TYP_names Tidal systems - Endorheic or Glaciated          0    1.75      0.485      3.60 0.00432    
# 4 FIN_TYP_names Lagoons - Fjords and fjaerds                    0    1.69      0.326      5.19 0.00000234 
# 5 FIN_TYP_names Small deltas - Fjords and fjaerds               0    1.26      0.320      3.94 0.00130    
# 6 FIN_TYP_names Tidal systems - Fjords and fjaerds              0    1.71      0.311      5.49 0.000000382
# 7 FIN_TYP_names Small deltas - Lagoons                          0   -0.428     0.108     -3.94 0.00121    
# 8 FIN_TYP_names Tidal systems - Small deltas                    0    0.444     0.125      3.56 0.00532    


# Without POSIDONIA!
# term          contrast                               null.value estimate std.error statistic adj.p.value
# <chr>         <chr>                                       <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 FIN_TYP_names Fjords and fjaerds - Arheic                 0   -1.41      0.320     -4.40 0.000156   
# 2 FIN_TYP_names Lagoons - Fjords and fjaerds                0    1.65      0.311      5.31 0.00000112 
# 3 FIN_TYP_names Small deltas - Fjords and fjaerds           0    1.16      0.305      3.80 0.00202    
# 4 FIN_TYP_names Tidal systems - Fjords and fjaerds          0    1.63      0.296      5.50 0.000000561
# 5 FIN_TYP_names Small deltas - Lagoons                      0   -0.492     0.112     -4.41 0.000150   
# 6 FIN_TYP_names Tidal systems - Small deltas                0    0.471     0.124      3.80 0.00205    



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
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 10), size = 2) +
  geom_text(data = cstocks_data_box, aes(x = abb_species, y = max_stock + 20, label = letters), size = 3.5) +
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
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 10), size = 2) +
  # stat_summary(fun.data = a_fun, geom = "text", position = position_nudge(y = 2)) +
  geom_text(data = cstocks_data_box, aes(x = FIN_TYP_names, y = max_stock + 20, label = letters), size = 3.5) +
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

cowplot::plot_grid(plot_sp, 
          plot_geomorph + ylab(NULL),
          ncol = 2, nrow = 1, labels = "AUTO", align = "hv", rel_widths = c(4,3))
# ggsave(filename = "~/Desktop/panelModelSpeciesSIZED_NOPOSI2.pdf", width = 185, height = 130, units = "mm")


# R squared, to report in the paper
MuMIn::r.squaredGLMM(mfinal)
#       R2m       R2c
# 0.3689767 0.7704745 Sí posi oceanica
# 0.3327874 0.7057586 No posi oceanica







