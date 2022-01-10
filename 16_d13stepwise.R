# # # # # #
# Analysing d13C with stepwise linear models
# Jordi F. Pagès
# July 2021
# Working from CEAB
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

# 0.1. Loading cstocks and soil d13C data ----
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
  # filter(Species != "Posidonia oceanica") %>% # To check the effect of Species without Posidonia oceanica
  # select(Latitude, Species, Source, depth, cstocks, CoreID_unique) %>% 
  select(CoreID_unique, Latitude, Species, FIN_TYP_names, Nearest_distance, 
         Lat_Zone, Lat_Zone_Posi, Source, depth, cstocks, d13C_surface) %>% 
  filter_all(all_vars(!is.na(.)))  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
cstocks_data$FIN_TYP_names <- droplevels(cstocks_data$FIN_TYP_names)


# 0.2. Loading traits data, to get the leaf d13C ----
traits <- read_csv(file = "Seagrass_traits_nov2020+d13soil.csv")

# According to what we talked with H. Kennedy, we complete surf.d13C for Z.noltii with the value from Z.marina.
traits$surf.d13C[which(traits$Species == "Zostera noltii")] <- -19.7


# 0.3. Now we must join the L.d13C in the traits data set with the corresponding species in cstocks_data ----
cstocks_data <- cstocks_data %>% 
  left_join(traits %>% select(Species, L.d13C), by = "Species") %>%  
  mutate(diff_LvsSurfd13C = L.d13C - d13C_surface) %>% 
  mutate(Species = as_factor(Species)) %>% 
  filter(!is.na(diff_LvsSurfd13C))

cstocks_data$Species <- droplevels(cstocks_data$Species)

# Data exploration
ggplot(cstocks_data, aes(x = reorder(FIN_TYP_names, diff_LvsSurfd13C, median), y = diff_LvsSurfd13C)) +
  geom_boxplot()

ggplot(cstocks_data, aes(x = reorder(Species, diff_LvsSurfd13C, median), y = diff_LvsSurfd13C)) +
  geom_boxplot() + 
  coord_flip()

# ggplot(cstocks_data, aes(x = sqrt(diff_LvsSurfd13C+abs(min(diff_LvsSurfd13C)+0.1)))) +
#   geom_histogram()
# shapiro.test(sqrt(cstocks_data$diff_LvsSurfd13C))
# 
# shapiro.test(cstocks_data$diff_LvsSurfd13C)
# shapiro.test(abs(cstocks_data$diff_LvsSurfd13C))
# shapiro.test(log(cstocks_data$diff_LvsSurfd13C + abs(min(cstocks_data$diff_LvsSurfd13C))+0.1))
             
# # # # # # # # # # # # # # # # # # # # # # #
# STEPWISE LINEAR REGRESSION ----
# # # # # # # # # # # # # # # # # # # # # # #

m1 <- lm(diff_LvsSurfd13C ~ FIN_TYP_names + Species, data = cstocks_data)
step(m1)
Anova(m1)

# Weights?
m1.gls <- gls(diff_LvsSurfd13C ~ FIN_TYP_names + Species, data = cstocks_data)
m1.gls.w <- gls(diff_LvsSurfd13C ~ FIN_TYP_names + Species, 
              weights = varIdent(form = ~1|FIN_TYP_names),
              data = cstocks_data)
anova(m1.gls, m1.gls.w)

m1.gls.w2 <- gls(diff_LvsSurfd13C ~ FIN_TYP_names + Species, 
                 weights = varIdent(form = ~1|Species),
                 data = cstocks_data)
anova(m1.gls, m1.gls.w2)

m1.gls.w3 <- gls(diff_LvsSurfd13C ~ FIN_TYP_names + Species, 
                 weights = varIdent(form = ~1|FIN_TYP_names*Species),
                 data = cstocks_data)
anova(m1.gls.w2, m1.gls.w3) # Both weights needed.


# Random effects?
m1.gls.w3.re1 <- lme(diff_LvsSurfd13C ~ FIN_TYP_names + Species, 
                 weights = varIdent(form = ~1|FIN_TYP_names*Species),
                 random = ~1|Source,
                 control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                 data = cstocks_data)
anova(m1.gls.w3, m1.gls.w3.re1) # Random needed

m1.gls.w3.re2 <- lme(diff_LvsSurfd13C ~ FIN_TYP_names + Species, 
                     weights = varIdent(form = ~1|FIN_TYP_names*Species),
                     random = ~1|Source/CoreID_unique,
                     control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                     data = cstocks_data)
anova(m1.gls.w3.re1, m1.gls.w3.re2)
# CoreID_unique NOT needed

mfinal <- lme(diff_LvsSurfd13C ~ FIN_TYP_names + Species, 
              weights = varIdent(form = ~1|Species*FIN_TYP_names),
              random = ~1|Source,
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# Response: diff_LvsSurfd13C
#                 Chisq Df Pr(>Chisq)    
# FIN_TYP_names  86.821  4  < 2.2e-16 ***
# Species       270.820 12  < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

mcheck(mfinal)
boxplot(resid(mfinal) ~ cstocks_data$FIN_TYP_names)
boxplot(resid(mfinal, type = "pearson") ~ cstocks_data$FIN_TYP_names)
boxplot(resid(mfinal) ~ cstocks_data$Species)
boxplot(resid(mfinal, type = "pearson") ~ cstocks_data$Species)

# Fitted versus response
plot(fitted(mfinal) ~ cstocks_data$diff_LvsSurfd13C, col="darkgrey", 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")


# Multiple comparisons with library(multcompar) FOR SPECIES!
multcompar_species <- glht(mfinal, linfct=mcp(Species="Tukey"))
summary_multcompar_species <- broom::tidy(summary(multcompar_species))
letters_multcompar_species <- broom::tidy(cld(multcompar_species))
# Significant multiple comparisons
summary_multcompar_species %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf) #%>% 
# A tibble: 19 × 7
# term    contrast                                     null.value estimate std.error statistic adj.p.value
# <chr>   <chr>                                             <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 Species Halodule wrightii - Zostera marina                    0    -7.49     1.95      -3.84    5.91e- 3
# 2 Species Posidonia oceanica - Zostera marina                   0    -6.05     1.54      -3.94    4.16e- 3
# 3 Species Posidonia sinuosa - Zostera marina                    0    -6.81     1.23      -5.54    1.04e- 6
# 4 Species Halophila stipulacea - Halodule wrightii              0     9.97     1.98       5.04    2.56e- 5
# 5 Species Halophila stipulacea - Posidonia oceanica             0     8.54     1.49       5.73    2.51e- 7
# 6 Species Halodule uninervis - Halophila stipulacea             0    -3.76     0.708     -5.32    9.14e- 6
# 7 Species Thalassia hemprichii - Halophila stipulacea           0    -4.53     0.606     -7.48    1.39e-12
# 8 Species Enhalus acoroides - Halophila stipulacea              0    -4.67     0.413    -11.3     0       
# 9 Species Cymodocea rotundata - Halophila stipulacea            0    -6.32     1.62      -3.91    4.48e- 3
# 10 Species Posidonia australis - Halophila stipulacea            0    -5.52     0.998     -5.53    1.48e- 6
# 11 Species Posidonia sinuosa - Halophila stipulacea              0    -9.29     1.01      -9.20    0       
# 12 Species Halophila ovalis - Halophila stipulacea               0    -3.67     1.05      -3.50    2.05e- 2
# 13 Species Posidonia sinuosa - Halodule uninervis                0    -5.53     0.958     -5.77    8.72e- 7
# 14 Species Posidonia sinuosa - Thalassodendron ciliatum          0    -6.10     1.80      -3.38    2.99e- 2
# 15 Species Posidonia sinuosa - Thalassia hemprichii              0    -4.76     1.13      -4.21    1.18e- 3
# 16 Species Posidonia sinuosa - Enhalus acoroides                 0    -4.62     1.04      -4.45    5.23e- 4
# 17 Species Posidonia sinuosa - Posidonia australis               0    -3.77     0.624     -6.05    9.64e- 8
# 18 Species Halophila ovalis - Posidonia sinuosa                  0     5.62     0.857      6.55    9.79e-10

# Multiple comparisons with library(multcompar) FOR FIN_TYP_NAMES... we have to include it for a second.
multcompar_FIN_TYP_names <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
summary_multcompar_FIN_TYP_names <- broom::tidy(summary(multcompar_FIN_TYP_names))
letters_multcompar_FIN_TYP_names <- broom::tidy(cld(multcompar_FIN_TYP_names))
# Significant multiple comparisons
summary_multcompar_FIN_TYP_names %>%
  filter(adj.p.value<0.05) %>%
  print(n = Inf)
# A tibble: 4 × 7
# term          contrast                null.value estimate std.error statistic adj.p.value
# <chr>         <chr>                        <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
# 1 FIN_TYP_names Lagoons - Arheic                 0    10.0      1.14       8.76    0       
# 2 FIN_TYP_names Small deltas - Arheic            0     4.18     0.810      5.16    1.56e- 6
# 3 FIN_TYP_names Small deltas - Lagoons           0    -5.82     0.935     -6.23    2.53e- 9
# 4 FIN_TYP_names Tidal systems - Lagoons          0    -7.08     1.07      -6.59    2.19e-10


# # # # # # # # # # # # # # # # # # # #
# 1 BOXPLOTS - C stocks by SPECIES ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data %>% 
  group_by(Species) %>% 
  summarise(n = n(),
            mean_stock = mean(diff_LvsSurfd13C, na.rm = T),
            median_stock = median(diff_LvsSurfd13C, na.rm = T),
            std.error = std.error(diff_LvsSurfd13C, na.rm = T), 
            max_stock = max(diff_LvsSurfd13C, na.rm = T)) %>% 
  left_join(letters_multcompar_species, by = "Species") %>% 
  mutate(genus = str_split(string = Species, pattern = " ", simplify = TRUE)[,1],
         abb_genus = str_c(str_sub(genus, 1, 1), "."),
         epitetsp = str_split(string = Species, pattern = " ", simplify = TRUE)[,2],
         abb_species = str_c(abb_genus, epitetsp, sep = " "))

plot_sp <- cstocks_data %>%
  left_join(letters_multcompar_species, by = "Species") %>% 
  mutate(genus = str_split(string = Species, pattern = " ", simplify = TRUE)[,1],
         abb_genus = str_c(str_sub(genus, 1, 1), "."),
         epitetsp = str_split(string = Species, pattern = " ", simplify = TRUE)[,2],
         abb_species = str_c(abb_genus, epitetsp, sep = " ")) %>% 
  ggplot(aes(x = reorder(abb_species, diff_LvsSurfd13C, FUN = median), y = diff_LvsSurfd13C)) +
  geom_boxplot(aes(fill = abb_species, colour = abb_species)) +
  scale_colour_d3("category20") +
  scale_fill_d3("category20") +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 2), size = 2) +
  geom_text(data = cstocks_data_box, aes(x = abb_species, y = max_stock + 3, label = letters), size = 4) +
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
# ggsave(file = "Figs/Final_Figures_March2021/diff13C_by_Species.pdf")



# # # # # # # # # # # # # # # # # # # #
# 2 BOXPlotting FIN_TYP_names ----
# # # # # # # # # # # # # # # # # # # #

cstocks_data_box <- cstocks_data %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(n = n(),
            mean_stock = mean(diff_LvsSurfd13C, na.rm = T),
            median_stock = median(diff_LvsSurfd13C, na.rm = T),
            std.error = std.error(diff_LvsSurfd13C, na.rm = T), 
            max_stock = max(diff_LvsSurfd13C, na.rm = T)) %>% 
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names")


plot_geomorph <- cstocks_data %>%
  # cstocks_tidy %>%
  filter(!is.na(FIN_TYP_names)) %>% 
  filter(!is.na(diff_LvsSurfd13C)) %>% 
  # filter(Nearest_distance < mean(Nearest_distance)) %>%
  left_join(letters_multcompar_FIN_TYP_names, by = "FIN_TYP_names") %>% 
  ggplot(aes(x = reorder(FIN_TYP_names, diff_LvsSurfd13C, FUN = median), y = diff_LvsSurfd13C)) +
  geom_boxplot(fill = "#B2B2B2", colour = "#B2B2B2") +
  stat_summary(fun = median, geom = "crossbar", width = 0.8, lwd = 0.2) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 2), size = 2) +
  # stat_summary(fun.data = a_fun, geom = "text", position = position_nudge(y = 2)) +
  geom_text(data = cstocks_data_box, aes(x = FIN_TYP_names, y = max_stock + 3, label = letters), size = 4) +
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
# ggsave(file = "Figs/Final_Figures_March2021/d13C_vs_FINtype.pdf")

library(cowplot)
plot_grid(plot_sp, 
          plot_geomorph + ylab(NULL),
          ncol = 2, nrow = 1, labels = "AUTO", align = "hv", rel_widths = c(4,3))
# ggsave(filename = "Figs/Final_Figures_March2021/panelModeld13COK.pdf", width = 185, height = 130, units = "mm")


# # # # 



## # # 

# Checking relationship between C provenance and Cstocks



# Full linear model
m0 <- lm(log(cstocks) ~ Species + FIN_TYP_names + Lat_Zone + diff_LvsSurfd13C, data = cstocks_data)
Anova(m0)

# Do we need SPECIES weights?
m0.gls <- gls(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, data = cstocks_data)
m0.gls.w <- gls(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
                weights = varIdent(form = ~1|Species), 
                data = cstocks_data)
anova(m0.gls, m0.gls.w) # We do need SPECIES weights!

# Do we need FIN_TYP_names weights?
m0.gls.w2 <- gls(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
                 weights = varIdent(form = ~1|FIN_TYP_names), 
                 data = cstocks_data)
anova(m0.gls, m0.gls.w2) # We do need FIN_TYP_names weights!

# Do we need both type of weights?
m0.gls.w3 <- gls(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
                 weights = varIdent(form = ~1|Species*FIN_TYP_names), 
                 data = cstocks_data)
AIC(m0.gls, m0.gls.w, m0.gls.w2, m0.gls.w3)
anova(m0.gls.w, m0.gls.w3)
# No need to add both. Species is enough.

# Therefore, we need weights for SPECIES only

# Do we need random effects?
m0.w.lme <- lme(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
                random = ~1|Source,
                weights = varIdent(form = ~1|Species), 
                control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                data = cstocks_data)
anova(m0.w.lme, m0.gls.w) # Source is needed as random intercept

m0.w.lme1 <- lme(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
                 random = ~1|CoreID_unique,
                 weights = varIdent(form = ~1|Species), 
                 control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                 data = cstocks_data)
anova(m0.w.lme1, m0.gls.w) # CoreID_unique is needed as random intercept

m0.w.lme2 <- lme(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
                 random = ~1|Source/CoreID_unique,
                 weights = varIdent(form = ~1|Species), 
                 control = lmeControl(maxIter = 1000, msMaxIter = 1000),
                 data = cstocks_data)
anova(m0.w.lme, m0.w.lme2) # Simple model is better 

# m0.w.lme3 <- lme(log(cstocks) ~ Species*depth + FIN_TYP_names + Lat_Zone,
#                  random = ~1 + CoreID_unique|Source,
#                  weights = varIdent(form = ~1|Species),
#                  control = lmeControl(maxIter = 1000, msMaxIter = 1000),
#                  data = cstocks_data_both_depths)
# anova(m0.w.lme, m0.w.lme3) # Nested model is better than random intercept and slope.

AIC(m0.w.lme, m0.w.lme1, m0.w.lme2)
#           df      AIC
# m0.w.lme  51 662.2859This is the best model.
# m0.w.lme1 51 611.6578
# m0.w.lme2 52 596.2468   

# Model selection 
m1.w.lme.full <- update(m0.w.lme, method = "ML")
m1.w.lme.sub <- update(m1.w.lme.full, .~. -Species)
anova(m1.w.lme.full, m1.w.lme.sub) # Species CANNOT be dropped.

m1.w.lme.sub2 <- update(m1.w.lme.full, ~. -FIN_TYP_names)
anova(m1.w.lme.full, m1.w.lme.sub2) # FIN_TYP_names CANNOT be dropped... 

m2.w.lme.sub <- update(m1.w.lme.full, ~. -diff_LvsSurfd13C)
anova(m1.w.lme.full, m2.w.lme.sub) # diff_LvsSurfd13C CAN be dropped... 

mfinal <- lme(log(cstocks) ~ Species + FIN_TYP_names + diff_LvsSurfd13C, 
              random = ~1|Source,
              weights = varIdent(form = ~1|Species),
              method = "REML",
              control = lmeControl(maxIter = 1000, msMaxIter = 1000),
              data = cstocks_data)
Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# 
# Response: log(cstocks)
#                   Chisq Df Pr(>Chisq)    
# Species          44.317 12  1.348e-05 ***
# FIN_TYP_names    21.921  4  0.0002078 ***
# diff_LvsSurfd13C 11.487  1  0.0007010 ***




mcheck(mfinal)
boxplot(resid(mfinal) ~ cstocks_data$FIN_TYP_names)
boxplot(resid(mfinal, type = "pearson") ~ cstocks_data$FIN_TYP_names)
boxplot(resid(mfinal) ~ cstocks_data$Species)
boxplot(resid(mfinal, type = "pearson") ~ cstocks_data$Species)

# Fitted versus response
plot(fitted(mfinal) ~ log(cstocks_data$cstocks), col="darkgrey", 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")


a <- visreg::visreg(mfinal, xvar = "diff_LvsSurfd13C", trans = exp, gg = T, partial = T)
b <- ggplot_build(a)
b2 <- b$data[[2]]
b2$Species <- cstocks_data$Species
visreg::visreg(mfinal, xvar = "diff_LvsSurfd13C", gg = T, trans = exp, partial = F, rug = F, line=list(col="black", size = 0.5)) +
  geom_point(data = b2, aes(x = x, y = y, colour = Species)) +
  # stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20")#,
                  # alpha = 0.5,
                  # guide = guide_legend(override.aes = list(size = 3,
                                                           # alpha = 1))) +
  xlab(bquote(Delta^13~'C')) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))
