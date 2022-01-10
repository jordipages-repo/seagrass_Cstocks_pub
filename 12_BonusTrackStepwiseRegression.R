# # # # # # #
# Analysing seagrass trait data set provided by Hilary Kennedy
# Jordi F. Pagès
# 05-05-2021
# University of Barcelona 
# (currently working from home/CEAB)
# # # # # # 

# # # # # # # # #
# LIBRARIES ----
# # # # # # # # #

library(tidyverse)
library(tidylog)
library(ggsci)
library(cowplot)

source("HighstatLib.r")
source("mcheck_function.R")

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
# We only include important predictors according to VIP PLSR (those >1)
cstocks_traits <- cstocks_data %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Species %in% SpeciesList) %>%
  left_join(traits_full, by = "Species") %>%
  select(Species,
         cstocks,
         Source,
         CoreID_unique,
         N.leaves.shoot,
         BG.biomass,
         L.lifespan,
         ABG.biomass,
         FIN_TYP_names,
         L.lignin,
         L.breakforce,
         L.carbon)


# # # # # # # # # # # # # # # # # # # # # # #
# STEPWISE LINEAR REGRESSION ----
# # # # # # # # # # # # # # # # # # # # # # #

corvif(
  as.matrix(
    cstocks_traits %>% 
      select(-Species, -cstocks, -Source, -CoreID_unique, -FIN_TYP_names)
  )
)

# It suggests us to drop ABG.biomass first (highest VIF)
# Then it suggests to drop BG.biomass, but since it's one of our most interesting variables, I drop the second highest VIF (L.lifespan)

corvif(
  as.matrix(
    cstocks_traits %>% 
      select(-Species, -cstocks, -Source, -CoreID_unique, -FIN_TYP_names, -ABG.biomass, -L.lifespan)
  )
)

# Now VIF values are <10 as suggested by Montgomery, D.C. & Peck, E.A. (1992) Introduction to Linear Regression Analysis. Wiley, NewYork.

m1 <- lm(log(cstocks) ~ N.leaves.shoot + BG.biomass + FIN_TYP_names + L.lignin + L.breakforce + L.carbon, data = cstocks_traits)
step(m1)
m2 <- lm(log(cstocks) ~ BG.biomass + FIN_TYP_names + L.breakforce, data = cstocks_traits)
car::Anova(m2)
visreg::visreg(m2)

library(nlme)
# let's add weights to FIN_TYP_names
m3 <- gls(log(cstocks) ~ N.leaves.shoot + BG.biomass + FIN_TYP_names + L.lignin + L.breakforce + L.carbon, data = cstocks_traits)
m3w <- gls(log(cstocks) ~ N.leaves.shoot + BG.biomass + FIN_TYP_names + L.lignin + L.breakforce + L.carbon, 
          weights =  varIdent(form = ~1|FIN_TYP_names), 
          data = cstocks_traits)
anova(m3, m3w) # weights are needed

# Let's add random effects
m4 <- lme(log(cstocks) ~ N.leaves.shoot + BG.biomass + FIN_TYP_names + L.lignin + L.breakforce + L.carbon,
          weights =  varIdent(form = ~1|FIN_TYP_names), 
          random = ~1|Source,
          data = cstocks_traits)
anova(m3w, m4) # Random effect Source is needed

m5 <- lme(log(cstocks) ~ N.leaves.shoot + BG.biomass + FIN_TYP_names + L.lignin + L.breakforce + L.carbon, 
          weights =  varIdent(form = ~1|FIN_TYP_names), 
          random = ~1|Source/CoreID_unique,
          data = cstocks_traits)
anova(m4, m5) # Random effect CoreID_unique nested in Source is needed

# Start selection of the fixed part
m5.2 <- update(m5, method = "ML")
m5.2sub <- update(m5.2, .~. -L.carbon)
anova(m5.2, m5.2sub) # L.carbon can be dropped

m6 <- update(m5.2sub, .~. -L.breakforce)
anova(m5.2sub, m6) # L.breakforce can be dropped

m7 <- update(m6, .~. -L.lignin)
anova(m6, m7) # L.lignin can be dropped

m8 <- update(m7, .~. -FIN_TYP_names)
anova(m7, m8) # FIN_TYP_names cannot be dropped. (on the limit)

m9 <- update(m7, .~. -BG.biomass)
anova(m7, m9) # BG.biomass cannot be dropped

m10 <- update(m7, .~. -N.leaves.shoot)
anova(m7, m10) # N.leaves.shoot can be dropped

mfinal <- lme(log(cstocks) ~ BG.biomass + FIN_TYP_names, 
              weights = varIdent(form = ~1|FIN_TYP_names),
              random = ~1|Source/CoreID_unique,
              data = cstocks_traits)

car::Anova(mfinal)
# Analysis of Deviance Table (Type II tests)
# Response: log(cstocks)
#                Chisq Df   Pr(>Chisq)    
# BG.biomass    83.337  1    < 2e-16 ***
# FIN_TYP_names 10.746  6    0.09655 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(mfinal)

# Model validation
mcheck(mfinal) 
plot(fitted(mfinal) ~ log(cstocks_traits$cstocks), col="darkgrey", 
     xlab="Y (response)", ylab="Fitted Values")
abline(a=0, b=1, col="red")


# MULTIPLE COMPARISONS! 
# Multiple comparisons with library(multcompar) FOR SPECIES!
library(multcomp)
multcompar <- glht(mfinal, linfct=mcp(FIN_TYP_names="Tukey"))
multcompar_summary <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(summary(multcompar), level = 0.07))
# Significant multiple comparisons
multcompar_summary %>%
  filter(adj.p.value<0.07) %>%
  print(n = Inf) #%>% 


# # # # # # # # # # # # # # # # # # #
# Plotting significant effects ----
# # # # # # # # # # # # # # # # # # #

# BELOWGROUND BIOMASS 
a <- visreg::visreg(mfinal, xvar = "BG.biomass", trans = exp, gg = T, partial = T)
b <- ggplot_build(a)
b2 <- b$data[[2]]
b2$Species <- cstocks_traits$Species
p.BG.biomass <- visreg::visreg(mfinal, xvar = "BG.biomass", gg = T, trans = exp, partial = F, rug = F, line=list(col="black", size = 0.5)) +
  geom_point(data = b2, aes(x = x, y = y, colour = Species)) +
  # stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab(bquote('Belowground biomass (g DW' ~m^-2* ')')) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))
# ggsave(filename = "Figs/Final_Figures_March2021/PLSR/cstocksvsBG.biomassVISREG2.pdf", width = 88, height = 80, units = "mm")


# GEOMORPHOLOGY - not worth plotting this, since effects are marginal.
cstocks_box <- cstocks_traits %>% 
  left_join(letters_multcompar, by = "FIN_TYP_names") %>% 
  group_by(FIN_TYP_names) %>% 
  summarise(max = max(cstocks),
            median = median(cstocks),
            letters = first(letters))
  
p.geomorphology <- ggplot(data = cstocks_traits, aes(x = reorder(FIN_TYP_names, cstocks, median), y = cstocks)) +
  geom_boxplot(fill = "#e25528", colour = "#e25528") + 
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_nudge(y = 20), size = 3) +
  # scale_colour_d3("category20") +
  # scale_fill_d3("category20") +
  geom_text(data = cstocks_box, aes(x = reorder(FIN_TYP_names, median, identity), y = max + 40, label = letters, )) +
  xlab("") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  # coord_flip() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
# ggsave(filename = "Figs/cstocksvsGeomorphology.pdf")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Plotting other important variables according to PLSR ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Number of leaves per shoot
p.N.leaves.shoot <- ggplot(data = cstocks_traits, aes(x = N.leaves.shoot, y = cstocks)) +
  geom_jitter(aes(colour = Species), width = 0.05) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab("Number of leaves per shoot") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))


# Leaf lifespan
p.L.lifespan <- ggplot(data = cstocks_traits, aes(x = L.lifespan, y = cstocks)) +
  geom_jitter(aes(colour = Species), width = 3) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab("Leaf lifespan (days)") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))
  

# ABG.biomass
p.ABG.biomass <- ggplot(data = cstocks_traits, aes(x = ABG.biomass, y = cstocks)) +
  geom_jitter(aes(colour = Species), width = 3) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab(bquote('Mean aboveground biomass (g DW' ~m^-2* ')')) +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))


# L.lignin %DW
p.L.lignin <- ggplot(data = cstocks_traits, aes(x = L.lignin, y = cstocks)) +
  geom_jitter(aes(colour = Species), width = 0.3) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab("Leaf lignin content (% DW)") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))

  
# L.breakforce %DW
p.L.breakforce <- ggplot(data = cstocks_traits, aes(x = L.breakforce, y = cstocks)) +
  geom_jitter(aes(colour = Species), width = 0.3) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab("Leaf breaking force (N)") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))


# L.carbon %DW
p.L.carbon <- ggplot(data = cstocks_traits, aes(x = L.carbon, y = cstocks)) +
  geom_jitter(aes(colour = Species), width = 0.3) +
  stat_summary(fun = median, geom = "point", shape = 15, size = 1) +
  scale_colour_d3("category20", 
                  alpha = 0.5,
                  guide = guide_legend(override.aes = list(size = 3,
                                                           alpha = 1))) +
  xlab("Leaf carbon content (% DW)") +
  ylab(bquote('Carbon stock (Mg C' ~ha^-1* ')')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 14))


# All plots together
plot.grid_top <- plot_grid(p.BG.biomass  + theme(legend.position = "none"), 
                           p.N.leaves.shoot + ylab(NULL) + coord_cartesian(ylim = c(0,300)) + theme(legend.position = "none"), 
                           p.L.lifespan + ylab(NULL) + coord_cartesian(ylim = c(0,300)) + theme(legend.position = "none"),
                           p.ABG.biomass + coord_cartesian(ylim = c(0,300)) + theme(legend.position = "none"),
                           p.L.lignin + coord_cartesian(ylim = c(0,300)) + ylab(NULL) + theme(legend.position = "none"),
                           ncol = 3, nrow = 2, labels = "AUTO", align = "hv")

legend <- get_legend(
  # create some space to the left of the legend
  p.BG.biomass + theme(legend.box.margin = margin(0, 0, 0, 12),
                       legend.title = element_blank(),
                       text = element_text(size = 12))
)

plot.grid_bottom <- plot_grid(p.L.breakforce + coord_cartesian(ylim = c(0,300)) + theme(legend.position = "none"),
                              p.L.carbon + coord_cartesian(ylim = c(0,300)) + ylab(NULL) + theme(legend.position = "none"),
                              legend,
                              labels = c("F", "G"),
                              rel_widths = c(1, 1, 1),
                              nrow = 1)

plot_grid(plot.grid_top, plot.grid_bottom, ncol = 1, rel_heights = c(2,1), align = "hv")
ggsave(filename = "Figs/Final_Figures_March2021/plotgridNOU.pdf", width = 10, height = 8)
