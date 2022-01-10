# # # # # #
# Analysing the data set provided by Hilary Kennedy
# Jordi F. Pagès
# 20-03-2019
# University of Barcelona
# # # # # # 


# # # # # # # # #
# LIBRARIES --------------------------------------------------------------------------------------------------------------------------------------
# # # # # # # # #

library(nlme)
library(car)
library(visreg)
library(multcomp)
# library(lme4)


# # # # # # # # # # # # # # # # # # # # # 
# Loading the data clean and functions  #
# # # # # # # # # # # # # # # # # # # # # 

source("1_DataImport&Corrections_CarbonReview.R")
source('mcheck_function.R')


# We have to convert all predictors that are nominal into factors!!!
cstocks <- cstocks %>% 
  mutate(Species = factor(Species),
         Meadow_type = factor(Meadow_type),
         Country = factor(Country),
         Source = factor(Source),
         Genus = factor(Genus),
         Published = factor(Published))

cstocks_traits <- cstocks_traits %>% 
  mutate(Species = factor(Species),
         Meadow_type = factor(Meadow_type),
         Country = factor(Country),
         Source = factor(Source),
         Genus = factor(Genus),
         Published = factor(Published))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Data analysis: Modelling using SPECIES as a predictor ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# RESPONSE VARIABLE: Stock_0_20cm --------------------------------------------------------------------------------------------------------
monosp.cstocks020 <- cstocks %>%
  filter(Meadow_type == 'monospecific') %>% 
  mutate(logStock_0_20cm = log(Stock_0_20cm)) %>%
  select(Species, Published, Source, logStock_0_20cm, Stock_0_20cm, Latitude, Longitude) %>% 
  na.omit(.)  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
monosp.cstocks020$Species <- droplevels(monosp.cstocks020$Species) 

# Row 219 is a clear outlier (checking plots previously)
# Step 1. Plain linear model
m1 <- lm(logStock_0_20cm ~ Species, data = monosp.cstocks020[-219,])
plot(m1)
mcheck(m1) # A bit of heterogeneity is obvious, and qqplot is quite good
boxplot(resid(m1)~monosp.cstocks020$Species[-219]) # Heterogeneity is very obvious at the species level

# Steps 2-4. Plain linear model using gls() + and gls with variance structure specified. 
# Variance structure both refers to heterogeneity and random effects
# Dealing with heterogeneity
m1 <- gls(logStock_0_20cm ~ Species, data = monosp.cstocks020[-219,])
m2 <- gls(logStock_0_20cm ~ Species, weights = varIdent(form = ~1|Species), data = monosp.cstocks020[-219,])
anova(m1, m2)
AIC(m1,m2)
# Weights are needed p-val < 0.0001
mcheck(m2)
boxplot(resid(m2, type = "normalized")~monosp.cstocks020$Species[-219]) # Heterogeneity is solved, at species level.

# Dealing with random effects
m3 <- lme(logStock_0_20cm ~ Species, random = ~1|Source, 
          weights = varIdent(form = ~1|Species), 
          data = monosp.cstocks020[-219,], control = lmeControl(maxIter = 100, msMaxIter = 100))
anova(m2, m3)
# Random effect SOURCE is needed p-val < 0.0001
# Since the variable SOURCE already incorporates information on whether a paper is published or not,
# it's better not to include both SOURCE and PUBLISHED variables. They are redundant.
# We will just include SOURCE.
anova(m3) # Not sure why Anova(m3) does not work
#             numDF denDF  F-value p-value
# (Intercept)     1   419 569.9840  <.0001
# Species        23   419   7.3049  <.0001
summary(m3) # the variance due to the random effect is big.
# Amphibolis antarctica used as base level for contrasts.


# MULTIPLE COMPARISONS! 
# Initially not working because of unused levels in the factor Species (a legacy of the subset)
length(coef(m3)) # 25 coefficients, but before, Species had 44 levels. 

# If you don't drop the levels it appears an error when you run glht(). Let´s see if we have the correct nº of levels...
length(levels(monosp.cstocks020$Species)) # 25, correct!

# Multiple comparisons
multcompar <- glht(m3, linfct=mcp(Species="Tukey"))
summary_multcompar <- broom::tidy(summary(multcompar))
letters_multcompar <- broom::tidy(cld(multcompar))

# Significant multiple comparisons
summary_multcompar %>% 
  filter(p.value<0.05) %>% 
  print(n = Inf)


# Model validation
mcheck(m3) # Heterogeneity is much better... but is there some autocorrelation in the residuals?
boxplot(resid(m3, type = "normalized")~monosp.cstocks020$Species[-219]) # Heterogeneity is solved, at species level.

# I HAVE TRIED TO MODEL CORRELATION STRUCTURE, BUT IT DOESN'T REALLY GET BETTER,... WE STILL HAVE THE
# PROBLEM OF 'UNDERDISPERSION', IE. LOWER RESIDUALS THAN EXPECTED AT THE BEGINNING OF QQPLOT



# Data visualisation of model 0-20 cm results ------------------------------------------------------------------------------------------------------------
left_join(monosp.cstocks020, letters_multcompar, by = c('Species' = 'lhs'))

monosp.cstocks020[-219,] %>%
  left_join(letters_multcompar, by = c('Species' = 'lhs')) %>% 
  mutate(resid = resid(m3),
         fitted = fitted(m3)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, fitted, FUN = median), y = exp(fitted)), fill = "#9FDA3AFF") +
  geom_text(aes(x = Species, y = 100, label = letters)) +
  xlab("Species") +
  ylab("Fitted values for carbon stock 0-20 cm") +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y = element_text(face = "italic"))
# ggsave("Figs/Model_cstocks0_20cmVSspecies.pdf")



# RESPONSE VARIABLE: Stock_20_50cm --------------------------------------------------------------------------------------------------------------------------------------
monosp.cstocks2050 <- cstocks %>%
  filter(Meadow_type == 'monospecific') %>% 
  select(Species, Published, Source, Stock_20_50cm, Latitude, Longitude) %>% 
  na.omit(.)

# Step 1. Plain linear model
m1 <- lm(log(Stock_20_50cm) ~ Species, data = monosp.cstocks2050)
plot(m1)
mcheck(m1) # A bit of heterogeneity is obvious, and qqplot is quite good
boxplot(resid(m1)~monosp.cstocks2050$Species) # Heterogeneity is very obvious.

# Steps 2-4. Plain liner model using gls() + and gls with variance structure specified. 
# Variance structure both refers to heterogeneity and random effects
# Dealing with heterogeneity
m1 <- gls(log(Stock_20_50cm) ~ Species, data = monosp.cstocks2050)
m2 <- gls(log(Stock_20_50cm) ~ Species, weights = varIdent(form = ~1|Species), data = monosp.cstocks2050)
anova(m1, m2)
AIC(m1,m2)
# Weights are needed p-val = 0.0223
mcheck(m2)
boxplot(resid(m2, type = "normalized")~monosp.cstocks2050$Species) 
# Heterogeneity is a bit better, at species level.

# Dealing with random effects
m3 <- lme(log(Stock_20_50cm) ~ Species, random = ~1|Source, 
          weights = varIdent(form = ~1|Species), 
          control = lmeControl(maxIter = 100, msMaxIter = 100, returnObject = TRUE),
          data = monosp.cstocks2050)
anova(m2, m3)
# Random effect SOURCE is needed p-val < 0.0001
# Since the variable SOURCE already incorporates information on whether a paper is published or not,
# it's better not to include both SOURCE and PUBLISHED variables. They are redundant.
# We will just include SOURCE.
anova(m3) # Not sure why Anova(m2) does not work
#               numDF denDF   F-value p-value
# (Intercept)     1   155 1478.1463  <.0001
# Species        19   155   15.0024  <.0001
summary(m3)
# Amphibolis antarctica used as base level for contrasts.

# Model validation
mcheck(m3) # There's still a lot of heterogeneity in the residuals
boxplot(resid(m3, type = "normalized")~monosp.cstocks2050$Species) # Heterogeneity is NOT solved

# THERE DOESN'T SEEM TO BE ANY PROBLEM OF SPATIAL AUTOCORRELATION FOR THIS MODEL, ACCORDING TO THE 
# SEMIVARIOGRAM


# Data visualisation of model 0-50 cm results ------------------------------------------------------------------------------------------------------------
monosp.cstocks2050 %>% 
  mutate(resid = resid(m3),
         fitted = fitted(m3)) %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(Species, fitted, FUN = median), y = exp(fitted)), fill = "#9FDA3AFF") +
  theme_bw() +
  xlab("Species") +
  ylab("Fitted values for carbon stock 0-50 cm") +
  theme(axis.text.y = element_text(face = "italic")) +
  coord_flip()
# ggsave("Figs/Model_cstocks0_50cmVSspecies.pdf")




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Data analysis: Modelling using multiple predictors --------------------------------------------------------------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# RESPONSE VARIABLE: Stock_0_20cm --------------------------------------------------------------------------------------------------------
monosp.cstocks020 <- cstocks %>%
  filter(Meadow_type == 'monospecific') %>% 
  mutate(logStock_0_20cm = log(Stock_0_20cm),
         logStock_20_50cm = log(Stock_20_50cm)) %>%
  select(Species, Published, Source, logStock_0_20cm, Stock_0_20cm, Stock_20_50cm, logStock_20_50cm, Latitude, Longitude, d13C_0_20cm) %>% 
  na.omit(.)  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
monosp.cstocks020$Species <- droplevels(monosp.cstocks020$Species) 

# Step 1. Plain linear model
m1 <- lm(logStock_0_20cm ~ Species + logStock_20_50cm + d13C_0_20cm, data = monosp.cstocks020[-c(88,100),])
step(m1)

m1.step <- lm(logStock_0_20cm ~ logStock_20_50cm + d13C_0_20cm, data = monosp.cstocks020[-c(88,100),])
  
plot(m1.step)
mcheck(m1.step) # Perfect

# Dealing with random effects
m1.gls <- gls(logStock_0_20cm ~ Species + logStock_20_50cm + d13C_0_20cm, data = monosp.cstocks020[-c(88,100),])
m2 <- lme(logStock_0_20cm ~ Species + logStock_20_50cm + d13C_0_20cm, random = ~1|Source, 
          data = monosp.cstocks020[-c(88,100),], control = lmeControl(maxIter = 100, msMaxIter = 100))
anova(m1.gls, m2)
# Random effect SOURCE is not needed

m.final <- lm(logStock_0_20cm ~ logStock_20_50cm + d13C_0_20cm, data = monosp.cstocks020[-c(88,100),])
Anova(m.final)
# Anova Table (Type II tests)
# Response: logStock_0_20cm
#                   Sum Sq  Df F value    Pr(>F)    
# logStock_20_50cm 38.079   1 239.277 < 2.2e-16 ***
# d13C_0_20cm       4.551   1  28.599 4.421e-07 ***
# Residuals        18.779 118                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m.final)
# Call:
#   lm(formula = logStock_0_20cm ~ logStock_20_50cm + d13C_0_20cm, 
#      data = monosp.cstocks020[-c(88, 100), ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.97958 -0.23060 -0.01284  0.25038  1.04539 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.327531   0.211511   6.276 5.95e-09 ***
#   logStock_20_50cm 0.741297   0.047923  15.469  < 2e-16 ***
#   d13C_0_20cm      0.050754   0.009491   5.348 4.42e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3989 on 118 degrees of freedom
# Multiple R-squared:  0.7018,	Adjusted R-squared:  0.6968 
# F-statistic: 138.9 on 2 and 118 DF,  p-value: < 2.2e-16

# Model validation
mcheck(m.final) # Great validation

# Visreg
# pdf(file = "Figs/visreg_Stock0_20_allSpecies.pdf")
visreg(m.final, xvar = "logStock_20_50cm", ylab = "log(Stock 0-20 cm)", xlab = "log(Stock 20-50 cm)")
visreg(m.final, xvar = "d13C_0_20cm", ylab = "log(Stock 0-20 cm)", xlab = "d13C 0-20 cm")
# dev.off()


# RESPONSE VARIABLE: Stock_20_50cm --------------------------------------------------------------------------------------------------------
monosp.cstocks2050 <- cstocks %>%
  filter(Meadow_type == 'monospecific') %>% 
  mutate(logStock_0_20cm = log(Stock_0_20cm),
         logStock_20_50cm = log(Stock_20_50cm)) %>%
  select(Species, Published, Source, logStock_0_20cm, Stock_0_20cm, Stock_20_50cm, logStock_20_50cm, Latitude, Longitude, d13C_20_50cm) %>% 
  na.omit(.)  

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
monosp.cstocks2050$Species <- droplevels(monosp.cstocks2050$Species) 

# Step 1. Plain linear model
m1 <- lm(logStock_20_50cm ~ Species + logStock_0_20cm + d13C_20_50cm, data = monosp.cstocks2050)
step(m1)

m1.step <- lm(logStock_20_50cm ~ logStock_0_20cm + d13C_20_50cm, data = monosp.cstocks2050)

plot(m1.step)
mcheck(m1.step) # OK, MORE OR LESS...

# Dealing with random effects
m1.gls <- gls(logStock_20_50cm ~ Species + logStock_0_20cm + d13C_20_50cm, data = monosp.cstocks2050)
m2 <- lme(logStock_20_50cm ~ Species + logStock_0_20cm + d13C_20_50cm, data = monosp.cstocks2050,
          random = ~1|Source, control = lmeControl(maxIter = 100, msMaxIter = 100))
anova(m1.gls, m2)
# Random effect SOURCE is needed, but on the limit.

m.final <- lm(logStock_20_50cm ~ logStock_0_20cm + d13C_20_50cm, data = monosp.cstocks2050)
Anova(m.final)
# Anova Table (Type II tests)
# Response: logStock_0_20cm
#                   Sum Sq  Df F value    Pr(>F)    
# logStock_20_50cm 38.079   1 239.277 < 2.2e-16 ***
# d13C_0_20cm       4.551   1  28.599 4.421e-07 ***
# Residuals        18.779 118                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m.final)
# Call:
#   lm(formula = logStock_0_20cm ~ logStock_20_50cm + d13C_0_20cm, 
#      data = monosp.cstocks020[-c(88, 100), ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.97958 -0.23060 -0.01284  0.25038  1.04539 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.327531   0.211511   6.276 5.95e-09 ***
#   logStock_20_50cm 0.741297   0.047923  15.469  < 2e-16 ***
#   d13C_0_20cm      0.050754   0.009491   5.348 4.42e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3989 on 118 degrees of freedom
# Multiple R-squared:  0.7018,	Adjusted R-squared:  0.6968 
# F-statistic: 138.9 on 2 and 118 DF,  p-value: < 2.2e-16

# Model validation
mcheck(m.final) # Great validation

# Visreg
# pdf(file = "Figs/visreg_Stock0_20_allSpecies.pdf")
visreg(m.final, xvar = "logStock_0_20cm", ylab = "log(Stock 20-50 cm)", xlab = "log(Stock 0-20 cm)")
visreg(m.final, xvar = "d13C_0_20cm", ylab = "log(Stock 0-20 cm)", xlab = "d13C 0-20 cm")
# dev.off()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# Data analysis: Modelling using multiple predictors and TRAITS --------------------------------------------------------------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# RESPONSE VARIABLE: Stock_0_20cm --------------------------------------------------------------------------------------------------------
monosp.cstocks020 <- cstocks_traits %>%
  filter(Meadow_type == 'monospecific') %>% 
  mutate(logStock_0_20cm = log(Stock_0_20cm),
         logStock_0_50cm = log(Stock_0_50cm)) %>%
  select(Species, Published, Source, logStock_0_20cm, Stock_0_20cm, 
         Stock_0_50cm, logStock_0_50cm, Latitude, Longitude, d13C_0_20cm,
         Mean_aboveground_biomass, Mean_belowground_biomass, Root_shoot_ratio) %>% 
  na.omit(.)  

ggplot(monosp.cstocks020) +
  geom_point(aes(x = Mean_belowground_biomass, y = Stock_0_20cm, colour = Species)) +
  geom_smooth(aes(x = Mean_belowground_biomass, y = Stock_0_20cm))

# Otherwise the model might try to calculate coefficients for levels that no longer exist.
monosp.cstocks020$Species <- droplevels(monosp.cstocks020$Species) 


# IF WE PUT SPECIES OR STOCK0-50CM AS PREDICTORS, NONE OF THE TRAITS GET PICKED UP IN THE FINAL MODEL
# Step 1. Plain linear model
m1 <- lm(logStock_0_20cm ~  d13C_0_20cm + Mean_aboveground_biomass +
           Mean_belowground_biomass, data = monosp.cstocks020)
step(m1)

m1.step <- lm(logStock_0_20cm ~ d13C_0_20cm + Mean_belowground_biomass, data = monosp.cstocks020)

plot(m1.step)
mcheck(m1.step) # Perfect

# Dealing with random effects
m1.gls <- gls(logStock_0_20cm ~ Species + logStock_0_50cm + d13C_0_20cm, data = monosp.cstocks020[-c(88,100),])
m2 <- lme(logStock_0_20cm ~ Species + logStock_0_50cm + d13C_0_20cm, random = ~1|Source, 
          data = monosp.cstocks020[-c(88,100),], control = lmeControl(maxIter = 100, msMaxIter = 100))
anova(m1.gls, m2)
# Random effect SOURCE is not needed

m.final <- lm(logStock_0_20cm ~ logStock_0_50cm + d13C_0_20cm, data = monosp.cstocks020[-c(88,100),])
Anova(m.final)
# Anova Table (Type II tests)
# Response: logStock_0_20cm
#                 Sum Sq  Df F value    Pr(>F)    
# logStock_0_50cm 49.763   1 827.719 < 2.2e-16 ***
# d13C_0_20cm      1.644   1  27.337 7.486e-07 ***
# Residuals        7.094 118                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m.final)
# Call:
#   lm(formula = logStock_0_20cm ~ logStock_0_50cm + d13C_0_20cm, 
#      data = monosp.cstocks020[-c(88, 100), ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.59108 -0.18161  0.00468  0.17974  0.66834 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -0.038920   0.157200  -0.248    0.805    
# logStock_0_50cm  0.920275   0.031987  28.770  < 2e-16 ***
# d13C_0_20cm      0.030873   0.005905   5.229 7.49e-07 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2452 on 118 degrees of freedom
# Multiple R-squared:  0.8874,	Adjusted R-squared:  0.8854 
# F-statistic: 464.8 on 2 and 118 DF,  p-value: < 2.2e-16

# Model validation
mcheck(m.final) # Great validation

# Visreg
# pdf(file = "Figs/visreg_Stock0_20_allSpecies.pdf")
visreg(m.final, xvar = "logStock_0_50cm", ylab = "log(Stock 0-20 cm)", xlab = "log(Stock 0-50 cm)")
visreg(m.final, xvar = "d13C_0_20cm", ylab = "log(Stock 0-20 cm)", xlab = "d13C 0-20 cm")
# dev.off()





