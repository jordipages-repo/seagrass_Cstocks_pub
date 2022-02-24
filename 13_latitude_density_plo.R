# # # # # # #
# Analysing seagrass trait data set provided by Hilary Kennedy
# Jordi F. Pagès
# 12-05-2021
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


# # # # # # # # # # # # # # # # #
# Loading carbon stock data ----
# # # # # # # # # # # # # # # # #

source("01_DataImport&Corrections_CarbonReview.R")
rm(list = c("cstocks_tidy", "cstocks")) # Because we want to use, only, the 20 cm stocks data, not carbon densities.


# # # # # # # # # # # # # # # # #
# Plotting latitude ----
# # # # # # # # # # # # # # # # #

ggplot(data = cstocks_tidy20Stocks, aes(y = Latitude)) +
  geom_histogram(binwidth = 5, fill = "#0063ad") +
  # geom_density(adjust = 1/3, fill = "#0063ad", colour = "#0063ad") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(-55,70)) +
  # geom_freqpoly() +
  # coord_flip() +
  # coord_flip(xlim = c(-36,59)) +
  theme_classic() +
  theme(text = element_text(size=60),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave(filename = "Figs/Final_Figures_March2021/latitude_density3.png")


