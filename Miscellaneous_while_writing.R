visreg::visreg(mfinal, xvar = "BG.biomass", gg = T, trans = exp, partial = F, rug = F, line=list(col="black", size = 0.5)) +
  geom_point(data = b2, aes(x = x, y = y, colour = Species)) +
  geom_line(data = linia, aes(x = BGB, y = `0ver 20cm depth`), colour = "red") +
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
ggsave(filename = "~/Desktop/asiaticPaperLine.pdf")  

  
linia <- read_csv("~/Desktop/Copy of calculation.csv")


cstocks_tidy20Stocks %>% 
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
  group_by(FIN_TYP_names) %>% 
  summarise(n = n()) %>%
  mutate(percent = (n/sum(n))*100) %>% 
  arrange(desc(percent)) %>% 
  print(n = Inf)














library(formattable)
c_table <- cstocks_data %>%
  # cstocks_tidy %>%
  filter(!is.na(Species)) %>% 
  filter(!is.na(cstocks)) %>% 
  mutate(genus = str_split(string = Species, pattern = " ", simplify = TRUE)[,1],
         abb_genus = str_c(str_sub(genus, 1, 1), "."),
         epitetsp = str_split(string = Species, pattern = " ", simplify = TRUE)[,2],
         abb_species = str_c(abb_genus, epitetsp, sep = " ")) %>%
  group_by(Species) %>% 
  summarise(Min = min(cstocks),
            Max = max(cstocks),
            `Max/Min` = max(cstocks)/min(cstocks),
            Mean = mean(cstocks),
            SD = sd(cstocks),
            n = n(),
            CV = percent((SD / Mean), digits = 2)) %>% 
  filter(n > 10) %>% 
  arrange(desc(CV)) %>% 
  mutate(Min = round(Min, 2),
         Max = round(Max, 2),
         `Max/Min` = round(`Max/Min`, 0),
         SD = round(SD, 2),
         Mean = round(Mean, 2),
         CV = round(CV, 4))


formattable(c_table, 
            align = c("l","c","c","c","c","c","c","c"),
            list(CV = color_tile("transparent", "lightblue"),
                 `Max/Min` = color_tile("transparent", "lightseagreen"),
                 Species = formatter("span", style = ~style(font.style = "italic"))))


#####
# After having run the final stepwise model. To decompose variance of fixed effects:
a <- partR2(mfinal_lme4, data = cstocks_data, R2_type = "marginal", nboot = 100, partvars = c("Species", "FIN_TYP_names"))
forestplot(a)

b <- partR2(mfinal_lme4, data = cstocks_data, R2_type = "conditional", nboot = 100, partvars = c("Species", "FIN_TYP_names"))
forestplot(b)

