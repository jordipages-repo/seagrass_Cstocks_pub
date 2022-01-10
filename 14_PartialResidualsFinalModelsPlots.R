visreg::visreg(mfinal, xvar = "Species", gg = TRUE, trans = exp, partial = T) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(face = "oblique",
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))


a <- visreg::visreg(mfinal, xvar = "Species", gg = TRUE, trans = exp, partial = T)

b <- ggplot_build(a)

b2 <- b$data[[18]]

pepe <- cstocks_data %>% 
  group_by(Species) %>% 
  summarise(n = n())
pepe$Species <- as.character(pepe$Species)

b2$Species <- c(rep(pepe$Species[1], pepe$n[1]),
                rep(pepe$Species[2], pepe$n[2]),
                rep(pepe$Species[3], pepe$n[3]),
                rep(pepe$Species[4], pepe$n[4]),
                rep(pepe$Species[5], pepe$n[5]),
                rep(pepe$Species[6], pepe$n[6]),
                rep(pepe$Species[7], pepe$n[7]),
                rep(pepe$Species[8], pepe$n[8]),
                rep(pepe$Species[9], pepe$n[9]),
                rep(pepe$Species[10], pepe$n[10]),
                rep(pepe$Species[11], pepe$n[11]),
                rep(pepe$Species[12], pepe$n[12]),
                rep(pepe$Species[13], pepe$n[13]),
                rep(pepe$Species[14], pepe$n[14]),
                rep(pepe$Species[15], pepe$n[15]),
                rep(pepe$Species[16], pepe$n[16]),
                rep(pepe$Species[17], pepe$n[17]))

ggplot(b2, aes(x = reorder(Species, y, FUN = median), y = y)) + 
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(face = "oblique",
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))






visreg::visreg(mfinal, xvar = "FIN_TYP_names", gg = TRUE, trans = exp, partial = T) + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))
