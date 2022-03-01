posi <- plot_grid(plot_sp + xlab(NULL), 
          plot_geomorph + ylab(NULL) + xlab(NULL),
          ncol = 2, nrow = 1, labels = "AUTO", align = "hv", rel_widths = c(4,3))

no_posi <- plot_grid(plot_sp2, 
                     plot_geomorph2 + ylab(NULL),
                     ncol = 2, nrow = 1, labels = "AUTO", align = "hv", rel_widths = c(4,3))



plot_grid(plot_sp + xlab(NULL),
          plot_geomorph + ylab(NULL) + xlab(NULL),
          plot_sp2,
          plot_geomorph2 + ylab(NULL),
          ncol = 2, nrow = 2, labels = "AUTO", align = "hv", rel_widths = c(4,3))
ggsave(filename = "~/Desktop/Fig2_big.pdf", width = 185, height = 185, units = "mm")

fct_relevel(name,
            "`GeomorphEndorheic or Glaciated`", "`GeomorphFjords and fjaerds`", "GeomorphKarst","`GeomorphTidal systems`",
            "L.ratioNP", "`GeomorphSmall deltas`", "L.product.rate", "Ratio.AGBGbiomass", "L.phosphorus", 
            "Lat_ZoneTropical", "R.helong.rate", "Lat_ZoneTemperate", "R.internode.length", "L.width", "L.ratioCP", "L.fibre",
            "L.turnover", "diff_LvsSurfd13C", "L.length", "R.diameter", "L.nitrogen", "R.plastochrone.interv", "L.plastochrone.interv", "L.ratioCN",
            "L.carbon", "L.breakforce", "GeomorphLagoons", "L.lignin", "ABG.biomass", "L.lifespan", "BG.biomass", "N.leaves.shoot")) %>% 
  