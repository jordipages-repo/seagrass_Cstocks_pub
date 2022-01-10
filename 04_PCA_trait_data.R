# # # # # #
# Analysing the data set provided by Hilary Kennedy
# Jordi F. Pagès
# 19-03-2019
# CEAB
# # # # # # 

source("1_DataImport&Corrections_CarbonReview.R")

library(vegan)

# Variables that will go into the PCA.
env <- cstocks_traits %>% 
  select(Stock_0_20cm, d13C_0_20cm, Mean_aboveground_biomass, Mean_aboveground_production, Mean_belowground_production,
         Mean_belowground_biomass) %>% 
  na.omit()

# PCA
which(is.na(env))

boxplot(scale(env))
env.pca <- rda(env, scale=T)
biplot(env.pca, type="t", choice=c(1,2))
biplot(env.pca, type="t", choice=c(1,3))
screeplot(env.pca)
summary(env.pca)

# To make bubble plots
rescale <- function(x, maxsize=6){
  (x-min(x))/(diff(range(x)))*maxsize + 0.1
}

# We make biplots with bubbles
pdf(file = "Desktop/pca.tot.pdf")
biplot(env.pca, type="t", choice=c(1,2), col = c("transparent", "steelblue"), main = "Sites")
col.llocs <- as.factor(strtrim(rownames(env), width=6))
points(env.pca, display = "sites", pch = 19, col = col.llocs)
legend("bottomright", legend = c("TollBridge", "GlaslynCob", "MorfaHarlech", "Talsarnau", "Curian", "Fairbourne", "Unknown", "DyfiCCW"), 
       pch = 19, col = unique(col.llocs), bty = "n", cex = 0.7)

biplot(env.pca,  type="t", col = c("transparent", "steelblue"), main = "%Clay")
points(env.pca, display = "sites", pch = 19, cex = rescale(marsh[, "Perc_clay"]), col = rgb(0,0.5,0.1,0.3))
llegenda <- c(max(round(marsh[, "Perc_clay"], digits=1)), as.numeric(round(quantile(x=(marsh[, "Perc_clay"]), 0.5), digits=1)), min(round(marsh[, "Perc_clay"], digits=1)))
llegenda <- paste(llegenda, "%")
mida <- c(round(max(rescale(marsh[, "Perc_clay"])), digits=1), as.numeric(round(quantile(x=rescale(marsh[, "Perc_clay"]), 0.5), digits=1)), round(min(rescale(marsh[, "Perc_clay"])), digits=1))
legend("bottomright", legend = llegenda, pch = 19, col = rgb(0,0.5,0.1,0.3), pt.cex=mida, cex = 0.7, bty = "n")

biplot(env.pca,  type="t", col = c("transparent", "steelblue"), main = "%Sand")
points(env.pca, display = "sites", pch = 19, cex = rescale(marsh[, "Perc_sand"]), col = rgb(0,0.5,0.1,0.3))
llegenda <- c(max(round(marsh[, "Perc_sand"], digits=1)), as.numeric(round(quantile(x=(marsh[, "Perc_sand"]), 0.5), digits=1)), min(round(marsh[, "Perc_sand"], digits=1)))
llegenda <- paste(llegenda, "%")
mida <- c(round(max(rescale(marsh[, "Perc_sand"])), digits=1), as.numeric(round(quantile(x=rescale(marsh[, "Perc_sand"]), 0.5), digits=1)), round(min(rescale(marsh[, "Perc_sand"])), digits=1))
legend("bottomright", legend = llegenda, pch = 19, col = rgb(0,0.5,0.1,0.3), pt.cex=mida, cex = 0.7, bty = "n")

biplot(env.pca,  type="t", col = c("transparent", "steelblue"), main = "Bulk density")
points(env.pca, display = "sites", pch = 19, cex = rescale(marsh[, "BD_0_10"]), col = rgb(0,0.5,0.1,0.3))
llegenda <- c(max(round(marsh[, "BD_0_10"], digits=1)), as.numeric(round(quantile(x=(marsh[, "BD_0_10"]), 0.5), digits=1)), min(round(marsh[, "BD_0_10"], digits=1)))
llegenda <- paste(llegenda, "g cm-3")
mida <- c(round(max(rescale(marsh[, "BD_0_10"])), digits=1), as.numeric(round(quantile(x=rescale(marsh[, "BD_0_10"]), 0.5), digits=1)), round(min(rescale(marsh[, "BD_0_10"])), digits=1))
legend("bottomright", legend = llegenda, pch = 19, col = rgb(0,0.5,0.1,0.3), pt.cex=mida, cex = 0.7, bty = "n")

biplot(env.pca,  type="t", col = c("transparent", "steelblue"), main = "Erosion Rate")
points(env.pca, display = "sites", pch = 19, cex = rescale(marsh[, "Erosion_rate"]), col = rgb(0,0.5,0.1,0.3))
llegenda <- c(max(round(marsh[, "Erosion_rate"], digits=1)), as.numeric(round(quantile(x=(marsh[, "Erosion_rate"]), 0.5), digits=1)), min(round(marsh[, "Erosion_rate"], digits=1)))
llegenda <- paste(llegenda, "% min-1")
mida <- c(round(max(rescale(marsh[, "Erosion_rate"])), digits=1), as.numeric(round(quantile(x=rescale(marsh[, "Erosion_rate"]), 0.5), digits=1)), round(min(rescale(marsh[, "Erosion_rate"])), digits=1))
legend("bottomright", legend = llegenda, pch = 19, col = rgb(0,0.5,0.1,0.3), pt.cex=mida, cex = 0.7, bty = "n")

biplot(env.pca,  type="t", col = c("transparent", "steelblue"), main = "Elevation")
points(env.pca, display = "sites", pch = 19, cex = rescale(marsh[, "Elevation"]), col = rgb(0,0.5,0.1,0.3))
llegenda <- c(max(round(marsh[, "Elevation"], digits=1)), as.numeric(round(quantile(x=(marsh[, "Elevation"]), 0.5), digits=1)), min(round(marsh[, "Elevation"], digits=1)))
llegenda <- paste(llegenda, "m")
mida <- c(round(max(rescale(marsh[, "Elevation"])), digits=1), as.numeric(round(quantile(x=rescale(marsh[, "Elevation"]), 0.5), digits=1)), round(min(rescale(marsh[, "Elevation"])), digits=1))
legend("bottomright", legend = llegenda, pch = 19, col = rgb(0,0.5,0.1,0.3), pt.cex=mida, cex = 0.7, bty = "n")

biplot(env.pca,  type="t", col = c("transparent", "steelblue"), main = "Salinity")
points(env.pca, display = "sites", pch = 19, cex = rescale(marsh[, "Salinity"]), col = rgb(0,0.5,0.1,0.3))
llegenda <- c(max(round(marsh[, "Salinity"], digits=1)), as.numeric(round(quantile(x=(marsh[, "Salinity"]), 0.5), digits=1)), min(round(marsh[, "Salinity"], digits=1)))
llegenda <- paste(llegenda, "PSU")
mida <- c(round(max(rescale(marsh[, "Salinity"])), digits=1), as.numeric(round(quantile(x=rescale(marsh[, "Salinity"]), 0.5), digits=1)), round(min(rescale(marsh[, "Salinity"])), digits=1))
legend("bottomright", legend = llegenda, pch = 19, col = rgb(0,0.5,0.1,0.3), pt.cex=mida, cex = 0.7, bty = "n")

dev.off()

