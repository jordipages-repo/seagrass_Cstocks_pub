# # # # # #
# Analysing the data set provided by Hilary Kennedy
# Jordi F. Pagès
# 12-03-2019
# University of Barcelona
# # # # # # 


# # # # # # # # #
# LIBRARIES ----
# # # # # # # # #

library(tidyverse)
library(tidylog)
library(stringr)
library(forcats)
# library(RColorBrewer)
# library(ggmap)

# Function to calculate std. error.
std.error <- function(x, na.rm = TRUE){
  sd(x, na.rm = TRUE)/sqrt(length(x))
}

# Function to calculate n in ggplots
n_fun <- function(x){
  return(data.frame(y = max(x), label = paste0("(", length(x), ")")))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loading the main data set and checking for errors ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cstocks <- read_csv(file = "StockSheet_mono+multispecific_July2019.csv")
glimpse(cstocks)

cstocks <- cstocks %>% 
  mutate(
    CoreID_unique = 1:length(.$Species),
    Species = recode(Species, 
                     Huninervis = "Halodule uninervis",
                     Paustralis = "Posidonia australis",
                     `Thalasia hemprichii` = "Thalassia hemprichii",
                     `Thalassia hemimprichii` = "Thalassia hemprichii",
                     `Z. muelleri` = "Zostera muelleri",
                     `Z. nigricaulis` = "Zostera nigricaulis",
                     `Amphibolis grifficiae` = "Amphibolis griffithii",
                     `Thalassodendrum ciliatum` = "Thalassodendron ciliatum",
                     `Zostera nigracaulis` = "Zostera nigricaulis"),
    Genus = ifelse(Meadow_type == "multispecific", 
                   "multispecific", 
                   str_split(string = Species, pattern = " ", simplify = TRUE)[,1]),
    Country = recode(Country,
                     Keyna = "Kenya",
                     Bahama = "Bahamas",
                     `Abu Dhabi` = "United Arab Emirates"),
    Source = recode(Source,
                    `Serano et al., 2019` = "Serrano et al., 2019",
                    `Serrano unpublished` = "Serrano unpublished data",
                    `Serrano O unpublished data` = "Serrano unpublished data",
                    `Kennedy unpblished data` = "Kennedy unpublished data",
                    `Kennedy unpublished` = "Kennedy unpublished data",
                    `Miyajima, unpublished data` = "Miyajima unpublished data",
                    `Miyajima, T., In prep` = "Miyajima unpublished data",
                    `Miyajima, T In prep` = "Miyajima unpublished data",
                    `Miyajima in prep` = "Miyajima unpublished data",
                    `Miyajima, T., unpublished data` = "Miyajima unpublished data",
                    `Peter Macready unpublished data` = "Macreadie unpublished data",
                    `Arias-Ortiz et al., unpublished` = "Arias-Ortiz unpublished data",
                    `PARI, unpublished data` = "PARI unpublished data"
                    ),
    Published = ifelse(str_detect(Source, "unpublished"), "No", "Yes"),
    Latitude = str_replace(Latitude,  "[°]", ""),
    Longitude = str_replace(Longitude,  "[°]", "")
    ) 

cstocks$Latitude <- as.numeric(cstocks$Latitude)
cstocks$Longitude <- as.numeric(cstocks$Longitude)

# cstocks <- cstocks %>% 
#   mutate(Latitude = round(Latitude, 2),
#          Longitude = round(Longitude, 2))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# TIDYING THE DATA SET TO HAVE ALL DATA TOGETHER ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# We first add a column to identify posidonia oceanica observations
# We then gather the data set to be able to plot all depths in facets.
cstocks_tidy <- cstocks %>% 
  mutate(Posi = ifelse(Meadow_type == "multispecific", "zMixed", ifelse(Species == "Posidonia oceanica", "Posi", "qMonospecific")),
         affinity = ifelse(abs(Latitude) <= 23, "tropical", ifelse(abs(Latitude) > 23 & Posi == "Posi", "Posi", "temperate")),
         Cdensity_0_20 = Stock_0_20cm/20,
         Cdensity_20_50 = Stock_20_50cm/30) %>% 
  select(-Stock_20_50cm_estimate, -Stock_0_50cm) %>% 
  gather(key = depth, value = cstocks, Cdensity_0_20:Cdensity_20_50) %>% 
  select(CoreID_unique, Latitude, Longitude, Species, Meadow_type, Country, Source, Genus, 
         Published, Posi, affinity, depth, cstocks, d13C_surface) %>% 
  mutate(depth = recode(depth, 
                        `Cdensity_0_20` = "0-20 cm",
                        `Cdensity_20_50` = "20-50 cm")) %>% 
  mutate(Species = factor(Species),
         Meadow_type = factor(Meadow_type),
         Country = factor(Country),
         Source = factor(Source),
         Genus = factor(Genus),
         Published = factor(Published),
         depth = factor(depth))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# TIDYING THE DATA SET TO HAVE ALL DATA TOGETHER (FOR 20 cm STOCKS!!!) ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# We first add a column to identify posidonia oceanica observations
# We then gather the data set to be able to plot all depths in facets.
cstocks_tidy20Stocks <- cstocks %>% 
  mutate(Posi = ifelse(Meadow_type == "multispecific", "zMixed", ifelse(Species == "Posidonia oceanica", "Posi", "qMonospecific")),
         affinity = ifelse(abs(Latitude) <= 23, "tropical", ifelse(abs(Latitude) > 23 & Posi == "Posi", "Posi", "temperate")),
         Cdensity_0_20 = Stock_0_20cm/20,
         Cdensity_20_50 = Stock_20_50cm/30,
         OK_20Stock_0_20 = Stock_0_20cm,
         OK_20Stock_20_50 = (Stock_20_50cm/30)*20) %>% 
  select(-starts_with("Stock")) %>% 
  gather(key = depth, value = cstocks, starts_with("OK_20Stock")) %>% 
  select(CoreID_unique, Latitude, Longitude, Species, Meadow_type, Country, Source, Genus, 
         Published, Posi, affinity, depth, cstocks, d13C_surface) %>% 
  mutate(depth = recode(depth, 
                        `OK_20Stock_0_20` = "0-20 cm",
                        `OK_20Stock_20_50` = "20-50 cm")) %>% 
  mutate(Species = factor(Species),
         Meadow_type = factor(Meadow_type),
         Country = factor(Country),
         Source = factor(Source),
         Genus = factor(Genus),
         Published = factor(Published),
         depth = factor(depth))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Loading the geomorphic typology data set,  checking for errors and join to cstocks ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

geomorph <- read_csv(file = "Seagrass_sites_wTypology.csv")
glimpse(geomorph)
# 
# # We want to join cstocks data set with geomorph using lat + long as a key.
# # But we have 576 observation in the cstocks data set and 573 in the geomorph data set. 
# # We'll have to check which lat + long combinations are in cstocks, but not in geomorph.
# length(unique(cstocks_tidy$Latitude))
# length(unique(geomorph$Latitude))
# # We don't have a geomorph row for each latitude... there are 3 geomorph rows missing, compared to cstocks.
# 
# # Checking potential join mismatches
# # anti_join() return all rows from x without a match in y.
# anti_cstocks <- cstocks_tidy %>% 
#   anti_join(geomorph, by = c("Latitude", "Longitude")) %>% 
#   select(Latitude, Longitude)
# # There are 12 combinations of lat+long in cstocks that do not have a geomorph. typology
# 
# anti_geomorph <- geomorph %>% 
#   anti_join(cstocks_tidy, by = c("Latitude", "Longitude")) %>% 
#   select(Latitude, Longitude, FIN_TYP)
# # There are 3 combinations of lat+long in geomorph that do not have a cstock
# 
# left_join(anti_cstocks,anti_geomorph, by = c("Latitude", "Longitude"))
# 
# # This tells me that some rows that to the naked eye seem to be the same, they're not. 
# # And that's because Rstudio doesn't show us all the decimal places when you use view()
# # Checking the raw excel file, I've seen that lat and lon have much more decimals than that shown in Rstudio. 
# 
# # Checking again potential mismatches with anti_joins
# # anti_join() return all rows from x without a match in y.
# anti_join(cstocks_tidy, geomorph, by = c("Latitude", "Longitude"))
# # According to this antijoin, there are 5 rows in cstocks that do not have an equivalent in geomorph.
# 
# 
# # Checking keys for cstocks_tidy
# cstocks_tidy %>% 
#   count(Latitude, Longitude, CoreID_unique) %>% 
#   filter(n > 1)
# # n=2 because we have 2 depths for each CoreID_unique.
# # If we include depth, then we find the unique key.
# cstocks_tidy %>% 
#   count(Latitude, Longitude, CoreID_unique, depth) %>% 
#   filter(n > 1)
# # The key for cstocks_tidy is multivariate and includes all the above variables. 
# # A single observation is identified by Latitude, Longitude, CoreID_unique, depth.
# # THIS IS GOOD
# 
# # THE PROBLEM MUST BE WITH GEOMORPH DATA SET.
# # Checking keys for geomorph
# geomorph %>% 
#   count(Latitude, Longitude) %>% 
#   filter(n > 1)
# # Lots of n>1. This means there are duplicates!!!! AGHH!!!
# # We need to check for duplicates, and if they're really duplicates, delete all of them except one.
# # Then the left join will work!

unique_geomorph <- geomorph %>% 
  group_by(Latitude, Longitude) %>% 
  summarise(n = n(),
            Meadow_type = first(Meadow_type),
            FIN_TYP = first(FIN_TYP),
            Lat_Zone = first(Lat_Zone),
            Nearest_distance = first(`nearestDist (m)`)) %>% 
  select(-Meadow_type, -n)

# Now we've finally got the good geomorph data set.
length(unique(cstocks_tidy$Latitude))
length(unique_geomorph$Latitude)
# We now have 1 more geomorph data point compared to cstocks.


# Checking potential join mismatches, now that we have the correct geomorph data set
# anti_join() return all rows from x without a match in y.
a <- cstocks_tidy %>% 
  anti_join(unique_geomorph, by = c("Latitude", "Longitude")) #%>% 
  # select(Latitude, Longitude)
print(a)
# There are 12 combinations of lat+long in cstocks that do not have a geomorph. typology

# CoreID_unique  Latitude  Longitude Species    Meadow_type  Country  Source                              Genus Published Posi   affinity depth cstocks
# <int>    <dbl>     <dbl> <fct>      <fct>        <fct>    <fct>                               <fct> <fct>     <chr>  <chr>    <fct>   <dbl>
#             1     49.1    -126.  Zostera m… monospecific USA      Postlethwaite et al., 2018 Low blu… Zost… Yes       qMono… tempera… 0-20…   0.41 
#            46    -20.3     -40.3 Halodule … monospecific Brazil   Howard J.L. Creed,J.C.  Aguiar, M.… Halo… Yes       qMono… tropical 0-20…   0.895
#           183     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 0-20…   0.41 
#           184     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 0-20…   0.615
#           355    -32.2     116.  Posidonia… monospecific Austral… Serrano, O; Ricart, AM; Lavery, PS… Posi… Yes       qMono… tempera… 0-20…   0.665
#           576     NA        NA   Zostera n… monospecific Austral… Serrano unpublished data            Zost… No        qMono… NA       0-20…   1.04 
#             1     49.1    -126.  Zostera m… monospecific USA      Postlethwaite et al., 2018 Low blu… Zost… Yes       qMono… tempera… 20-5…  NA    
#            46    -20.3     -40.3 Halodule … monospecific Brazil   Howard J.L. Creed,J.C.  Aguiar, M.… Halo… Yes       qMono… tropical 20-5…   1.30 
#           183     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 20-5…   0.303
#           184     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 20-5…   0.427
#           355    -32.2     116.  Posidonia… monospecific Austral… Serrano, O; Ricart, AM; Lavery, PS… Posi… Yes       qMono… tempera… 20-5…   0.153
#           576     NA        NA   Zostera n… monospecific Austral… Serrano unpublished data            Zost… No        qMono… NA       20-5…   0.937

# Checking what happened with the 12 lat*long combinations that do not have a geomorph.
# Halodule emarginata doesn't have typology data because of an error in Hilary's data set (which still had the º unit in Latitude)
# No Latitude for Z. nigricaulis
# Rounding issue for Z. marina, Enhalus acoroides and Posidonia sinuosa.

# The easiest way to solve this rounding issue is to touch the data itself
cstocks_tidy <- cstocks_tidy %>% 
  mutate(Latitude = if_else(Latitude == 22.61, 22.60583333, Latitude)) %>% # To solve the rounding error for Enhalus 
  mutate(Latitude = if_else(Latitude == 49.13070, 49.13064, Latitude),
         Longitude = if_else(Longitude == -125.55836, -125.55835, Longitude)) %>% # To solve rounding issues for Zostera marina      
  mutate(Latitude = if_else(Latitude == -32.15990, -32.1599015, Latitude),
         Longitude = if_else(Longitude == 115.6717469, 115.671747, Longitude))  # To solve rounding issue for P. sinuosa     

# Checking mismatches again...
cstocks_tidy %>% 
  anti_join(unique_geomorph, by = c("Latitude", "Longitude"))
# Perfect, only Halodule emarginata (for which we really don't have typology data and Z. nigricans, for which we lack coordinates)
# ALL GOOD NOW!

unique_geomorph %>% 
  anti_join(cstocks_tidy, by = c("Latitude", "Longitude")) %>% 
  select(Latitude, Longitude, FIN_TYP)
# There are no combinations of lat+long in geomorph that do not have a cstock!
# ALL PERFECTLY MATCHED NOW!

# Now we've checked the mismatches we can proceed with the left_join
# Left join keeps all rows in x! Since we're interested in cstocks, cstocks must be x.
cstocks_tidy <- left_join(x = cstocks_tidy, y = unique_geomorph, by = c("Latitude", "Longitude"))

# Let's check if we have 4 combinations of lat+long that have FIN_TYP as NA, since not present in geomorph.
cstocks_tidy %>% 
  filter(is.na(FIN_TYP)) #%>%
  # write_csv(file = "cores_with_missing_geomorph_data.csv")

# We do have these 4 NA. Which makes sense.



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# And now we also add geomorph to cstocks_tidy20Stocks -----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Checking potential join mismatches, now that we have the correct geomorph data set
# anti_join() return all rows from x without a match in y.
a <- cstocks_tidy20Stocks %>% 
  anti_join(unique_geomorph, by = c("Latitude", "Longitude")) #%>% 
# select(Latitude, Longitude)
print(a)
# There are 12 combinations of lat+long in cstocks that do not have a geomorph. typology

# CoreID_unique  Latitude  Longitude Species    Meadow_type  Country  Source                              Genus Published Posi   affinity depth cstocks
# <int>    <dbl>     <dbl> <fct>      <fct>        <fct>    <fct>                               <fct> <fct>     <chr>  <chr>    <fct>   <dbl>
#             1     49.1    -126.  Zostera m… monospecific USA      Postlethwaite et al., 2018 Low blu… Zost… Yes       qMono… tempera… 0-20…   0.41 
#            46    -20.3     -40.3 Halodule … monospecific Brazil   Howard J.L. Creed,J.C.  Aguiar, M.… Halo… Yes       qMono… tropical 0-20…   0.895
#           183     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 0-20…   0.41 
#           184     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 0-20…   0.615
#           355    -32.2     116.  Posidonia… monospecific Austral… Serrano, O; Ricart, AM; Lavery, PS… Posi… Yes       qMono… tempera… 0-20…   0.665
#           576     NA        NA   Zostera n… monospecific Austral… Serrano unpublished data            Zost… No        qMono… NA       0-20…   1.04 
#             1     49.1    -126.  Zostera m… monospecific USA      Postlethwaite et al., 2018 Low blu… Zost… Yes       qMono… tempera… 20-5…  NA    
#            46    -20.3     -40.3 Halodule … monospecific Brazil   Howard J.L. Creed,J.C.  Aguiar, M.… Halo… Yes       qMono… tropical 20-5…   1.30 
#           183     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 20-5…   0.303
#           184     22.6      39.3 Enhalus a… monospecific Saudi A… Serrano, O., Almahasheer, H., Duar… Enha… Yes       qMono… tropical 20-5…   0.427
#           355    -32.2     116.  Posidonia… monospecific Austral… Serrano, O; Ricart, AM; Lavery, PS… Posi… Yes       qMono… tempera… 20-5…   0.153
#           576     NA        NA   Zostera n… monospecific Austral… Serrano unpublished data            Zost… No        qMono… NA       20-5…   0.937

# Checking what happened with the 12 lat*long combinations that do not have a geomorph.
# Halodule emarginata doesn't have typology data because of an error in Hilary's data set (which still had the º unit in Latitude)
# No Latitude for Z. nigricaulis
# Rounding issue for Z. marina, Enhalus acoroides and Posidonia sinuosa.

# The easiest way to solve this rounding issue is to touch the data itself
cstocks_tidy20Stocks <- cstocks_tidy20Stocks %>% 
  mutate(Latitude = if_else(Latitude == 22.61, 22.60583333, Latitude)) %>% # To solve the rounding error for Enhalus 
  mutate(Latitude = if_else(Latitude == 49.13070, 49.13064, Latitude),
         Longitude = if_else(Longitude == -125.55836, -125.55835, Longitude)) %>% # To solve rounding issues for Zostera marina      
  mutate(Latitude = if_else(Latitude == -32.15990, -32.1599015, Latitude),
         Longitude = if_else(Longitude == 115.6717469, 115.671747, Longitude))  # To solve rounding issue for P. sinuosa     

# Checking mismatches again...
cstocks_tidy20Stocks %>% 
  anti_join(unique_geomorph, by = c("Latitude", "Longitude"))
# Perfect, only Halodule emarginata (for which we really don't have typology data and Z. nigricans, for which we lack coordinates)
# ALL GOOD NOW!

unique_geomorph %>% 
  anti_join(cstocks_tidy20Stocks, by = c("Latitude", "Longitude")) %>% 
  select(Latitude, Longitude, FIN_TYP)
# There are no combinations of lat+long in geomorph that do not have a cstock!
# ALL PERFECTLY MATCHED NOW!

# Now we've checked the mismatches we can proceed with the left_join
# Left join keeps all rows in x! Since we're interested in cstocks, cstocks must be x.
cstocks_tidy20Stocks <- left_join(x = cstocks_tidy20Stocks, y = unique_geomorph, by = c("Latitude", "Longitude"))

# Let's check if we have 4 combinations of lat+long that have FIN_TYP as NA, since not present in geomorph.
cstocks_tidy20Stocks %>% 
  filter(is.na(FIN_TYP)) #%>%
# write_csv(file = "cores_with_missing_geomorph_data.csv")

# We do have these 4 NA. Which makes sense.



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Loading the plant traits data set, checking for errors and join to cstocks data ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #







# Clear those data sets that I won't need.
rm(list = c("a",
            # "traits", 
            "geomorph",
            "unique_geomorph"))


