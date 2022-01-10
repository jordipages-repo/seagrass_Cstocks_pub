cstocks_tidy20Stocks %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Meadow_type == "monospecific") %>%
  group_by(CoreID_unique) %>% 
  summarise(n = n()) %>%
  filter(n == 1)


cstocks_tidy20Stocks %>% 
  filter(!is.na(cstocks)) %>% 
  filter(Meadow_type == "monospecific") %>%
  group_by(CoreID_unique) %>% 
  summarise(n = n()) %>%
  filter(n == 1)


cstocks_data %>%
  filter(!is.na(cstocks)) %>% 
  mutate(genus = str_split(string = Species, pattern = " ", simplify = TRUE)[,1],
         abb_genus = str_c(str_sub(genus, 1, 1), "."),
         epitetsp = str_split(string = Species, pattern = " ", simplify = TRUE)[,2],
         abb_species = str_c(abb_genus, epitetsp, sep = " ")) %>% 
  group_by(genus) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n/sum(n)*100) %>%
  print(n = Inf)
  ggplot(x = genus, y = percent) +
  geom_bar(aes(x = genus, y = percent), stat = "identity") 

  
  cstocks %>%
    mutate(Species = factor(Species) %>%  fct_infreq() %>% fct_rev()) %>%
    filter(Meadow_type == "monospecific") %>%
    group_by(Species) %>% 
    summarise(n = n()) %>% 
    mutate(percent = 100*(n/sum(n))) %>% 
    ggplot() +
    geom_bar(aes(x = Species, y = percent), stat = "identity") +
    geom_text(aes(x = Species, y = percent, label = str_c(round(percent), "%")), nudge_y = 1, size = 3) +
    ylab("%") +
    coord_flip() +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.text.align = 0,
          text = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(face = "italic"))
