# Sample maps

library(ggplot2)
library(sf)
library(readr)
library(dplyr)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

(sites <- data.frame(longitude = c(-80.144005, -80.109), 
                     latitude = c(26.479005, 26.83)))

# Load site data
popmap <- read_tsv("F:/PopFinder/Data/empirical data/lesp/popmap.txt")
popmap %>%
  filter(!pop %in% c("coa", "cog", "fun"))

sites <- popmap %>% 
  group_by(pop) %>%
  filter(row_number()==1) %>%
  filter(!is.na(pop)) %>%
  dplyr::select(pop,x,y) %>%
  filter(!pop %in% c("coa", "cog", "fun"))

p <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = sites, aes(x = x, y = y), size = 4, 
             shape = 21, fill = "#008080") +
  coord_sf(xlim = c(-70, 15), ylim = c(40, 70), expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_classic()
ggsave("F:/PopFinder/Presentation/images/lesp_map.png")
