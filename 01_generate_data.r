library(tidyverse)
library(factoextra)
library(ggfortify)
library(broom)
library(ggpubr)
library(here)
library(igraph)

# load data -----------------------------------------------------------
load_web <- function(file_name) {
  web <- read_csv(file_name, col_names = FALSE)
  web <- as.matrix(web)
  web[is.na(web)] <- 0
  web[web != 0] <- 1
  web
}

df <- read_csv("data_references.csv") %>%
  mutate(ID = str_c(ID, ".csv")) %>%
  select(-Reference, -`Type of data`) %>%
  mutate(
    web_type =
      case_when(
        str_detect(ID, "A_HP") ~ "antagonistic",
        str_detect(ID, "A_PH") ~ "antagonistic",
        str_detect(ID, "FW_") ~ "antagonistic",
        str_detect(ID, "M_AF") ~ "mutualistic",
        str_detect(ID, "M_PA") ~ "mutualistic",
        str_detect(ID, "M_PL") ~ "mutualistic",
        str_detect(ID, "M_SD") ~ "mutualistic",
      )
  ) %>%
  na.omit()

environment <- raster::getData("worldclim", var = "bio", res = 2.5)
environment <- environment[[c(1, 4, 7, 12, 15)]]
names(environment) <- c("Temp_mean", "Temp_var", "Temp_diff", "Prec_mean", "Prec_var")

points <- sp::SpatialPoints(data.frame(
  x = df$Longitude,
  y = df$Latitude
), proj4string = environment@crs)

data_environemnt <- raster::extract(environment, points) %>%
  as_tibble() %>%
  cbind.data.frame(sp::coordinates(points)) %>%
  rename(
    Longitude = x,
    Latitude = y
  ) %>%
  distinct()

df <- df %>%
  left_join(data_environemnt, by = c("Latitude", "Longitude")) %>%
  rename(Interaction_type = `Type of interactions`)

rm(points, environment, data_environemnt)

# compute statistics ------------------------------------------------------
source("toolbox.r")

data_info <- df %>% 
  mutate(stats = map(ID,~compute_stats_all(., Nrand = 1))) %>% 
  unnest(stats) %>% 
  rename(interaction_type = web_type,
         web_type = Interaction_type,
         type = data_type,
         var1 = v1,
         var2 = v2,
         var3 = v3)

rm(list=setdiff(ls(), "data_info"))
