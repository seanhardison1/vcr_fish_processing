library(tidyverse)
library(sf)
library(raster)
library(readxl)
library(magrittr)
library(plotly)
library(lubridate)


# Official 2019 and 2021 fish sampling coordinates from VCR-LTER
vcr_coord <- read_excel(here::here("data/fishsampling_coordinates.xls")) %>% 
  dplyr::rename(latitude_official = POINT_Y,
                longitude_official = POINT_X)


# Import and reorganize key so that it joins the 2019 data 
# with the official coordinates and synoptic coordinates
loc <- read_excel(here::here("data/synoptic_locs2019.xlsx")) %>% 
  {. ->> key} %>% 
  mutate(fish_site = ifelse(str_detect(fish_site, "-C"), str_replace_all(fish_site,
                                                                          "CB", "CI"),fish_site),
         fish_site = ifelse(fish_site == "FISH-CI6", "FISH-CI6X", fish_site)) %>% 
  right_join(.,vcr_coord, by = c("fish_site" = "Name")) %>% 
  filter(!is.na(fish_site)) %>% 
  mutate(site = str_replace_all(site, " ", "-")) %>% 
  dplyr::select(synoptic_reference = site,
                synoptic_latitude = latitude, 
                synoptic_longitude = longitude,
                fish_latitude = latitude_official,
                fish_longitude = longitude_official,
                fish_site,
                fish_site_syn) %>% 
  mutate(fish_site_syn = ifelse(fish_site_syn == "sb6", 
                                "sb6x",
                                ifelse(fish_site_syn == "sc6", 
                                       "sc6x",fish_site_syn)))

# 2019 data----
# Read in raw fish sampling data from 2019, clean it up, and join it 
# to processed location data from above
fish2019 <- read_excel(here::here("data/2019.fish.seining.10-22-2019.xlsx")) %>% 
  {. ->> og_fish} %>% 
  dplyr::rename(site = SITE) %>% 
  mutate(TIME = format(.$TIME, "%H:%M:%S")) %>% 
  mutate(site = str_to_lower(site),
         site = str_remove_all(site, "-| ")) %>% 
  mutate(site = ifelse(str_detect(site, "ci"),
                       str_replace(site, "ci", "cb"),
                       site)) %>% 
  left_join(.,loc, by = c("site" = "fish_site_syn")) %>% 
  mutate(site_comment = NA,
         observers = NA)

# visualize fish locations and their synoptic reference sites
ggplot(loc) +
  geom_point(aes(y = fish_latitude, x = fish_longitude)) +
  geom_point(aes(y = synoptic_latitude, x = synoptic_longitude),
             color = "red", alpha = 0.25)

# Expand counts such that each row is a different fish
fish2019_updated <- NULL
for (i in 1:nrow(fish2019)){
  df <- fish2019 %>% dplyr::slice(i)
  if (!is.na(df$COUNT_Total)){
    ndf <- df[rep(row.names(df), df$COUNT_Total), ]
    ndf$COUNT_Total <- 1
  } else {
    ndf <- df
  }
  assign('fish2019_updated', rbind(fish2019_updated, ndf))
}

# drop columns and make all column names lowercase
fish2019_final <-
  fish2019_updated %>% 
  dplyr::select(-site, -COUNT_Total, -CODE) %>% 
  rename_with(str_to_lower, DATE:fish_site) %>% 
  dplyr::rename(depth = depth_cm) %>% 
  mutate(do = mean(c_across(c(`do_top_mg/l`,`do_bottom_mg/l`)), na.rm = T),
         temperature = mean(c_across(c(temp_top_c, temp_bottom_c)), na.rm = T),
         salinity = mean(c_across(c(sal_top_permil, sal_bottom_permil)), na.rm = T),
         conductivity = mean(c_across(c(cond_top_ms, cond_bottom_ms)), na.rm = T)) %>% 
  dplyr::select(-`do_top_mg/l`:-cond_bottom_ms) %>% 
  dplyr::rename(length = length_cm,
                length_category = cat_cm) %>% 
  mutate(common = str_to_lower(common),
         common = str_replace_all(common, "atlantic", "Atlantic"),
         common = ifelse(common == "halfbeak", "American halfbeak", common))

# 2021 data----
# Process 2021 data so that it has the same columns as 2019
fish_2021 <- read_csv(here::here("data/FishQueryTable_10.27.2021.csv")) %>%
  dplyr::rename(date = sampleDate,
                time = sampleTime,
                fish_longitude = coordW,
                fish_latitude = coordN,
                fish_site = site,
                common = speciesName) %>% 
  mutate(time = format(.$time, "%H:%M:%S"),
         date = as.Date(date, "%m/%d/%Y"),
         fish_longitude = fish_longitude * -1,
         fish_longitude = ifelse(fish_site == "HI_1", -75.72969, fish_longitude),
         fish_site = str_to_upper(paste0("FISH-",str_remove_all(fish_site, "_"))),
         fish_site = ifelse(fish_site == "FISH-SB1" & year(date) == 2021,
                             "FISH-SB1-1",
                             fish_site),
         common = str_to_lower(common),
         common = str_replace_all(common, "atlantic", "Atlantic"),
         common = ifelse(common == "silversides", "Atlantic silverside", 
                         ifelse(common == "croaker", "Atlantic croaker",
                                ifelse(common == "american halfbeak", "American halfbeak", common))),
         length_category = NA) %>% 
  dplyr::rename(do = dissolvedOxygen,
                temperature = Temperature,
                depth = waterDepth,
                length = Length,
                site_comment = siteComment) 


# Expand counts such that each row is a different fish
fish2021_updated <- NULL
for (i in 1:nrow(fish_2021)){
  df <- fish_2021 %>% dplyr::slice(i)
  if (!is.na(df$Count)){
    ndf <- df[rep(row.names(df), df$Count), ]
    ndf$Count <- 1
  } else {
    ndf <- df
  }
  assign('fish2021_updated', rbind(fish2021_updated, ndf))
}

# final 2021 output
fish2021_final <- 
  fish2021_updated %>% 
  rename_with(str_to_lower, date:modDate) %>% 
  dplyr::select(-count, -dateid:-moddate) %>% 
  left_join(.,loc %>% 
              dplyr::select(synoptic_reference, synoptic_latitude, 
                            synoptic_longitude, fish_site))


# Bind 2021 and 2019 data----
fish_final <- bind_rows(fish2019_final,
          fish2021_final) %>% 
  mutate(year = year(date),
         depth = ifelse(year == 2019,
                        depth/100,
                        depth))

write.csv(fish_final, file = here::here("output/vcr_fish_sampling.csv"),
          row.names = F)