# NGA_Microzoo_v2.R
# author: Bia Dias
# email: biadsdias@gmail.com
# Purpose: pull down NGA Microzooplankton biomass data from Google Drive. Here I am compiling Susanne's data from microzoo and bellow 20 size classes. We will generate raw data size bin plots. 
# Load libraries
librarian::shelf(scales,tidyverse, googledrive, janitor, here)
source('R/functions.R')

#Geometric means of size bins

size_0_5_2_5 <- c(0.5,2.5)
size_0_5_2_5_mean <- exp(mean(log(size_0_5_2_5)))
size_2_5_5 <- c(2.5,5)
size_2_5_5_mean <- exp(mean(log(size_2_5_5)))
size_5_10 <- c(5,10)
size_5_10_mean <- exp(mean(log(size_5_10)))
size_10_20 <- c(10,20)
size_10_20_mean <- exp(mean(log(size_10_20)))
size_20_40 <- c(20,40)
size_20_40_mean <- exp(mean(log(size_20_40)))
size_40_80 <- c(40,80)
size_40_80_mean <- exp(mean(log(size_40_80)))

size_bins <- c(size_0_5_2_5_mean,size_2_5_5_mean,size_5_10_mean,size_10_20_mean,size_20_40_mean,size_40_80_mean)

size_0_5_20 <- c(0.5,20)
size_0_5_20_mean <- exp(mean(log(size_0_5_20)))

# define file path (change for your own purpose)
path <- "C:/Users/Bia/Dropbox/A_GWA/LTER_WG_pelagic/lterwg-pelagic-comm/raw_data/"

#identify all csv files on the google drive
raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk"), type = "csv") 

#identify the ID of the files I want
mz_NGA_ids <-
  googledrive::drive_ls(
    googledrive::as_id(
      "https://drive.google.com/drive/u/2/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk"
    )
  ) %>%
  dplyr::filter(
    name %in% c(
      "NGA_signature_microzoop_summary_abundance_carbon_biomass_2011_2022.csv",
      "NGA_signature_microzoop_abundance_2011_2022.csv",
      "NGA_signature_microzoop_carbon_biomass_2011_2022.csv",
      "NGA_signature_microzoop_carbon_biomass_2011_2022_v2.csv"
    )
  )

#download Microzooplankton time series from google drive
for(i in 1:nrow(mz_NGA_ids)){
  googledrive::drive_download(file = mz_NGA_ids[i, ]$id, overwrite = T,
                              path = file.path(path, mz_NGA_ids[i, ]$name))
}

# read csv file on local computer
# We have decided to use carbon biomass for the size spectra analysis
mz.raw <- read.csv( file = file.path(path,mz_NGA_ids[1, ]$name)) %>% 
  clean_names()

mz.raw.v2 = mz.raw %>%
  dplyr::mutate(
    year = lubridate::year(date_time_utc),
    month = lubridate::month(date_time_utc)
  )


mz.raw.v3 <- mz.raw.v2 %>%
  mutate("3.162278" = rowSums(select(., contains("20_um")))) %>%
  mutate("28.284271" = rowSums(select(., contains("20_39_um")))) %>%
  mutate("56.568542" = rowSums(select(., contains(
    c("40_59_um", "59_um")
  )))) %>% 
  pivot_longer(
    "28.284271":"56.568542",
    names_to = "size_bins",
    values_to = "carbon_um_ug_l"
  ) %>%
  select(c(
    cruise,
    station,
    date_time_utc,
    year,
    month,
    size_bins,
    carbon_um_ug_l
  ))

colnames(mz.raw.v3) <- c(
  "cruise",
  "station",
  "date_time",
  "year",
  "month",
  "size_bins",
  "carbon_um_ug_l"
)


# Export locally
write.csv(mz.raw.v3, "raw_data/nga_mz_above20.csv", row.names = FALSE)

# Export clean dataset to Drive
googledrive::drive_upload(
  media = "raw_data/nga_mz_above20.csv",
  overwrite = T,
  path = googledrive::as_id(
    "https://drive.google.com/drive/u/2/folders/11uAUqgoRt0aHOtY65z5zkVgUM_4ct51z"
  )
)


#-------------------------------------------------------------
#-------------------------------------------------------------
#Suzanne data for 20 under size bins_epifluorescence microscopy data. 
Micro20under_ids <-
  googledrive::drive_ls(
    googledrive::as_id(
      "https://drive.google.com/drive/u/2/folders/1v376Oi6PLVv672xdPivF8XP_4dSTJDZA"
    ),
    type = "csv"
  )%>%
  dplyr::filter(grepl("biomass", name)
  )

for(i in 1:nrow(Micro20under_ids)) {
  googledrive::drive_download(
    file = Micro20under_ids[i,]$id,
    overwrite = T,
    path = file.path(path, Micro20under_ids[i,]$name)
  )
}

# Use setwd() if your directory 


mz_under20 <-
  dir(path = "raw_data",
      full.names = T,
      pattern = "less_") %>%
  map_dfr(read_csv) %>%
  clean_names()

mz.under20 <- mz_under20 %>%
  dplyr::mutate(year = lubridate::year(date_time),
                month = lubridate::month(date_time)) %>%
  mutate("1.118034" = rowSums(select(., contains(
    c("syn_ug_c_l", "0_2_5_ug")
  )))) %>%
  mutate("3.535534" = rowSums(select(., contains("2_5_5_ug_c_l")))) %>%
  mutate("7.071068" = rowSums(select(., contains("5_10_ug_c_l")))) %>%
  mutate("14.142136" = rowSums(select(., contains(
    c("10_15_ug_c_l", "15_20_ug_c_l")
  )))) %>%
  pivot_longer(
    c("1.118034":"14.142136"),
    names_to = "size_bins",
    values_to = "carbon_um_ug_l"
  ) %>%
  select(c(
    cruise,
    station,
    date_time,
    year,
    month,
    size_bins,
    carbon_um_ug_l
  ))


#All files merged
# Export locally
write.csv(mz.under20, "raw_data/nga_mz_under20.csv", row.names=FALSE)

# Export clean dataset to Drive
googledrive::drive_upload(media = "raw_data/nga_mz_under20.csv", overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/2/folders/11uAUqgoRt0aHOtY65z5zkVgUM_4ct51z"))

#-------------------------------------------------------------
#-------------------------------------------------------------
#Gwenn data for flow cytometry data (0.5-2.5 size bins). 

Microflow_ids <-
  googledrive::drive_ls(
    googledrive::as_id(
      "https://drive.google.com/drive/u/2/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk"
    ),
    type = "csv"
  )%>%
  dplyr::filter(grepl("abd_biomass", name)
  )

for(i in 1:nrow(Microflow_ids)) {
  googledrive::drive_download(
    file = Microflow_ids[i,]$id,
    overwrite = T,
    path = file.path(path, Microflow_ids[i,]$name)
  )
}

# Use setwd() if your directory 

mz_flow2019 <- read.csv("raw_data/nga_SKQ2019-15S_FCM_abd_biomass.csv") %>%
  mutate(Date_Time = as.character(Date_Time))
mz_flow2020_1 <- read.csv("raw_data/nga_SKQ2020-10S_FCM_abd_biomass.csv")%>%
  mutate(Date_Time = as.character(Date_Time))
mz_flow2020_2 <- read.csv("raw_data/nga_SKQ2020-12S_FCM_abd_biomass.csv")%>%
  mutate(Date_Time = as.character(Date_Time))
mz_flow2021_1 <- read.csv("raw_data/nga_SKQ2021-06S_FCM_abd_biomass.csv") %>% 
  mutate(Date_Time= stringr::str_replace(Date_Time, "[:space:][:digit:]{1,}:[:digit:]{1,}", "")) %>%
  separate_wider_delim(Date_Time, delim = "/", names = c("month", "day", "year")) %>%
  mutate(year = paste0("20", year),
         Date_Time = as.character(as.Date(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"))) %>%
  relocate(Date_Time, .after = Station) %>%
  select(-month, -day, -year)
mz_flow2021_2 <- read.csv("raw_data/nga_SKQ2021-10S_FCM_abd_biomass.csv")%>%
  mutate(Date_Time = as.character(Date_Time))
mz_flow2021_3 <- read.csv("raw_data/nga_TGX2021-09_FCM_abd_biomass.csv")%>% 
  mutate(Date_Time= stringr::str_replace(Date_Time, "[:space:][:digit:]{1,}:[:digit:]{1,}", "")) %>%
  separate_wider_delim(Date_Time, delim = "/", names = c("month", "day", "year")) %>%
  mutate(year = paste0("20", year),
         Date_Time = as.character(as.Date(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"))) %>%
  relocate(Date_Time, .after = Station) %>%
  select(-month, -day, -year)


mz_flow <-
  bind_rows(
    mz_flow2019,
    mz_flow2020_1,
    mz_flow2020_2,
    mz_flow2021_1,
    mz_flow2021_2,
    mz_flow2021_3
  ) %>%
  clean_names()

mz.flow <- mz_flow %>%
  dplyr::mutate(year = lubridate::year(as.Date(date_time)),
                month = lubridate::month(as.Date(date_time))) %>%
  mutate("1.118034" = rowSums(select(., contains("ug_c_l")))) %>%
  pivot_longer(
    c("1.118034"),
    names_to = "size_bins",
    values_to = "carbon_um_ug_l"
  ) %>%
  select(c(
    cruise,
    station,
    date_time,
    year,
    month,
    size_bins,
    carbon_um_ug_l
  ))


#---------------------------------------------------------------------------
# Export locally
write.csv(mz.flow, "raw_data/nga_mz_flow.csv", row.names=FALSE)

# Export clean dataset to Drive
googledrive::drive_upload(media = "raw_data/nga_mz_flow.csv", overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/2/folders/11uAUqgoRt0aHOtY65z5zkVgUM_4ct51z"))



#----------------------------------------------------------------------------
nga_mz_id <-
  googledrive::drive_ls(
    googledrive::as_id(
      "https://drive.google.com/drive/u/2/folders/11uAUqgoRt0aHOtY65z5zkVgUM_4ct51z"
    ),
    type = "csv"
  )

for(i in 1:nrow(nga_mz_id)) {
  googledrive::drive_download(
    file = nga_mz_id[i,]$id,
    overwrite = T,
    path = file.path(path, nga_mz_id[i,]$name)
  )
}

nga_mz <-
  dir(path = "raw_data",
      full.names = T,
      pattern = "nga_mz_") %>%
  map_dfr(read_csv) %>%
  clean_names() 

nga_mz_v1 <- nga_mz %>% 
  mutate(date_time= stringr::str_replace(date_time, "[:space:][:digit:]{1,}:[:digit:]{1,}:[:digit:]{1,}", "")) 

  
nga_mz_v1 <- nga_mz_v1[
  with(nga_mz_v1, order(date_time, size_bins)),
]

# We needed to assign bin widths to normalization of the size spectra. 
# Here I make binkey date-frame with the size_bins from nga_mz_v1 file, the nga_mz_v1$size_bins cbind to the correspondent bin width.
binkey <- cbind(unique(nga_mz_v1$size_bins), c(20,40,2,2.5,5,10))
binkey <- as.data.frame(binkey)
colnames(binkey) <- c("size_bins","bin_width")

#left join to create the bin width for the nga_mz_v2 file.
nga_mz_v2 <- left_join(nga_mz_v1, binkey)

#Normalizing the size spectra bins with normalized biomass size spectra (NBSS). 
nga_mz_v2 <- nga_mz_v2 %>% 
  mutate(NBSS=carbon_um_ug_l/bin_width)

# Jack function for creating breaks. 
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)}}
    
#ggplot the data. 
#for future reference we will need to have only the data from 2018-present
#for future, separating the 4 domains (coast-GAK1, Shelf- GAK1-GAK9, Oceanic GAK10-GAK15, PWS)

nga_mz_plot_bystation <- ggplot(data = nga_mz_v2, aes(
  x = size_bins,
  y = NBSS,
  group = interaction(date_time, station)
)) +
  geom_point() +
  geom_path() +
  scale_x_continuous(
    trans = "log10",
    breaks = base_breaks(),
    limits = c(5 * 10 ^ -1, 10 ^ 5),
    labels = scales::math_format(.x)
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = base_breaks(),
    limits = c(10 ^ -4, 10 ^ 4),
    labels = scales::math_format(.x)
  ) +
  facet_wrap( ~ year+station) +
  theme_classic() +
  xlab(expression(paste("Equivalent spherical diameter (", mu, "m)"))) +
  ylab(expression(
    paste(
      "Normalized biomass size spectrum (mg C m" ^ "-3",
      " ",
      mu,
      "m" ^ "-1",
      ")"
    )
  ))

ggsave("NGA_FIGURES/nga_mz_plot.pdf",nga_mz_plot, width = 10, height = 4.6, bg="transparent")

ggsave("NGA_FIGURES/nga_mz_plot_bystation.pdf",nga_mz_plot_bystation, width = 6, height = 10, bg="transparent") #work on this to print better!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# saving the final data

# Export locally
write.csv(nga_mz_v2, "raw_data/nga_mz_v2.csv", row.names=FALSE)

# Export clean dataset to Drive
googledrive::drive_upload(media = "raw_data/nga_mz_v2.csv", overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/2/folders/1VRTbP3VP6eM65m5gp4tv2QkZ-4sqv75Q"))

googledrive::drive_upload(media = "NGA_FIGURES/nga_mz_plot_bystation.pdf", overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/2/folders/1R-3mzYjNWzs2LqQtde3sdEdNfagjTuqz"))

googledrive::drive_upload(media = "NGA_FIGURES/nga_mz_plot.pdf", overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/2/folders/1R-3mzYjNWzs2LqQtde3sdEdNfagjTuqz"))
