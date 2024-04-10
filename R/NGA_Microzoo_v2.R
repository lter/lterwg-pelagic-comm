# NGA_Microzoo_v2.R
# author: Bia Dias
# email: biadsdias@gmail.com
# Purpose: pull down NGA Microzooplankton biomass data from Google Drive. Here I am compiling Susanne's data from microzoo and bellow 20 size classes. We will generate raw data size bin plots. 
# Load libraries
librarian::shelf(tidyverse, googledrive, janitor, here)
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



#---------------------------------------------------------------------------


df <- df_0 %>%
  mutate(new_column = sum(select(contains("10-5"))))

