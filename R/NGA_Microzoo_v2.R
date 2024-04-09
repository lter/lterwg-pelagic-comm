# NGA_Microzoo_v2.R
# author: Bia Dias
# email: biadsdias@gmail.com
# Purpose: pull down NGA Microzooplankton biomass data from Google Drive. Here I am compiling Susanne's data from microzoo and bellow 20 size classes. We will generate raw data size bin plots. 
# Load libraries
librarian::shelf(tidyverse, googledrive, janitor, here)
source('R/functions.R')

# define file path (change for your own purpose)
path <- "C:/Users/Bia/Dropbox/A_GWA/LTER_WG_pelagic/lterwg-pelagic-comm/raw_data/"

#identify all csv files on the google drive
raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk"), type = "csv") 

#identify the ID of the files I want
mz_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk")) %>%
  dplyr::filter(name %in% c("NGA_signature_microzoop_summary_abundance_carbon_biomass_2011_2022.csv","NGA_signature_microzoop_abundance_2011_2022.csv"))

#download Microzooplankton time series from google drive
for(i in 1:nrow(mz_NGA_ids)){
  googledrive::drive_download(file = mz_NGA_ids[i, ]$id, overwrite = T,
                              path = file.path(path, mz_NGA_ids[i, ]$name))
}

# read csv file on local computer
# We have decided to use carbon biomass for the size spectra analysis
mz.raw <- read.csv( file = file.path(path,mz_NGA_ids[2, ]$name)) %>% 
  clean_names()

#-------------------------------------------------------------
#-------------------------------------------------------------
#Susanne data for 20 under size bins. 
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
  dir(path = ,
             recursive = TRUE,
             pattern = "less_") %>%
  map_dfr(read_csv) %>%
  bind_rows() %>% 
  clean_names()

#All files merged
write.csv(mz_under20, "mz_under20.csv", row.names=FALSE)

