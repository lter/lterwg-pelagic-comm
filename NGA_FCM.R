# NGA_FCM.R
# author: Gwenn Hennon
# email: gmhennon@alaska.edu
# pull down flow cytometry data from NGA raw files on google drive, clean up the tables (remove missing values etc..), calculate biomass from abundance based on literature values of C per cell, integrate by depth?, average integrated picoplankton biomass by region and season, make plots to visualize, write output .csv files with average picoplankton biomass and upload them to google drive. 

# Load libraries
librarian::shelf(tidyverse, googledrive)

# define file path
path <- "~/Desktop/NGA-LTER/proc_data"

#identify all csv files on the google drive
raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk"), type = "csv") 
#identify the ID of the files I want
fcm_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/10KgTVkAMODzgvfaf6c-XI-eelV_PtUSk")) %>%
  dplyr::filter(name %in% c("NGA-LTER_FCM_abundance_2019_2021.csv"))

#download FCM time series from google drive
for(i in 1:nrow(fcm_NGA_ids)){
  googledrive::drive_download(file = fcm_NGA_ids[i, ]$id, overwrite = T,
                              path = file.path(path, fcm_NGA_ids[i, ]$name))
}

# read csv file on local computer
fcm.raw <- read.csv( file = file.path(path,fcm_NGA_ids[1, ]$name))

#Date/time not in a standard format so need to tell R how to read it
fcm.raw$Date_Time <- strptime(fcm.raw$date,format = "%m/%d/%y") 

#Remove rows with NAs? currently this data set doesn't have any, but may be necessary in future?

#check out the data
#make boxplots of cell conc. data
boxplot(cbind(fcm.raw$picoeukaryotes..cells.mL.,fcm.raw$Synechococcus..cell.mL., fcm.raw$heterotrophic.bacteria..cells.mL.), names=c("Picoeuk","Synecho", "Hetbac") , ylab="cells/mL", xlab="Picoplankton",)

# Calculate biomass of picoplankton population in each sample
# Picoeukaryotes average