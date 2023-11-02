# NGA_Microzoo.R
# author: Gwenn Hennon
# email: gmhennon@alaska.edu
# Purpose: pull down NGA Microzooplankton biomass data from Google Drive, clean up the data and calculate summary stats for further analysis

# Load libraries
librarian::shelf(tidyverse, googledrive)

# define file path
path <- "~/Desktop/NGA-LTER/proc_data"

#identify all csv files on the google drive
raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M"), type = "csv") 
#identify the ID of the files I want
mz_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M")) %>%
  dplyr::filter(name %in% c("NGA_signature_microzoop_summary_abundance_carbon_biomass_2011_2022.csv","NGA_signature_microzoop_abundance_2011_2022.csv"))

#download Microzooplankton time series from google drive
for(i in 1:nrow(mz_NGA_ids)){
  googledrive::drive_download(file = mz_NGA_ids[i, ]$id, overwrite = T,
                              path = file.path(path, mz_NGA_ids[i, ]$name))
}

# read csv file on local computer
mz.raw <- read.csv( file = file.path(path,mz_NGA_ids[1, ]$name))

#Date/time not in a standard format so need to tell R how to read it
mz.raw$Date_Time_.UTC. <- gsub("Z","",mz.raw$Date_Time_.UTC.)
mz.raw$Date_Time_.UTC. <- strptime(mz.raw$Date_Time_.UTC.,format = "%Y-%m-%dT%H:%M:%S", tz="GMT") 

#remove rows with NA for a timestamp and for bottle flags
mz <- mz.raw[!is.na(mz.raw$Date_Time_.UTC.),]
mz <- mz[!is.na(mz$Bottle_Number_Flag),]

#check out data to see if it should be log transformed
hist(mz$Total_MZ_Carbon_Biomass_.ug.L.)
hist(log(mz$Total_MZ_Carbon_Biomass_.ug.L.))

#log transform total Microzoo biomass
mz$Log_Tot_MZ_Biomass <- log(mz$Total_MZ_Carbon_Biomass_.ug.L.)

# Compute areal means for each cruise by: shelf and slope stations
#average by shelf v. slope
cruises <- unique(mz$Cruise)
shelf.stations <- c("GAK1", "GAK2", "GAK3", "GAK4", "GAK5","GAK6","GAK7", "GAK8","GAK9")
slope.stations <- c("GAK10","GAK11","GAK12","GAK13","GAK14","GAK15")
mz.area <- as.data.frame(matrix(NA, nrow=length(cruises), ncol=4))
colnames(mz.area) <- c("Cruise", "Date_Time", "Shelf_ave","Slope_ave")
mz.area$Date_Time <- Sys.time()
index = 1
for(c in cruises){
  sub <- subset(mz, mz$Cruise == c)
  sub.shelf <- sub[which(sub$Station %in% shelf.stations),]
  sub.slope <- sub[which(sub$Station %in% slope.stations),]
  mz.area$Cruise[index] <-c
  mz.area$Date_Time[index] <- mean(sub$Date_Time_.UTC., na.rm=T)
  mz.area$Shelf_ave[index] <- mean(sub.shelf$Log_Tot_MZ_Biomass, na.rm=T)
  mz.area$Slope_ave[index] <- mean(sub.slope$Log_Tot_MZ_Biomass, na.rm=T)
  index= index+1
}

#get the month to split into seasons to remove inter-annual signals
Month <-as.numeric(format(mz.area$Date_Time, format="%m"))
#subset by fall and spring
mz.fall <- mz.area[which(Month > 7),]
mz.spring <- mz.area[which(Month < 6),]

#test to look for linear correlations
plot(mz.spring$Date_Time, mz.spring$Shelf_ave)
points(mz.spring$Date_Time, mz.spring$Slope_ave, col="blue")
lm.shelf <-lm(mz.spring$Shelf_ave~mz.spring$Date_Time)
summary(lm.shelf)
lm.slope <-lm(mz.spring$Slope_ave~mz.spring$Date_Time)
summary(lm.slope)

colnames(mz.area) <- c("Cruise", "Date_Time","logBiomass_shelf","logBiomass_slope")

write.csv(mz.area, file= file.path(path,"NGA_Microzoo_arealmean.csv"))
