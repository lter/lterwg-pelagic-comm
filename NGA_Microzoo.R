# NGA_Microzoo.R
# author: Gwenn Hennon
# email: gmhennon@alaska.edu
# Purpose: pull down NGA Microzooplankton biomass data from Google Drive, clean up the data and calculate summary stats for further analysis

# Load libraries
librarian::shelf(tidyverse, googledrive)

# define file path
path <- "~/Desktop/NGA-LTER/proc_data"

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

#function for calculating a smooth of parameter with k= window size between 5 and big number (credit Tom Kelly)
sma = function(parameter, k = 100) {
  # if (k %% 2 == 0) {
  #   k = k - 1
  # }
  sma = rep(NA, length(parameter)) ## need this vector to be the same length as parameter
  for (i in 1:length(parameter)) {
    if (i < k/2 | i > length(parameter) - k/2) {
      sma[i] = NA
    } else {
      series = c((i - floor(k/2)):(i + floor(k/2)))
      #message(k, ' = k\t', length(series), ' = series')
      sma[i] = mean(parameter[series], na.rm = TRUE)
    }
  }
  
  sma
}

#perform five year smooth
mz.spring$Shelf_5yr_mean <- sma(mz.spring$Shelf_ave, k=5)
mz.spring$Slope_5yr_mean <- sma(mz.spring$Slope_ave, k=5)
mz.fall$Shelf_5yr_mean <- sma(mz.fall$Shelf_ave, k=5)
mz.fall$Slope_5yr_mean <- sma(mz.fall$Slope_ave, k=5)

#test to look for linear correlations
lm.shelf <-lm(mz.spring$Shelf_ave~mz.spring$Date_Time)
summary(lm.shelf)
lm.slope <-lm(mz.spring$Slope_ave~mz.spring$Date_Time)
summary(lm.slope)

#make plots
plot(mz.spring$Date_Time, mz.spring$Shelf_ave, col="green", ylim=c(0,5), xlab="Date", ylab="mean log10 Microzoop Biomass", main="NGA-LTER GAK Line Microzoops")
points(mz.spring$Date_Time, mz.spring$Slope_ave, col="blue")
lines(mz.spring$Date_Time, mz.spring$Shelf_5yr_mean, col="green")
lines(mz.spring$Date_Time, mz.spring$Slope_5yr_mean, col="blue")
legend("topright", c("shelf", "slope"), col=c("green","blue"), lty=1, pch=1)
mtext(side=1,line=-1, adj=0, text= paste( "Shelf SD =",round(sd(mz.spring$Shelf_5yr_mean, na.rm=T), digits=2)))
mtext(side=1,line=-1, text= paste( "Slope SD =",round(sd(mz.spring$Slope_5yr_mean, na.rm=T), digits=2)))

#rename columns for clarity
colnames(mz.area) <- c("Cruise", "Date_Time","logBiomass_shelf","logBiomass_slope")

#write csv file out
write.csv(mz.area, file= file.path(path,"NGA_Microzoo_arealmean.csv"))

