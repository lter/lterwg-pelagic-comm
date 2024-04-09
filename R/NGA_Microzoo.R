# NGA_Microzoo.R
# author: Gwenn Hennon, Bia Dias
# email: gmhennon@alaska.edu
# Purpose: pull down NGA Microzooplankton biomass data from Google Drive, clean up the data and calculate summary stats for further analysis

# Load libraries
librarian::shelf(tidyverse, googledrive, janitor, here)
source('R/functions.R')

# define file path (change for your own purpose)
path <- "~//LTER_WG_pelagic"

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

#Date/time not in a standard format so need to tell R how to read it
mz.raw$date_time_utc <- gsub("Z","",mz.raw$date_time_utc)
mz.raw$date_time_utc <- strptime(mz.raw$date_time_utc,format = "%Y-%m-%dT%H:%M:%S", tz="GMT")

#remove rows with NA for a timestamp and for bottle flags
mz <- mz.raw[!is.na(mz.raw$date_time_utc),]
mz <- mz[!is.na(mz$bottle_number_flag),]

#check out data to see if it should be log transformed
hist(mz$total_mz_carbon_biomass_ug_l)
hist(log(mz$total_mz_carbon_biomass_ug_l))

#log transform total Microzoo biomass
mz$log_tot_mz_biomass <- log(mz$total_mz_carbon_biomass_ug_l)

# Compute areal means for each cruise by: shelf and slope stations
#average by shelf v. slope
cruises <- unique(mz$cruise)
shelf.stations <- c("GAK1", "GAK2", "GAK3", "GAK4", "GAK5","GAK6","GAK7", "GAK8","GAK9")
slope.stations <- c("GAK10","GAK11","GAK12","GAK13","GAK14","GAK15")
mz.area <- as.data.frame(matrix(NA, nrow=length(cruises), ncol=4))
colnames(mz.area) <- c("cruise", "date_time", "shelf_mean","slope_mean")
mz.area$date_time <- Sys.time()
index = 1
for(c in cruises){
  sub <- subset(mz, mz$cruise == c)
  sub.shelf <- sub[which(sub$station %in% shelf.stations),]
  sub.slope <- sub[which(sub$station %in% slope.stations),]
  mz.area$cruise[index] <-c
  mz.area$date_time[index] <- mean(sub$date_time_utc, na.rm=T)
  mz.area$shelf_mean[index] <- mean(sub.shelf$log_tot_mz_biomass, na.rm=T)
  mz.area$slope_mean[index] <- mean(sub.slope$log_tot_mz_biomass, na.rm=T)
  index= index+1
}

#get the month to split into seasons to remove inter-annual signals
month <-as.numeric(format(mz.area$date_time, format="%m"))
#subset by fall and spring
mz.fall <- mz.area[which(month > 7),]
mz.spring <- mz.area[which(month < 6),]

#function for calculating a smooth of parameter with k= window size between 5 and big number (credit Tom Kelly)
sma = function(parameter, k = 100) {
   if (k %% 2 == 1) {
     k = (k - 1) / 2
   } else {
    stop('Please use an odd size window.')
  }
  sma = rep(NA, length(parameter)) ## need this vector to be the same length as parameter
  for (i in 1:length(parameter)) {
    if (i <= k | i >= length(parameter) - k + 1) {
      sma[i] = NA
    } else {
      series = c((i - k):(i + k))
      #message(k, ' = k\t', length(series), ' = series')
      sma[i] = mean(parameter[series], na.rm = TRUE)
    }
  }
  
  sma
}

#perform five year smooth
mz.spring$shelf_5yr_mean <- sma(mz.spring$shelf_mean, k=5)
mz.spring$slope_5yr_mean <- sma(mz.spring$slope_mean, k=5)
mz.fall$shelf_5yr_mean <- sma(mz.fall$shelf_mean, k=5)
mz.fall$slope_5yr_mean <- sma(mz.fall$slope_mean, k=5)

#test to look for linear correlations
lm.shelf <-lm(mz.spring$shelf_mean~mz.spring$date_time)
summary(lm.shelf)
lm.slope <-lm(mz.spring$slope_mean~mz.spring$date_time)
summary(lm.slope)
lm.fall.shelf <-lm(mz.fall$shelf_mean~mz.fall$date_time)
summary(lm.fall.shelf)
lm.fall.slope <-lm(mz.fall$slope_mean~mz.fall$date_time)
summary(lm.fall.slope)

#make plots
plot(
  mz.spring$date_time,
  mz.spring$shelf_mean,
  col = "#FFC20A",
  ylim = c(0, 5),
  xlab = "Date",
  ylab = "mean log10 Microzoop Biomass",
  main = "NGA-LTER GAK Line Spring"
)
points(mz.spring$date_time, mz.spring$slope_mean, col = "#0C7BDC")
lines(mz.spring$date_time, mz.spring$shelf_5yr_mean, col = "#FFC20A")
lines(mz.spring$date_time, mz.spring$slope_5yr_mean, col = "#0C7BDC")
legend(
  "topright",
  c("shelf", "slope"),
  col = c("#FFC20A", "#0C7BDC"),
  lty = 1,
  pch = 1
)
mtext(
  side = 1,
  line = -1,
  adj = 0,
  text = paste("Shelf SD =", round(
    sd(mz.spring$shelf_5yr_mean, na.rm = T), digits = 2
  ))
)
mtext(side = 1,
      line = -1,
      text = paste("Slope SD =", round(
        sd(mz.spring$slope_5yr_mean, na.rm = T), digits = 2
      )))

plot(
  mz.fall$date_time,
  mz.fall$shelf_mean,
  col = "#FFC20A",
  ylim = c(0, 5),
  xlab = "Date",
  ylab = "mean log10 Microzoop Biomass",
  main = "NGA-LTER GAK Line Fall"
)
points(mz.fall$date_time, mz.fall$slope_mean, col = "#0C7BDC")
lines(mz.fall$date_time, mz.fall$shelf_5yr_mean, col = "#FFC20A")
lines(mz.fall$date_time, mz.fall$slope_5yr_mean, col = "#0C7BDC")
legend(
  "topright",
  c("shelf", "slope"),
  col = c("#FFC20A", "#0C7BDC"),
  lty = 1,
  pch = 1
)
mtext(
  side = 1,
  line = -1,
  adj = 0,
  text = paste("Shelf SD =", round(
    sd(mz.fall$shelf_5yr_mean, na.rm = T), digits = 2
  ))
)
mtext(side = 1,
      line = -1,
      text = paste("Slope SD =", round(
        sd(mz.fall$slope_5yr_mean, na.rm = T), digits = 2
      )))

#rename columns for clarity
colnames(mz.area) <- c("Cruise", "Date_Time","logBiomass_shelf","logBiomass_slope")

#write csv file out
write.csv(mz.area, file= file.path(path,"NGA_Microzoo_arealmean.csv"))
write.csv(mz.spring, file=file.path(path,"NGA_Microzoo_Spring.csv"))
write.csv(mz.fall, file=file.path(path,"NGA_Microzoo_Fall.csv"))




