#NGA_chl_synthesis.R
#Clean up Seward line chlorophyll time series
#integrate chlorophyll by station, then average stations on shelf v. over the slope, split by season into fall and spring cruises
# author: Gwenn Hennon
# email: gmhennon@gmail.com

# Housekeeping (borrowed from Angel Chen's scripts)#
# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
install.packages("librarian")
librarian::shelf(tidyverse, googledrive)


# define file path
path <- "~/Desktop/NGA-LTER/proc_data"

#read in chlorophyll time series
#identify all csv files on the google drive
#raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M"), type = "csv") 
#identify the ID of the one file I want
raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M")) %>%
  dplyr::filter(name %in% c("Seward Line chl timeseries dataset 1997-2022.csv"))

#load raw chlorophyll data from google drive and then bring into R
googledrive::drive_download(file = raw_NGA_ids[1, ]$id, overwrite = T,
                                       path = file.path(path, raw_NGA_ids[1, ]$name))
chl.raw <- read.csv( file = file.path(path,raw_NGA_ids[1, ]$name))

#Make a new object and replace 'NR' with NA to make missing values more computer readable
chl = chl.raw
chl[chl == 'NR'] = NA

#remove any rows with NAs in date/time
chl <- chl[!is.na(chl$Date_Time),]

#make sure time is interpreted correctly by R
chl$Date_Time <- strptime(chl$Date_Time, format="%m/%d/%Y %H:%M")


#integrate Total chlorophyll for each cruise by station
#run a nested for loop by cruise
#for each station, sum the chlorophyll by depth
stations <- c("GAK1", "GAK2", "GAK3", "GAK4", "GAK5","GAK6","GAK7", "GAK8", "GAK9", "GAK10", "GAK11","GAK12","GAK13","GAK14","GAK15" )
cruises <- unique(chl$Cruise)
chl.int <- as.data.frame(matrix(NA, nrow=length(cruises)*length(stations), ncol=4))
colnames(chl.int) <- c("Cruise", "Station","Date_Time", "Int_Tot_Chl")
chl.int$Date_Time <- Sys.time()
index = 1
for(i in 1:length(cruises)){
  #print(c)
  c= cruises[i]
  sub <- subset(chl, chl$Cruise == c)
  for(j in 1:length(stations)){
    s= stations[j]
    sub2 <- subset(sub, sub$Station == s)
    #calculate trapezoidal integration of chlorophyll
    z <- sub2$Depth_.m.
    y <- as.numeric(sub2$Total_Chl_A._.ug.L.)
    dz <- diff(z)
    yy <- 0.5*(y[-1] +y[-length(y)])
    Int <- abs(sum(dz*yy))
    chl.int$Cruise[index] <- c
    chl.int$Date_Time[index] <-sub2[1,"Date_Time"]
    chl.int$Station[index] <- s
    chl.int$Int_Tot_Chl[index] <- Int
    index = index + 1
  }
}

#remove missing stations where there are not chl data
chl.int <- chl.int[which(chl.int$Int_Tot_Chl > 0),]

#Log transform the chlorophyll data
hist(log10(chl.int$Int_Tot_Chl))
chl.int$Log_Tot_Chl <- log10(chl.int$Int_Tot_Chl)

#average by shelf v. slope
shelf.stations <- c("GAK1", "GAK2", "GAK3", "GAK4", "GAK5","GAK6","GAK7", "GAK8","GAK9")
slope.stations <- c("GAK10","GAK11","GAK12","GAK13","GAK14","GAK15")
chl.area <- as.data.frame(matrix(NA, nrow=length(cruises), ncol=4))
colnames(chl.area) <- c("Cruise", "Date_Time", "Shelf_ave","Slope_ave")
chl.area$Date_Time <- Sys.time()
index = 1
for(c in cruises){
  sub <- subset(chl.int, chl.int$Cruise == c)
  sub.shelf <- sub[which(sub$Station %in% shelf.stations),]
  sub.slope <- sub[which(sub$Station %in% slope.stations),]
  chl.area$Cruise[index] <-c
  chl.area$Date_Time[index] <- mean(sub$Date_Time, na.rm=T)
  chl.area$Shelf_ave[index] <- mean(sub.shelf$Log_Tot_Chl, na.rm=T)
  chl.area$Slope_ave[index] <- mean(sub.slope$Log_Tot_Chl, na.rm=T)
  index= index+1
}

#get the month to split into seasons to remove inter-annual signals
Month <-as.numeric(format(chl.area$Date_Time, format="%m"))
#subset by fall and spring
chl.fall <- chl.area[which(Month > 7),]
chl.spring <- chl.area[which(Month < 6),]

#test to look for linear correlations
lm.shelf <-lm(chl.spring$Shelf_ave~chl.spring$Date_Time)
lm.slope <-lm(chl.spring$Slope_ave~chl.spring$Date_Time)

#Write outfiles
write.csv(chl.area, file= file.path(path,"Seward_Line_Chl.csv"))
write.csv(chl.spring, file= file.path(path,"Seward_Line_Chl_spring.csv"))
write.csv(chl.fall, file= file.path(path,"Seward_Line_Chl_fall.csv"))
