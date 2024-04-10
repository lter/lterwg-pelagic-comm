################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2024                 #############################
#############          Double Integration          #############################
## Alexandra Cabanelas #########################################################
################################################################################
# Double INT custom TAU
## Packages
library(tidyverse)
library(fasttime) #because regular posixt function takes too long in loop
library(gridExtra) 
library(cowplot)
#install.packages("devtools")
#remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE) 
#or pak::pkg_install("noaa-edab/ecodata")
library(ecodata) #coldpool and 
#https://github.com/NOAA-EDAB/ecodata
library(listviewer)
library(magrittr) #map_dfr
################################################################################
## Data
abu <- read.csv("raw/EcoMon_v3_8_wDateStrata.csv")

taxa_of_interest <- c("name1", "name2", "name3")

################################################################################
# add season to df
abu <- abu %>%
  mutate(season = case_when(between(month, 3, 5) ~ "spring",
                            between(month, 6, 8) ~ "summer",
                            between(month, 9, 11) ~ "fall",
                            TRUE ~ "winter"))
# add region to df
abu <- abu %>%
  mutate(Region = case_when(region == 1 ~ "MAB",
                            region == 2 ~ "SNE",
                            region == 3 ~ "GB",
                            region == 4 ~ "GOM",
                            TRUE ~ "Outside"))
################################################################################
#select cols of interest - change if wanna do more/different taxa

abu1 <- abu[ , c(11:14,23:25, 75, 79, 208, 235, 298:299)] #this selects ctypicus and cfinmarchicus
abu1 <- abu %>%
  select(year, month, taxa_of_interest)
#throwing in pseudocalanus, megan, euphkr, clupea harengus, ammodytes

# Looking just at 10m2 for now, but can change by selecting
# the other columns with the 100m3 values 

names(abu1)<-gsub("_10m2","",names(abu1)) #to get rid of _10m2 in colnames

#change data from wide to long
abu_long <- abu1 %>%
  pivot_longer(cols = c(ctyp:ammspp),
               names_to = "taxa", values_to = "abundance")

#remove rows with NANS = no sampling for zp or itchyo 
#both gears not always used
abu_long <- abu_long %>% filter(!is.nan(abundance))
################################################################################
# Define a list of taxa and regions - for the loop below
# not essential; used it for regions to exclude "outside" 
# can use it to select specific taxa; but we also filtered those out above in 
#line 34

#taxa_list <- c("ctyp", "calfin")  
region_list <- c("MAB", "SNE", "GOM", "GB")  #excluding "outside"
#season_list <- c("spring","summer","winter","fall")

################################################################################
#######   ----    Double Integration Calculations:
# 1) Biological TS
#     a. Log(x + min/2) transformation
#     b. Spatial average
#     c. Z-score transformation for plotting on same axes
# 2) Driver TS 
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 
# 3) Test correlations between bio TS + each driver TS 
################################################################################
################################################################################
################################################################################
# 1) Biological TS
#     a. Log(x + min/2) transformation
#     b. Spatial average
#     c. Z-score transformation for plotting on same axes

######### Loop to calculate anomalies

zscore_list <- list()

for (taxa in unique(abu_long$taxa)) {
  for (Region in region_list) { #excluding "outside" region 
    for (season in unique(abu_long$season)) {
      cat("Processing:", taxa, "-", Region, "-", season, "\n")
      
      abu_long_subset <- abu_long[abu_long$taxa == taxa & 
                                    abu_long$Region == Region & 
                                    abu_long$season == season, ]
      
      #print(dim(abu_long_subset))
      
      # Step 1a: Data transform
      abu_longTran <- abu_long_subset %>%
        #filter(taxa == taxa, Region == Region, season == season) %>%
        group_by(taxa, Region, season) %>%
        # Find the minimum non-zero value for data transform for each region and season
        mutate(min_nonzero = min(abundance[abundance != 0], na.rm = TRUE)) %>%
        ungroup() %>%
        #to find the log for each year, have to group by yr
        group_by(taxa, Region, season, year) %>%
        # Log transformation - common/base10
        mutate(LogM = log10(abundance + min_nonzero/2)) %>%
        ungroup()
      
      # Step 1b: Spatial average
      # 1b.1: Calculate the mean of LogM values by season, year, taxa, and region
      abu_longMean <- abu_longTran %>%
        group_by(season, year, taxa, Region) %>%
        summarize(mean_LogM = mean(LogM, na.rm = TRUE))
      
      # 1b.2: Calculate the standard deviation of LogM values by season, taxa, and region
      abu_longSD  <- abu_longMean %>%
        group_by(season, taxa, Region) %>% #not year, to get whole ts 
        summarize(Yc_sd1 = sd(mean_LogM, na.rm = TRUE),
                  Yc_mean1 = mean(mean_LogM, na.rm = TRUE)) 
      
      # Merge with abu_longMean to include year information in final df
      abu_longSD <- abu_longSD %>%
        left_join(abu_longMean, by = c("season", "taxa", "Region")) 
      
      # Step 1c: Calculate z-score/anomalies
      zscore1 <- abu_longSD %>%
        mutate(Anomaly_yr = (mean_LogM - Yc_mean1) / Yc_sd1)
      
      # Merge with abu_longTran to include the date column
      zscore2 <- zscore1 %>%
        left_join(abu_longTran %>% select(season, year, date), by = c("season", "year"))
      
      # Convert date column to POSIXct
      #zscore2$date <- fastPOSIXct(zscore2$date) #this attaches time to date
      
      zscore <- zscore2 %>%
        mutate(date = as.Date(fastPOSIXct(date))) %>%
        distinct(season,taxa,Region,year, .keep_all = T)
      
      # Store zscore in the list
      zscore_list[[paste(taxa, Region, season, sep = "_")]] <- zscore
      
      # Convert date column to Date type to have df with all data 
      #abu_longTran$date <- as.Date(abu_longTran$date)
      
      #zscore <- bind_rows(zscore, abu_longTran) %>%
      #  select(!(month:LogM)) #%>%
      #drop_na(Yc_sd1)
    }
  }
}
zscore <- do.call(rbind, zscore_list) 
#write.csv(zscore, "output/zscore_anomalies_moretaxa.csv")
rm(zscore1, zscore2)

zscore_list <- zscore_list[-which(names(zscore_list) == "ammspp_SNE_fall")]

zscore %>%
  group_by(taxa, Region, season) %>%
  summarize(mean_Anomaly_yr = mean(Anomaly_yr, na.rm = TRUE),
            sd_Anomaly_yr = sd(Anomaly_yr, na.rm = TRUE)) %>%
  print(n = 50)

#checking the list output from the loop - making sure it looks right
listviewer::jsonedit(zscore_list)
purrr::map(zscore_list, "name") #to see taxa_region_season
purrr::map(zscore_list, 2) #to check that taxa were properly assigned 
purrr::map(zscore_list, 1) #can also pipe this 

# to extract as tibble
map_dfr(zscore_list, extract, c("season", "taxa", "Region", "Yc_sd1",
                                "Yc_mean1","year","mean_LogM","Anomaly_yr","date"))
#zscore_list %>% pluck("ctyp_MAB_spring") #to access specific 'groups'

###############################################################################
# 2) Driver TS 
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 

#standardize = function(x) {
#  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
#}


calculateIntegrations = function(data, tau, f = function(x) {mean(x, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) #norm
    
    data[[paste0(n,'Int')]] = NA
    data[[paste0(n,'DInt')]] = NA
    
    for (i in 2:nrow(data)) {
      k = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 #secs
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k])
      data[[paste0(n,'DInt')]][i] = f(data[[paste0(n,'Int')]][k])
    }
  }
  data
}

# tau values for each taxa
tau_values <- c(
  "ctyp" = 30,
  "calfin" = 365*1,
  "pseudo" = 90,
  "megan" = 365*2,
  "euphkr" = 365*1,
  "cluhar" = 365*5,
  "ammspp" = 365*4
)


#There are 86,400 seconds in a day (60 seconds/minute * 60 minutes/hour * 24 hours/day).
#multiplying tau by 86,400 gives the equivalent time window in seconds.
################################################################################
################################################################################

## Drivers - AMO, NAO, AO, GSI
#                               Driver Data 
################################################################################

##                  AMO - Atlantic Multidecadal Oscillation                   ##
amo <- read.csv("raw/amo_readyDI.csv", header = TRUE)
# https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.amo.dat
# 1973-2023

# Convert character to POSIXct
amo$time <- as.POSIXct(amo$time, format = "%d-%m-%y")
#amo$time <- fastPOSIXct(amo$time) #this attaches time to date
#amo$time <- as.Date(amo$time)

#amoInt <- calculateIntegrations(amo, tau_values)
#amoInt$time <- as.Date(amoInt$time)
#write.csv(amoInt, "amoIntegrations.csv")
#driver_amo <- amoInt[, c(1,3)]
#int1_amo <- amoInt[, c(1,4)]
#int2_amo <- amoInt[, c(1,5)]
#################

##                            NAO - North Atlantic Oscillation                ##
#https://www.ncei.noaa.gov/access/monitoring/nao/
#https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii
# 1973-2023
nao_month_long <- read.csv("raw/norm.nao.monthly.b5001.current.ascii.LONG.csv")
nao_month_long$time <- as.Date(paste(nao_month_long$yr, nao_month_long$month, "15", sep = "-"))
nao_month_long <- nao_month_long %>% filter(yr > 1972)
nao_month_long <- nao_month_long[, c("time", "nao", "yr", "month")]
nao_month_long1 <- nao_month_long[, c(1:2)]
nao_month_long1$time <- as.POSIXct(nao_month_long1$time, format = "%Y-%m-%d")
nao_month_long1 <- nao_month_long1[-611 ,]

nao_month_long <- read.csv(file.path("raw", "norm.nao.monthly.b5001.current.ascii.LONG.csv")) %>%
  mutate(time = as.Date(paste(yr, month, "15", sep = "-"))) %>%
  filter(yr > 1972) %>%
  select(time, nao, yr, month) %>%
  select(1,2) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d")) %>%
  filter(time =! as.POSIXct("2023-11-15"))


#naoInt <- calculateIntegrations(nao_month_long1)
#naoInt$time <- as.Date(naoInt$time)
#write.csv(naoInt, "output/naoIntegrations.csv")
#driver_nao <- naoInt[, c(1,3)]
#int1_nao <- naoInt[, c(1,4)]
#int2_nao <- naoInt[, c(1,5)]
#################

##                           AO - Arctic Oscillation                          ##
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml
# 1973-2023
ao <- read.csv("raw/ao.csv")
ao$date <- as.Date(paste(ao$year, ao$month, "15", sep = "-"))
ao <- ao %>% filter(year > 1972)
ao <- ao[, -c(1:2)]
ao$date <- as.POSIXct(ao$date, format = "%Y-%m-%d")
ao <- ao[, c("date", "ao")]
colnames(ao)[1] <- 'time'
ao <- ao[-611 ,]

#aoInt <- calculateIntegrations(ao)
#aoInt$time <- as.Date(aoInt$time)
#write.csv(aoInt, "output/aoIntegrations.csv")
#driver_ao <- aoInt[, c(1,3)]
#int1_ao <- aoInt[, c(1,4)]
#int2_ao <- aoInt[, c(1,5)]
#################

##                               GSI - Gulf Stream Index                      ##
# data(package="ecodata"); 
#gsi <- ecodata::gsi
#write.csv(gsi, "raw/gsi.csv") downloaded 26-MAR-2024
# 1954-2023

gsi <- read.csv("raw/gsi.csv")
gsi <- gsi %>% filter(Var == "gulf stream index") #there is also western gsi ..? 
gsi$Date <- as.Date(sprintf("%.2f.01", gsi$Time), format = "%Y.%m.%d")
gsi <- gsi[, -c(1, 3, 5)]
gsi$Time <- floor(gsi$Time)
gsi <- gsi %>% filter(Time > 1972)
gsi <- gsi[, -1]
gsi$Date <- as.POSIXct(gsi$Date, format = "%Y-%m-%d")
gsi <- gsi[, c("Date", "Value")]
colnames(gsi)[1] <- 'time'
colnames(gsi)[2] <- 'gsi'

#gsiInt = calculateIntegrations(gsi)
#gsiInt$time <- as.Date(gsiInt$time)

#write.csv(gsiInt, "output/gsiIntegrations.csv")
#driver_gsi <- gsiInt[, c(1,3)]
#int1_gsi <- gsiInt[, c(1,4)]
#int2_gsi <- gsiInt[, c(1,5)]

# Create a list of data frames, one for each tau value
amoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(amo, tau))
naoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(nao_month_long1, tau))
aoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(ao, tau))
gsiInt_list <- lapply(tau_values, function(tau) calculateIntegrations(gsi, tau))

#have to figure out a way to save these for next scripts...

listviewer::jsonedit(amoInt_list)


#par(mfrow = c(3,1))
#for (n in names(amo)[-1]) {
#  for (tau_name in names(amoInt_list)) {
#    plot(amo[,1], amoInt_list[[tau_name]][[n]], type = 'l', ylab = n, xlab = 'Date') # Plot each integration
#    mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(amo[,1], amoInt_list[[tau_name]][[paste0(n,'Int')]], type = 'l', ylab = paste0(n, ' Int'), xlab = 'Date') # Plot Int
    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(amo[,1], amoInt_list[[tau_name]][[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, ' DInt'), xlab = 'Date') # Plot DInt
    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#  }
#}

#par(mfrow = c(3,1))
#for (n in names(nao_month_long1)[-1]) {
#  for (tau_name in names(naoInt_list)) {
#    plot(nao_month_long1[,1], naoInt_list[[tau_name]][[n]], type = 'l', ylab = n, xlab = 'Date') # Plot each integration
#    mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(nao_month_long1[,1], naoInt_list[[tau_name]][[paste0(n,'Int')]], type = 'l', ylab = paste0(n, ' Int'), xlab = 'Date') # Plot Int
#    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(nao_month_long1[,1], naoInt_list[[tau_name]][[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, ' DInt'), xlab = 'Date') # Plot DInt
#    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#  }
#}

#par(mfrow = c(3,1))
#for (n in names(ao)[-1]) {
#  for (tau_name in names(aoInt_list)) {
#    plot(ao[,1], aoInt_list[[tau_name]][[n]], type = 'l', ylab = n, xlab = 'Date') # Plot each integration
#    mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(ao[,1], aoInt_list[[tau_name]][[paste0(n,'Int')]], type = 'l', ylab = paste0(n, ' Int'), xlab = 'Date') # Plot Int
#    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(ao[,1], aoInt_list[[tau_name]][[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, ' DInt'), xlab = 'Date') # Plot DInt
#    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#  }
#}

#par(mfrow = c(3,1))
#for (n in names(gsi)[-1]) {
#  for (tau_name in names(gsiInt_list)) {
#    plot(gsi[,1], gsiInt_list[[tau_name]][[n]], type = 'l', ylab = n, xlab = 'Date') # Plot each integration
#    mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(gsi[,1], gsiInt_list[[tau_name]][[paste0(n,'Int')]], type = 'l', ylab = paste0(n, ' Int'), xlab = 'Date') # Plot Int
#    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#    plot(gsi[,1], gsiInt_list[[tau_name]][[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, ' DInt'), xlab = 'Date') # Plot DInt
#    #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
#  }
#}


plotData <- function(data, int_list, tau_values) {
  par(mfrow = c(3,1))
  for (n in names(data)[-1]) {
    for (tau_name in names(int_list)) {
      plot(data[,1], int_list[[tau_name]][[n]], type = 'l', ylab = n, xlab = 'Date') # Plot each integration
      #mtext(paste("TAU:", tau_values[tau_name]), side = 3, line = -2, adj = 1) # Add TAU value annotation
      title(paste("TAU:", tau_values[tau_name]), outer = TRUE, line = -1)
      plot(data[,1], int_list[[tau_name]][[paste0(n,'Int')]], type = 'l', ylab = paste0(n, ' Int'), xlab = 'Date') # Plot Int
      plot(data[,1], int_list[[tau_name]][[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, ' DInt'), xlab = 'Date') # Plot DInt
    }
  }
}

# Now you can call this function for each of your datasets
plotData(amo, amoInt_list, tau_values)
plotData(nao_month_long1, naoInt_list, tau_values)
plotData(ao, aoInt_list, tau_values)
plotData(gsi, gsiInt_list, tau_values)


listviewer::jsonedit(amoInt_list)
listviewer::jsonedit(zscore_list)
################################################################################
################################################################################
# 3) Test correlations between bio TS + each driver TS 
#                                   Correlations 
# List to store the correlations for each combination
correlation_list <- list()

# Iterate over each group (taxa_region_season) in zscore_list
for (group_name in names(zscore_list)) {
  cat("Processing group:", group_name, "\n")
  
  # Extract the abundance data for the current group
  abundance_data <- zscore_list[[group_name]]
  
  # Convert abundance_data$date to Date format
  #abundance_data$date <- as.Date(abundance_data$date)
  
  # List to store correlations for each combination of taxa and driver
  group_correlations <- list()
  
  # Iterate over each atmospheric driver
  for (driver_name in c("amo", "nao", "ao", "gsi")) {  # Drivers
    cat("Calculating correlations for driver:", driver_name, "\n")
    
    # Extract the driver data
    driver_data_list <- switch(driver_name,
                               "amo" = amoInt_list,
                               "nao" = naoInt_list,
                               "ao" = aoInt_list,
                               "gsi" = gsiInt_list)
    
    # Iterate over each taxa and its corresponding tau value
    for (taxa_name in names(tau_values)) {
      # Check if the taxa_name matches the group_name
      if (taxa_name == strsplit(group_name, "_")[[1]][1]) {
        # Get the driver data for the current taxa
        driver_data <- driver_data_list[[taxa_name]]
        
        # Convert driver_data$time to Date format
        driver_data$time <- as.Date(driver_data$time)
        
        # Interpolated driver values - need to check it makes sense to use abundance data date 
        interpolated_driver <- approx(driver_data$time, driver_data[[2]], xout = abundance_data$date)$y
        interpolated_int1 <- approx(driver_data$time, driver_data[[3]], xout = abundance_data$date)$y
        interpolated_int2 <- approx(driver_data$time, driver_data[[4]], xout = abundance_data$date)$y
        
        # Calculate correlations
        correlations <- c(
          cor(abundance_data$Anomaly_yr, interpolated_driver, use = "complete.obs"),
          cor(abundance_data$Anomaly_yr, interpolated_int1, use = "complete.obs"),
          cor(abundance_data$Anomaly_yr, interpolated_int2, use = "complete.obs")
        )
        
        p_values <- c(
          cor.test(abundance_data$Anomaly_yr, interpolated_driver)$p.value,
          cor.test(abundance_data$Anomaly_yr, interpolated_int1)$p.value,
          cor.test(abundance_data$Anomaly_yr, interpolated_int2)$p.value
        )
        
        # Store correlations for the current combination of taxa, driver, and integrations
        group_correlations[[paste0(group_name, "_", driver_name, "_", taxa_name)]] <- correlations
      }
    }
  }
  
  # Store correlations for the current group
  correlation_list[[group_name]] <- group_correlations
}


correlation_list

###############################################################################
################################################################################

#                                       Plots

# List to store the final plots for each group and driver
final_plots_list <- list()

output_dir <- "C:/Users/alexc/OneDrive/Documents/MITWHOI/1. RESEARCH/NES-LTER/03_Pelagic Synthesis WG/DoubleIntegration/figures/biodriverTS_dInt_manyTax"

# Iterate over each group in correlation_list
for (group_name in names(correlation_list)) {
  cat("Processing group:", group_name, "\n")
  
  # Extract the abundance data for the current group
  abundance_data <- zscore_list[[group_name]]
  
  # List to store the final plots for each driver
  group_plots_list <- list()

# Iterate over each atmospheric driver
  for (driver_name in c("amo", "nao", "ao", "gsi")) {  # Drivers
  cat("Calculating correlations for driver:", driver_name, "\n")
  
  # Extract the pre-calculated correlations for the current driver
  correlations <- correlation_list[[group_name]]
  
  # Extract the correlation values for the current driver, taxa, and group combination
  #current_correlations <- correlations[[paste0(group_name, "_", driver_name, "_", taxa_name)]]
  current_correlations <- correlations[[paste0(group_name, "_", driver_name, "_", strsplit(group_name, "_")[[1]][1])]]
  
  # Extract the driver data
  driver_data_list <- switch(driver_name,
                             "amo" = amoInt_list,
                             "nao" = naoInt_list,
                             "ao" = aoInt_list,
                             "gsi" = gsiInt_list)
  
  # Get the driver data for the current taxa
  driver_data <- driver_data_list[[taxa_name]]
  
  # Extract taxa, Region, and season for the current group from zscore_list
  taxa <- abundance_data$taxa  
  Region <- abundance_data$Region  
  season <- abundance_data$season
  
  abundance_data$date <- as.POSIXct(abundance_data$date)
  driver_data$time <- as.POSIXct(driver_data$time)
  
  # Three plots for the current driver
  NormPlot <- ggplot() +
    geom_line(data = driver_data, aes(x = time, y = .data[[paste0(driver_name, "Norm")]]), 
              color = "red", size = 1.2) +
    geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr), 
              color = "blue", size = 1.2) +
    geom_text(aes(x = min(driver_data$time), y = 3,
                  label = paste("Correlation:", round(current_correlations[1], 4))),
              hjust = 0, vjust = 1, color = "black", size = 4.5) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank())
  
  Int1Plot <- ggplot() +
    geom_line(data = driver_data, aes(x = time, y = .data[[paste0(driver_name, "Int")]]), 
              color = "red", size = 1.2) +
    geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr), 
              color = "blue", size = 1.2) +
    geom_text(aes(x = min(driver_data$time), y = 3,
                  label = paste("Correlation:", round(current_correlations[2], 4))),
              hjust = 0, vjust = 1, color = "black", size = 4.5) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank())
  
  Int2Plot <- ggplot() +
    geom_line(data = driver_data, aes(x = time, y = .data[[paste0(driver_name, "DInt")]]),
              color = "red", size = 1.2) +
    geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr), 
              color = "blue", size = 1.2) +
    geom_text(aes(x = min(driver_data$time), y = 3,
                  label = paste("Correlation:", round(current_correlations[3], 4))),
              hjust = 0, vjust = 1, color = "black", size = 4.5) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 14, color = "black"))
  
  # Combine the three plots into one grid
  combined_plots <- plot_grid(NormPlot, Int1Plot, Int2Plot, 
                              ncol = 1, align = "v")
  
  # Add plot title 
  title <- ggdraw() +
    draw_label(label = bquote(bolditalic(.(taxa)) ~ 
                                bold(" - ") ~ .(toupper(driver_name)) ~ 
                                bold(" - ") ~ bold(.(Region)) ~ 
                                bold(" - ") ~ bold(.(season))), size = 12) +
    theme(plot.background = element_rect(fill = "white"))
  
  # Combine the title and the plots
  final_plot <- plot_grid(title, combined_plots, ncol = 1, rel_heights = c(0.1, 1))
  
  file_name <- paste0(driver_name, "_", group_name, ".png")
  
  # Save plot
  ggsave(file.path(output_dir, file_name), final_plot, 
         width = 8, height = 6, dpi = 300)
  
  # Store the final plot for the current driver
  group_plots_list[[driver_name]] <- final_plot
  }
  # Store the final plots for the current group
  final_plots_list[[group_name]] <- group_plots_list
}

final_plots_list[["ctyp_MAB_summer"]]
final_plots_list[["calfin_SNE_fall"]]







