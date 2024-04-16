################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2024                 #############################
#############          Double Integration          #############################
## Alexandra Cabanelas #########################################################
################################################################################
# still working on this - this script uses same tau for all taxa
## ------------------------------------------ ##################################
#            Packages -----
## ------------------------------------------ ##################################

library(tidyverse)
#library(fasttime) #because regular posixt function takes too long in loop
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

## ------------------------------------------ ##################################
#            Data -----
## ------------------------------------------ ##################################

abu <- read.csv("raw/EcoMon_v3_8_wDateStrata.csv")
################################################################################

# using 10m2 values; 100m3 data is available
names(abu)<-gsub("_10m2","",names(abu)) #to get rid of _10m2 in colnames

# Select taxa to analyze
#change if wanna do more/different taxa
#c.typicus, c.finmarchicus, pseudocalanus, megan, euphkr, clupea harengus, ammodytes
taxa_of_interest <- c("ctyp", "calfin")

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

# select cols of interest
abu1 <- abu %>%
  select(date, month, day, year, all_of(taxa_of_interest), season, Region)

#change data from wide to long
abu_long <- abu1 %>%
  pivot_longer(cols = c(ctyp:calfin), #need to adjust depending on taxa
               names_to = "taxa", values_to = "abundance")


#remove rows with NANS = no sampling for zp or itchyo 
#both gears not always used
abu_long <- abu_long %>% filter(!is.nan(abundance))

################################################################################
# list of regions - for the loop below
# not essential; used it for regions to exclude "outside" 
# can use it to select specific taxa; but we also filtered those out above in 
#line 30

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
      zscore2$date <- as.Date(zscore2$date, format = "%m/%d/%Y")
      
      zscore <- zscore2 %>%
        #mutate(date = as.Date(fastPOSIXct(date))) %>% this messed up the date
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
#write.csv(zscore, "output/zscore_anomalies.csv")
rm(zscore1, zscore2)

ggplot(zscore, aes(x=year, y=Anomaly_yr)) + 
  geom_line() + 
  facet_grid(taxa ~ season+Region)

#zscore_list <- zscore_list[-which(names(zscore_list) == "ammspp_SNE_fall")]

zscore %>%
  group_by(taxa, Region, season) %>%
  summarize(mean_Anomaly_yr = mean(Anomaly_yr, na.rm = TRUE),
            sd_Anomaly_yr = sd(Anomaly_yr, na.rm = TRUE)) %>%
  print(n = 60)

# checking the list output from the loop - making sure it looks right
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

calculateIntegrations = function(data, tau = 365*1, f = function(x) {mean(x, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) #norm
    
    data[[paste0(n,'Int')]] = NA
    data[[paste0(n,'DInt')]] = NA
    
    for (i in 2:nrow(data)) {
      k = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 #secs
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k])
      data[[paste0(n,'DInt')]][i] = f(data[[paste0(n,'Int')]][k])
    }
    
    par(mfrow = c(3,1))
    plot(data[,1], data[[n]], type = 'l', ylab = n, xlab = '') #original signal
    plot(data[,1], data[[paste0(n,'Int')]], type = 'l', ylab = paste0(n, 'Int'), xlab = '') #integration1
    plot(data[,1], data[[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, 'DInt')) #integration2
  }
  data
}

#There are 86,400 seconds in a day (60 seconds/minute * 60 minutes/hour * 24 hours/day).
#multiplying tau by 86,400 gives the equivalent time window in seconds.
################################################################################

################################################################################
## ------------------------------------------ ##################################
#            Driver Data -----
## ------------------------------------------ ##################################
## Drivers - AMO, NAO, AO, GSI       
################################################################################

## ------------------------------------------ ##################################
#            AMO - Atlantic Multidecadal Oscillation  -----
# https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.amo.dat
# 1973-2023

amo <- read.csv("raw/amo_readyDI.csv", header = TRUE)

# convert character to POSIXct
amo$time <- as.POSIXct(amo$time, format = "%d-%m-%y")
#amo$time <- fastPOSIXct(amo$time) #this attaches time to date
#amo$time <- as.Date(amo$time)

amoInt <- calculateIntegrations(amo)
amoInt$time <- as.Date(amoInt$time)
#write.csv(amoInt, "amoIntegrations.csv")
#driver_amo <- amoInt[, c(1,3)]
#int1_amo <- amoInt[, c(1,4)]
#int2_amo <- amoInt[, c(1,5)]
#################


## ------------------------------------------ ##################################
#            NAO - North Atlantic Oscillation  -----

#https://www.ncei.noaa.gov/access/monitoring/nao/
#https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii
# 1973-2023

#nao_month_long <- read.csv("raw/norm.nao.monthly.b5001.current.ascii.LONG.csv")
#nao_month_long$time <- as.Date(paste(nao_month_long$yr, nao_month_long$month, "15", sep = "-"))
#nao_month_long <- nao_month_long %>% filter(yr > 1972)
#nao_month_long <- nao_month_long[, c("time", "nao", "yr", "month")]
#nao_month_long1 <- nao_month_long[, c(1:2)]
#nao_month_long1$time <- as.POSIXct(nao_month_long1$time, format = "%Y-%m-%d")
#nao_month_long1 <- nao_month_long1[-611 ,]

nao_month_long <- read.csv(file.path("raw", "norm.nao.monthly.b5001.current.ascii.LONG.csv")) %>%
  mutate(time = as.Date(paste(yr, month, "15", sep = "-"))) %>%
  filter(yr > 1972 & 
           time < as.Date("2023-11-15")) %>%
  select(time, nao) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d")) 

naoInt <- calculateIntegrations(nao_month_long)
naoInt$time <- as.Date(naoInt$time)
#write.csv(naoInt, "output/naoIntegrations.csv")
#driver_nao <- naoInt[, c(1,3)]
#int1_nao <- naoInt[, c(1,4)]
#int2_nao <- naoInt[, c(1,5)]
#################


## ------------------------------------------ ##################################
#            AO - Arctic Oscillation   -----

# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml
# 1973-2023

#ao <- read.csv("raw/ao.csv")
#ao$date <- as.Date(paste(ao$year, ao$month, "15", sep = "-"))
#ao <- ao %>% filter(year > 1972)
#ao <- ao[, -c(1:2)]
#ao$date <- as.POSIXct(ao$date, format = "%Y-%m-%d")
#ao <- ao[, c("date", "ao")]
#colnames(ao)[1] <- 'time'
#ao <- ao[-611 ,]

ao <- read.csv(file.path("raw", "ao.csv")) %>%
  mutate(time = as.Date(paste(year, month, "15", sep = "-"))) %>% # 
  filter(year > 1972 & 
           time < as.Date("2023-11-15")) %>%
  select(time, ao) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d"))


aoInt <- calculateIntegrations(ao)
aoInt$time <- as.Date(aoInt$time)
#write.csv(aoInt, "output/aoIntegrations.csv")
#driver_ao <- aoInt[, c(1,3)]
#int1_ao <- aoInt[, c(1,4)]
#int2_ao <- aoInt[, c(1,5)]
#################


## ------------------------------------------ ##################################
#            GSI - Gulf Stream Index   -----

# data(package="ecodata"); 
# gsi <- ecodata::gsi
#write.csv(gsi, "raw/gsi.csv") downloaded 26-MAR-2024
# 1954-2023

#gsi <- read.csv("raw/gsi.csv")
#gsi <- gsi %>% filter(Var == "gulf stream index") #there is also western gsi ..? 
#gsi$Date <- as.Date(sprintf("%.2f.01", gsi$Time), format = "%Y.%m.%d")
#gsi <- gsi[, -c(1, 3, 5)]
#gsi$Time <- floor(gsi$Time)
#gsi <- gsi %>% filter(Time > 1972)
#gsi <- gsi[, -1]
#gsi$Date <- as.POSIXct(gsi$Date, format = "%Y-%m-%d")
#gsi <- gsi[, c("Date", "Value")]
#colnames(gsi)[1] <- 'time'
#colnames(gsi)[2] <- 'gsi'

# need to clean this a bit 
gsi <- read.csv(file.path("raw", "gsi.csv")) %>%
  filter(Var == "gulf stream index") %>% #there is also western gsi ..?
  mutate(Date = as.Date(sprintf("%.2f.01", Time), format = "%Y.%m.%d")) %>%
  select(Time, Value, Date) %>% 
  mutate(Time = floor(Time)) %>%
  filter(Time > 1972) %>%
  select(Date, Value) %>%
  mutate(Date = sub("\\d{2}$", "15", Date)) %>%  # change day to 15th
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d")) %>% #class POSIXct
  rename(time = Date, gsi = Value)

gsiInt = calculateIntegrations(gsi)
gsiInt$time <- as.Date(gsiInt$time) #class Date

#write.csv(gsiInt, "output/gsiIntegrations.csv")
#driver_gsi <- gsiInt[, c(1,3)]
#int1_gsi <- gsiInt[, c(1,4)]
#int2_gsi <- gsiInt[, c(1,5)]
################################################################################

################################################################################
## ------------------------------------------ ##################################
#            Correlations -----
## ------------------------------------------ ##################################
# 3) Test correlations between bio TS + each driver TS 


# list to store the correlations for each combination
correlation_list <- list()
p_value_list <- list()

# Iterate over each group (taxa_region_season) in zscore_list
for (group_name in names(zscore_list)) {
  cat("Processing group:", group_name, "\n")
  
  # Extract the abundance data for the current group
  abundance_data <- zscore_list[[group_name]]
  
  # list to store correlations and pvalues for each driver
  group_correlations <- list()
  group_pvalues <- list()
  
  # Iterate over each atmospheric driver
  for (driver_name in c("amo", "nao", "ao","gsi")) {  # drivers
    cat("Calculating correlations for driver:", driver_name, "\n")
    
    # Extract the driver data
    driver_data <- switch(driver_name,
                          "amo" = amoInt,
                          "nao" = naoInt,
                          "ao" = aoInt,
                          "gsi" = gsiInt)
    
    # Extract the required columns for interpolation
    driver_columns <- switch(driver_name,
                             "amo" = c("time", "amoNorm", "amoInt", "amoDInt"),
                             "nao" = c("time", "naoNorm", "naoInt", "naoDInt"),
                             "ao" = c("time", "aoNorm", "aoInt", "aoDInt"),
                             "gsi" = c("time", "gsiNorm", "gsiInt", "gsiDInt"))
    
    driver <- driver_data[, c(1,3)]
    int1 <- driver_data[, c(1,4)]
    int2 <- driver_data[, c(1,5)]
    
    # Interpolate driver values - need to check it makes sense to use abundance data date 
    interpolated_driver <- approx(driver$time, driver[[2]], xout = abundance_data$date)$y
    interpolated_int1 <- approx(int1$time, int1[[2]], xout = abundance_data$date)$y
    interpolated_int2 <- approx(int2$time, int2[[2]], xout = abundance_data$date)$y
    
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
    
    # Store correlations for the current driver
    group_correlations[[driver_name]] <- correlations
    group_pvalues[[driver_name]] <- p_values
  }
  
  # Store correlations for the current group
  correlation_list[[group_name]] <- group_correlations
  p_value_list[[group_name]] <- group_pvalues
}

# view correlation results
correlation_list

# Convert lists of lists to dataframes
correlation_df <- do.call(rbind, lapply(correlation_list, data.frame, stringsAsFactors = FALSE))
p_value_df <- do.call(rbind, lapply(p_value_list, data.frame, stringsAsFactors = FALSE))

# Add group column to correlation and p-value data frames
correlation_df$group <- rownames(correlation_df)
p_value_df$group <- rownames(p_value_df)

p_value_df <- p_value_df %>%
  mutate(Integration = case_when(
    grepl("\\.1$", group) ~ "Direct",
    grepl("\\.2$", group) ~ "1Int",
    grepl("\\.3$", group) ~ "2Int",
    TRUE ~ NA_character_  # Default value if none of the patterns match
  ))

split_group <- str_split(p_value_df$group, "_")

# Extracting taxa, region, and season
p_value_df$taxa <- sapply(split_group, function(x) x[1])
p_value_df$region <- sapply(split_group, function(x) x[2])
p_value_df$season <- sapply(split_group, function(x) x[3])

p_value_df$season <- gsub("\\..*", "", p_value_df$season)

p_value_long <- p_value_df %>%
  pivot_longer(cols = c(amo:gsi), names_to = "driver", values_to = "pvalue")



################################################################################
################################################################################

#                                       Plots

# list to store the final plots for each group and driver
final_plots_list <- list()

#for saving plots.. 
output_dir <- "directory"

# Iterate over each atmospheric driver
for (driver_name in c("amo", "nao", "ao","gsi")) {
  cat("Processing driver:", driver_name, "\n")
  
  # list to store the final plots for each group
  group_plots_list <- list()
  
  # Iterate over each group in zscore_list
  for (group_name in names(zscore_list)) {
    cat("Processing group:", group_name, "\n")
    
    # Extract the abundance data for the current group
    abundance_data <- zscore_list[[group_name]]
    
    # Extract the pre-calculated correlations for the current driver
    correlations <- correlation_list[[group_name]][[driver_name]]
    current_pvals <- p_value_list[[group_name]][[driver_name]]
    
    # Extract the driver data
    driver_data <- switch(driver_name,
                          "amo" = amoInt,
                          "nao" = naoInt,
                          "ao" = aoInt,
                          "gsi" = gsiInt)
    
    # Extract taxa, Region, and season for the current group from zscore_list
    taxa <- abundance_data$taxa  
    Region <- abundance_data$Region  
    season <- abundance_data$season
    
    # three plots for the current driver
    NormPlot <- ggplot() +
      geom_line(data = driver_data, aes(x = time, y = .data[[paste0(driver_name, "Norm")]]), 
                color = "red", size = 1.2) +
      geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr), 
                color = "blue", size = 1.2) +
      annotate("text", x = min(driver_data$time), y = 3,
               label = sprintf("rho == %.4f", correlations[1]),
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 4.5) +
      geom_text(aes(x = min(driver_data$time), y = 2.4,
                    label = paste("P-value:", round(current_pvals[1], 4))),
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
      annotate("text", x = min(driver_data$time), y = 3,
               label = sprintf("rho == %.4f", correlations[2]),
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 4.5) +
      geom_text(aes(x = min(driver_data$time), y = 2.4,
                    label = paste("P-value:", round(current_pvals[2], 4))),
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
      annotate("text", x = min(driver_data$time), y = 3,
               label = sprintf("rho == %.4f", correlations[3]),
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 4.5) +
      geom_text(aes(x = min(driver_data$time), y = 2.4,
                    label = paste("P-value:", round(current_pvals[3], 4))),
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
    
    # Combine the title and plots
    final_plot <- plot_grid(title, combined_plots, ncol = 1, rel_heights = c(0.05, 1))
    
    print(final_plot)
    
    # file name for saving the plot
    file_name <- paste0(driver_name, "_", group_name, ".png")
    
    # Save plot
    ggsave(file.path(output_dir, file_name), final_plot, 
           width = 8, height = 6, dpi = 300)
    
    
    # Store the final plot for the current group in the list
    group_plots_list[[group_name]] <- final_plot
  }
  
  # Store the final plots for all groups for the current driver in the list
  final_plots_list[[driver_name]] <- group_plots_list
}

# Access the plots 
#plot_to_examine <- final_plots_list[["amo"]][["ctyp_MAB_spring"]]
################################################################################
################################################################################

################################################################################
##            Simple correlations - Annual values - no integration            ##
################################################################################
#                                   Slopewater                                 #
#slopewater <- ecodata::slopewater
slopewater <- read.csv("raw/slopewater.csv", header = T)
# this is from ecodata v5.0.0; downloaded 12-MAR-24
# 1964-2022

#Percent total of water type observed in the deep Northeast Channel
# 150-200 meters water depth

# Specifies slopewater source 
# Labrador Slope Water (LSLW) or Warm Slope Water (WSW))
slopeLSLW <- slopewater %>% filter(Var == "LSLW")
slopeWSW <- slopewater %>% filter(Var == "WSW")


################################################################################
#                               Cold Pool Index                                #
#coldpool <- ecodata::cold_pool
#write.csv(coldpool, "raw/coldpool.csv")
coldpool <- read.csv("raw/coldpool.csv", header = T)
# this is from ecodata v5.0.0; downloaded 26-MAR-24
# 1959-2023

# The cold pool is an area of relatively cold bottom water that forms on the 
# in the Mid-Atlantic Bight.

##################################
simple_correlation_list <- list()

# Iterate over each group in zscore_list
for (group_name in names(zscore_list)) {
  cat("Processing group:", group_name, "\n")
  
  # Extract the abundance data for the current group
  abundance_data <- zscore_list[[group_name]]
  
  # Initialize a list to store correlations for each driver
  sim_group_correlations <- list()
  
  # Iterate over each atmospheric driver
  for (driver_name in c("slopew","coldpool")) {  # drivers
    cat("Calculating correlations for driver:", driver_name, "\n")
    
    # Extract the driver data
    driver_data <- switch(driver_name,
                          "slopew" = amoInt,
                          "coldpool" = naoInt)
    
    # Extract the required columns for interpolation
    driver_columns <- switch(driver_name,
                             "amo" = c("time", "amoNorm", "amoInt", "amoDInt"),
                             "nao" = c("time", "naoNorm", "naoInt", "naoDInt"),
                             "ao" = c("time", "aoNorm", "aoInt", "aoDInt"),
                             "gsi" = c("time", "gsiNorm", "gsiInt", "gsiDInt"))
    
    driver <- driver_data[, c(1,3)]
    
    # Interpolate driver values - need to check it makes sense to use abundance data date 
    interpolated_driver <- approx(driver$time, driver[[2]], xout = abundance_data$date)$y
    #interpolated_int1 <- approx(int1$time, int1[[2]], xout = abundance_data$date)$y
    #interpolated_int2 <- approx(int2$time, int2[[2]], xout = abundance_data$date)$y
    
    # Calculate correlations
    scorrelations <- c(
      cor(abundance_data$Anomaly_yr, interpolated_driver, use = "complete.obs"),
      #cor(abundance_data$Anomaly_yr, interpolated_int1, use = "complete.obs"),
      #cor(abundance_data$Anomaly_yr, interpolated_int2, use = "complete.obs")
    )
    
    # Store correlations for the current driver
    sim_group_correlations[[driver_name]] <- scorrelations
  }
  
  # Store correlations for the current group
  simple_correlation_list[[group_name]] <- sim_group_correlations
}

# View the correlation results
simple_correlation_list


### NEED TO ADD LAG ....??????

################################################################################
################################################################################



