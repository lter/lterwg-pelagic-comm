# PAL Trophic Amplification
# Jack Conroy
# jaconroy@ucsc.edu
# 2 Nov 2023
# updated 27 Jan 2025
# Prep data to test for trophic amplification at PAL using 2x2 m zoop biovolume,
# fish biovolume, and surface chl time series

## ------------------------------------------ ##
#            Housekeeping -----
## ------------------------------------------ ##

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive, zoo)

# Set site
site <- "PAL"

# Create necessary sub-folder(s)
dir.create(path = file.path("raw_data"), showWarnings = F)
dir.create(path = file.path("raw_data", site), showWarnings = F)

# Identify raw data files
# For example, here I'm pulling all the PAL csv files from Google Drive
# A new window will pop up asking you to select the appropriate Google Drive account
# For more help, see: https://nceas.github.io/scicomp.github.io/tutorials.html#using-the-googledrive-r-package
raw_PAL_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/14k1iUJv7a7ZXj34Y35Iy9pEKJnO53ONH"),
                                     type = "csv") 

# For each raw data file, download it into its own site folder
for(k in 1:nrow(raw_PAL_ids)){
  
  # Download file (but silence how chatty this function is)
  googledrive::with_drive_quiet(
    googledrive::drive_download(file = raw_PAL_ids[k, ]$id, overwrite = T,
                                path = file.path("raw_data", site, raw_PAL_ids[k, ]$name)) )
  
  # Print success message
  message("Downloaded file ", k, " of ", nrow(raw_PAL_ids))
}


# total zoop biovolume - north

zoopBV <- read.csv(file.path("raw_data", site, "PAL_2m_BV.csv"))

zoopNorth <- zoopBV %>%
  mutate(totalVol = as.numeric(totalVol)) %>% 
  drop_na(totalVol) %>% 
  subset(Rnd100GridLine > 300 & Rnd100GridLine < 700 &
           totalVol > 0)

zoopNorthClean <- aggregate(zoopNorth$totalVol,
                            by = list(Year = zoopNorth$Year,
                                      rnd100GridLine = zoopNorth$Rnd100GridLine,
                                      rnd020GridStation = zoopNorth$Rnd020GridStn),
                            FUN = mean)

colnames(zoopNorthClean)[4] <- "totalBV"

min(zoopNorthClean$totalBV)

hist(zoopNorthClean$totalBV)

zoopNorthClean$logZoopBV <- log10(zoopNorthClean$totalBV)

hist(zoopNorthClean$logZoopBV)

annualNorth_zoopBV <- aggregate(logZoopBV ~ Year, zoopNorthClean, mean)

plot(logZoopBV ~ Year, data = annualNorth_zoopBV, type = "b",
     main = "Total zooplankton biovolume - North")

zoop.run.mean.north <- rollmean(annualNorth_zoopBV$logZoopBV, k = 5)

sd(zoop.run.mean.north)


plot(annualNorth_zoopBV$Year, annualNorth_zoopBV$logZoopBV,
     type = "b",
     main = paste("PAL North total zoop biovolume, std. dev. = ", round(sd(zoop.run.mean.north), digits = 3)),
     ylim = c(0.8, 2.5),
     xlim = c(1990, 2020),
     ylab = "logZoop")

par(new = T)
plot(seq(1995, 2018),
     zoop.run.mean.north,
     type = "l",
     col = "red",
     lwd = 3,
     ylim = c(0.8, 2.5),
     xlim = c(1990, 2020),
     ylab = "",
     xlab = "")

# total zoop biovolume - south

zoopSouth <- zoopBV %>%
  mutate(totalVol = as.numeric(totalVol)) %>% 
  drop_na(totalVol) %>% 
  subset(Rnd100GridLine > 100 & Rnd100GridLine < 400 &
           totalVol > 0)

zoopSouthClean <- aggregate(zoopSouth$totalVol,
                            by = list(Year = zoopSouth$Year,
                                      rnd100GridLine = zoopSouth$Rnd100GridLine,
                                      rnd020GridStation = zoopSouth$Rnd020GridStn),
                            FUN = mean)

colnames(zoopSouthClean)[4] <- "totalBV"

min(zoopSouthClean$totalBV)

hist(zoopSouthClean$totalBV)

zoopSouthClean$logZoopBV <- log10(zoopSouthClean$totalBV)

hist(zoopSouthClean$logZoopBV)

annualSouth_zoopBV <- aggregate(logZoopBV ~ Year, zoopSouthClean, mean)

plot(logZoopBV ~ Year,
     data = annualSouth_zoopBV,

     type = "b",
     main = "Total zooplankton biovolume - South")

zoop.run.mean.south <- rollmean(annualSouth_zoopBV$logZoopBV, k = 5)

sd(zoop.run.mean.south)


plot(logZoopBV ~ Year,
     data = annualSouth_zoopBV,
     type = "b",
     main = paste("PAL South total zoop biovolume, std. dev. = ", round(sd(zoop.run.mean.south), digits = 3)),
     ylim = c(0.8, 2.5),
     xlim = c(1990, 2020),
     ylab = "logZoop")

par(new = T)
plot(zoop.run.mean.south ~ seq(1995, 2018),
     type = "l",
     col = "red",
     lwd = 3,
     ylim = c(0.8, 2.5),
     xlim = c(1990, 2020),
     ylab = "",
     xlab = "")


# fish biovolume - north

fishNorth <- zoopBV %>%
  mutate(fishVol = as.numeric(fishVol)) %>% 
  drop_na(fishVol) %>% 
  subset(Rnd100GridLine > 300 & Rnd100GridLine < 700)

fishNorthClean <- aggregate(fishNorth$fishVol,
                            by = list(Year = fishNorth$Year,
                                      rnd100GridLine = fishNorth$Rnd100GridLine,
                                      rnd020GridStation = fishNorth$Rnd020GridStn),
                            FUN = mean)

colnames(fishNorthClean)[4] <- "fishBV"

min(fishNorthClean$fishBV)

hist(fishNorthClean$fishBV)

nonZeroFishNorthClean <- subset(fishNorthClean$fishBV, fishNorthClean$fishBV > 0)
min(nonZeroFishNorthClean)

fishNorthClean$logFishBV <- log10(fishNorthClean$fishBV + 
                                    min(nonZeroFishNorthClean) / 2)

hist(fishNorthClean$logFishBV)

annualNorth_fishBV <- aggregate(logFishBV ~ Year, fishNorthClean, mean)

plot(logFishBV ~ Year,
     data = annualNorth_fishBV,
     type = "b",
     main = "PAL North fish biovolume")

fish.run.mean.north <- rollmean(annualNorth_fishBV$logFishBV, k = 5)

sd(fish.run.mean.north)


plot(logFishBV ~ Year,
     data = annualNorth_fishBV, type = "b",
     main = paste("PAL North fish biovolume, std. dev. = ", round(sd(fish.run.mean.north), digits = 3)),
     ylim = c(-2, -0.5),
     xlim = c(1990, 2020),
     ylab = "logFish")

par(new = T)
plot(fish.run.mean.north ~ seq(2011, 2018),
     type = "l",
     col = "red",
     lwd = 3,
     ylim = c(-2, -0.5),
     xlim = c(1990, 2020),
     ylab = "",
     xlab = "")


# fish biovolume - south

fishSouth <- zoopBV %>%
  mutate(fishVol = as.numeric(fishVol)) %>% 
  drop_na(fishVol) %>% 
  subset(Rnd100GridLine > 100 & Rnd100GridLine < 400)

fishSouthClean <- aggregate(fishSouth$fishVol,
                            by = list(Year = fishSouth$Year,
                                      rnd100GridLine = fishSouth$Rnd100GridLine,
                                      rnd020GridStation = fishSouth$Rnd020GridStn),
                            FUN = mean)

colnames(fishSouthClean)[4] <- "fishBV"

min(fishSouthClean$fishBV)

hist(fishSouthClean$fishBV)

nonZeroFishSouthClean <- subset(fishSouthClean$fishBV, fishSouthClean$fishBV > 0)
min(nonZeroFishSouthClean)

fishSouthClean$logFishBV <- log10(fishSouthClean$fishBV + 
                                    min(nonZeroFishSouthClean) / 2)

hist(fishSouthClean$logFishBV)

annualSouth_fishBV <- aggregate(logFishBV ~ Year, fishSouthClean, mean)

plot(logFishBV ~ Year,
     data = annualSouth_fishBV,
     type = "b",
     main = "PAL South fish biovolume")

fish.run.mean.south <- rollmean(annualSouth_fishBV$logFishBV, k = 5)

sd(fish.run.mean.south)


plot(logFishBV ~ Year,
     data = annualSouth_fishBV,
     type = "b",
     main = paste("PAL South fish biovolume, std. dev. = ", round(sd(fish.run.mean.south), digits = 3)),
     ylim = c(-2, -0.5),
     xlim = c(1990, 2020),
     ylab = "logFish")

par(new = T)
plot(fish.run.mean.south ~ seq(2011, 2018),
     type = "l",
     col = "red",
     lwd = 3,
     ylim = c(-2, -0.5),
     xlim = c(1990, 2020),
     ylab = "",
     xlab = "")


# chlorophyll a - north

chlNorthBottles <- read_csv(file.path("raw_data", site, "PAL_Chl_Cruise.csv")) %>% 
  dplyr::select(Year, YearEvent, RoundedGridLine, RoundedGridStation,
                Depth, Chlorophyll) %>% 
  na.omit() %>% 
  subset(RoundedGridLine > 300 & RoundedGridLine < 700 &
           RoundedGridStation < 400 & Depth >= 0 & Depth <= 5 &
           Chlorophyll > 0 & Chlorophyll <= 60)


chlNorthDepths <- aggregate(chlNorthBottles$Chlorophyll,
                            by = list(Year = chlNorthBottles$Year,
                                      YearEvent = chlNorthBottles$YearEvent,
                                      RoundedGridLine = chlNorthBottles$RoundedGridLine,
                                      RoundedGridStn = chlNorthBottles$RoundedGridStation,
                                      Depth = chlNorthBottles$Depth),
                            FUN = mean)

colnames(chlNorthDepths)[6] <- "Chl"

chlNorthCasts <- aggregate(chlNorthDepths$Chl,
                           by = list(Year = chlNorthDepths$Year,
                                     YearEvent = chlNorthDepths$YearEvent,
                                     RoundedGridLine = chlNorthDepths$RoundedGridLine,
                                     RoundedGridStn = chlNorthDepths$RoundedGridStn),
                           FUN = mean)

colnames(chlNorthCasts)[5] <- "Chl"

chlNorthStns <- aggregate(chlNorthCasts$Chl,
                          by = list(Year = chlNorthCasts$Year,
                                    RoundedGridLine = chlNorthCasts$RoundedGridLine,
                                    RoundedGridStn = chlNorthCasts$RoundedGridStn),
                          FUN = mean)

colnames(chlNorthStns)[4] <- "Chl"

hist(chlNorthStns$Chl)

chlNorthStns$logChl <- log10(chlNorthStns$Chl)

hist(chlNorthStns$logChl)

AnnualNorth_Chl <- aggregate(logChl ~ Year, chlNorthStns, mean)

plot(logChl ~ Year,
     data = AnnualNorth_Chl,
     pch = 16,
     type = 'b',
     lwd = 3, 
     main = "Surface chlorophyll a - North")

chl.run.mean.north <- rollmean(AnnualNorth_Chl$logChl, k = 5)

sd(chl.run.mean.north)

plot(logChl ~ Year,
     data = AnnualNorth_Chl,
     type = "b",
     main = paste("PAL North chl a, std. dev. = ", round(sd(chl.run.mean.north), digits = 3)),
     ylim = c(-0.5, 1),
     xlim = c(1990, 2020),
     ylab = "logChl")

par(new = T)
plot(chl.run.mean.north ~ seq(1995, 2018),
     type = "l",
     col = "red",
     lwd = 3,
     ylim = c(-0.5, 1),
     xlim = c(1990, 2020),
     ylab = "",
     xlab = "")


# chlorophyll a - south

chlSouthBottles <- read_csv(file.path("raw_data", site, "PAL_Chl_Cruise.csv")) %>% 
  dplyr::select(Year, YearEvent, RoundedGridLine, RoundedGridStation,
                Depth, Chlorophyll) %>% 
  na.omit() %>% 
  subset(RoundedGridLine > 100 & RoundedGridLine < 400 &
           RoundedGridStation < 400 & Depth >= 0 & Depth <= 5 &
           Chlorophyll > 0 & Chlorophyll <= 60)


chlSouthDepths <- aggregate(chlSouthBottles$Chlorophyll,
                            by = list(Year = chlSouthBottles$Year,
                                      YearEvent = chlSouthBottles$YearEvent,
                                      RoundedGridLine = chlSouthBottles$RoundedGridLine,
                                      RoundedGridStn = chlSouthBottles$RoundedGridStation,
                                      Depth = chlSouthBottles$Depth),
                            FUN = mean)

colnames(chlSouthDepths)[6] <- "Chl"

chlSouthCasts <- aggregate(chlSouthDepths$Chl,
                           by = list(Year = chlSouthDepths$Year,
                                     YearEvent = chlSouthDepths$YearEvent,
                                     RoundedGridLine = chlSouthDepths$RoundedGridLine,
                                     RoundedGridStn = chlSouthDepths$RoundedGridStn),
                           FUN = mean)

colnames(chlSouthCasts)[5] <- "Chl"

chlSouthStns <- aggregate(chlSouthCasts$Chl,
                          by = list(Year = chlSouthCasts$Year,
                                    RoundedGridLine = chlSouthCasts$RoundedGridLine,
                                    RoundedGridStn = chlSouthCasts$RoundedGridStn),
                          FUN = mean)

colnames(chlSouthStns)[4] <- "Chl"

hist(chlSouthStns$Chl)

chlSouthStns$logChl <- log10(chlSouthStns$Chl)

hist(chlSouthStns$logChl)

AnnualSouth_Chl <- aggregate(logChl ~ Year, chlSouthStns, mean)

plot(logChl ~ Year,
     data = AnnualSouth_Chl,
     pch = 16,
     type = 'b',
     lwd = 3, 
     main = "Surface chlorophyll a - South")

chl.run.mean.south <- rollmean(AnnualSouth_Chl$logChl, k = 5)

sd(chl.run.mean.south)

plot(logChl ~ Year,
     data = AnnualSouth_Chl,
     type = "b",
     main = paste("PAL South chl a, std. dev. = ", round(sd(chl.run.mean.south), digits = 3)),
     ylim = c(-0.5, 1),
     xlim = c(1990, 2020),
     ylab = "logChl")

par(new = T)
plot(chl.run.mean.south ~ seq(1995, 2018),
     type = "l",
     col = "red",
     lwd = 3,
     ylim = c(-0.5, 1),
     xlim = c(1990, 2020),
     ylab = "",
     xlab = "")


# compile time series to send to Stukel

# palTimeSeries <- as.data.frame(matrix(data = c(1993:2020,
#                                            AnnualNorth_Chl$logChl,
#                                            AnnualSouth_Chl$logChl,
#                                            annualNorth_zoopBV$logZoopBV,
#                                            annualSouth_zoopBV$logZoopBV,
#                                            rep(NA, 16), annualNorth_fishBV$logFishBV,
#                                            rep(NA, 16), annualSouth_fishBV$logFishBV),
#                                   nrow = 28, ncol = 7,
# dimnames = list(NULL, c("Year", "ChlNorth", "ChlSouth",
#                         "ZoopNorth", "ZoopSouth",
#                         "FishNorth", "FishSouth"))))
# 
# write_csv(palTimeSeries, "PAL time series.csv")
