# Jack Conroy
# 2 Nov 2023
# PAL 2x2 m zoop biovolume time series
# PAL surface chl time series

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

# total zoop biovolume

zoopBV <- read.csv(file.path("raw_data", site, "PAL_2m_BV.csv"))

keepCols <- colnames(zoopBV)[c(1, 5, 7, 62:65)]
taxaCols <- colnames(zoopBV)[c(62:65)]

bvOnly <- zoopBV %>%
  select(all_of(keepCols)) 

total <- zoopBV %>% drop_na(totalVol)

total$totalVol <- as.numeric(total$totalVol)

totalClean <- aggregate(total$totalVol, by = list(Year = total$Year,
                                              rnd100GridLine = total$Rnd100GridLine,
                                              rnd020GridStation = total$Rnd020GridStn),
                       FUN = mean)

colnames(totalClean)[4] <- "totalBV"

totalClean <- na.omit(totalClean)

min(totalClean$totalBV)

nonZeroTotalClean <- subset(totalClean$totalBV, totalClean$totalBV > 0)
min(nonZeroTotalClean)

totalClean$adjTotalBV <- log10(totalClean$totalBV + min(nonZeroTotalClean) / 2)
hist(totalClean$adjTotalBV)

annualTotalBV <- aggregate(adjTotalBV ~ Year, totalClean, mean)

plot(adjTotalBV ~ Year, data = annualTotalBV, type = "b",
     main = "Total zooplankton biovolume - full grid")

zoop.run.mean <- rollmean(annualTotalBV$adjTotalBV, k = 5)

plot(adjTotalBV ~ Year, data = annualTotalBV, type = "b",
     main = "Total zooplankton biovolume - full grid",
     ylim = c(1, 2.5), xlim = c(1990, 2020))
par(new = T)
plot(zoop.run.mean ~ seq(1995, 2018), type = "l",
     col = "red", lwd = 3,
     ylim = c(1, 2.5), xlim = c(1990, 2020),
     ylab = "", xlab = "")

sd(zoop.run.mean)


# fish biovolume

fish <- zoopBV %>% drop_na(fishVol)

fish$fishVol <- as.numeric(fish$fishVol)

fishClean <- aggregate(fish$fishVol, by = list(Year = fish$Year,
                                               rnd100GridLine = fish$Rnd100GridLine,
                                               rnd020GridStation = fish$Rnd020GridStn),
                       FUN = mean)

colnames(fishClean)[4] <- "fishBV"

fishClean <- na.omit(fishClean)

min(fishClean$fishBV)

nonZeroFishClean <- subset(fishClean$fishBV, fishClean$fishBV > 0)
min(nonZeroFishClean)

fishClean$adjFishBV <- log10(fishClean$fishBV + min(nonZeroFishClean) / 2)

annualFishBV <- aggregate(adjFishBV ~ Year, fishClean, mean)

plot(adjFishBV ~ Year, data = annualFishBV, type = "b",
     main = "Fish biovolume - full grid")

fish.run.mean <- rollmean(annualFishBV$adjFishBV, k = 5)

plot(adjFishBV ~ Year, data = annualFishBV, type = "b",
     main = "Fish biovolume - full grid",
     ylim = c(-2, -0.5), xlim = c(1990, 2020))
par(new = T)
plot(fish.run.mean ~ seq(2011, 2018), type = "l",
     col = "red", lwd = 3,
     ylim = c(-2, -0.5), xlim = c(1990, 2020),
     ylab = "", xlab = "")

sd(fish.run.mean)


# chlorophyll a

chl <- read.csv(file.path("raw_data", site, "chl_log10_gridavg.csv"))

plot(Chlorophyll_log10_mean ~ Year, data = chl, type = "b",
     main = "Surface chlorophyll a - full grid")
 
chl.run.mean <- rollmean(chl$Chlorophyll_log10_mean, k = 5)

plot(Chlorophyll_log10_mean ~ Year, data = chl, type = "b",
     main = "Surface chlorophyll a - full grid",
     ylim = c(-0.5, 1), xlim = c(1990, 2020))
par(new = T)
plot(chl.run.mean ~ seq(1995, 2018), type = "l",
     col = "red", lwd = 3,
     ylim = c(-0.5, 1), xlim = c(1990, 2020),
     ylab = "", xlab = "")

sd(chl.run.mean)

