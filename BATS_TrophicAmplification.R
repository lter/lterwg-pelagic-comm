# Jack Conroy
# 3 Nov 2023
# BATS zoop biomass time series
# BATS surface chl time series

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

# total zoop dry weight

zoopDW <- read.csv(file.path("raw_data", site, "BATSZoopBiomassMonthly.csv"))

min(zoopDW$avgDryWtMgM3TotYM)
min(na.omit(zoopDW$avgDryWtMgM3TotYM))

zoopDW$adjTotalDW <- log10(zoopDW$avgDryWtMgM3TotYM)
hist(zoopDW$adjTotalDW)

zoopDW$Date <- as.yearmon(paste(zoopDW$Year, zoopDW$Month), "%Y %m")

plot(adjTotalDW ~ Date, data = zoopDW, type = "b",
     main = "Total zooplankton dry weight - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.total <- rollapply(zoopDW$adjTotalDW, 12, function(x) mean(x, na.rm = T))

plot(adjTotalDW ~ Date, data = zoopDW, type = "b",
     main = "Total zooplankton dry weight - BATS",
     ylim = c(-0.2, 1), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.total ~ zoopDW$Date[6:345], type = "l",
     col = "red", lwd = 3,
     ylim = c(-0.2, 1), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.total)


# 0.2-0.5 mm zoop dry weight

min(zoopDW$avgDryWtMgM30200YM)
min(na.omit(zoopDW$avgDryWtMgM30200YM))

zoopDW$adj200DW <- log10(zoopDW$avgDryWtMgM30200YM)
hist(zoopDW$adj200DW)

plot(adj200DW ~ Date, data = zoopDW, type = "b",
     main = "0.2-0.5 mm zooplankton dry weight - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.200 <- rollapply(zoopDW$adj200DW, 12, function(x) mean(x, na.rm = T))

plot(adj200DW ~ Date, data = zoopDW, type = "b",
     main = "0.2-0.5 mm zooplankton dry weight - BATS",
     ylim = c(-1, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.200 ~ zoopDW$Date[6:345], type = "l",
     col = "red", lwd = 3,
     ylim = c(-1, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.200)


# 0.5-1 mm zoop dry weight

min(zoopDW$avgDryWtMgM30500YM)
min(na.omit(zoopDW$avgDryWtMgM30500YM))

zoopDW$adj500DW <- log10(zoopDW$avgDryWtMgM30500YM)
hist(zoopDW$adj500DW)

plot(adj500DW ~ Date, data = zoopDW, type = "b",
     main = "0.5-1 mm zooplankton dry weight - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.500 <- rollapply(zoopDW$adj500DW, 12, function(x) mean(x, na.rm = T))

plot(adj500DW ~ Date, data = zoopDW, type = "b",
     main = "0.5-1 mm zooplankton dry weight - BATS",
     ylim = c(-1, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.500 ~ zoopDW$Date[6:345], type = "l",
     col = "red", lwd = 3,
     ylim = c(-1, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.500)


# 1-2 mm zoop dry weight

min(zoopDW$avgDryWtMgM31000YM)
min(na.omit(zoopDW$avgDryWtMgM31000YM))

zoopDW$adj1000DW <- log10(zoopDW$avgDryWtMgM31000YM)
hist(zoopDW$adj1000DW)

plot(adj1000DW ~ Date, data = zoopDW, type = "b",
     main = "1-2 mm zooplankton dry weight - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.1000 <- rollapply(zoopDW$adj1000DW, 12, function(x) mean(x, na.rm = T))

plot(adj1000DW ~ Date, data = zoopDW, type = "b",
     main = "1-2 mm zooplankton dry weight - BATS",
     ylim = c(-1, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.1000 ~ zoopDW$Date[6:345], type = "l",
     col = "red", lwd = 3,
     ylim = c(-1, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.1000)


# 2-5 mm zoop dry weight

min(zoopDW$avgDryWtMgM32000YM)
min(na.omit(zoopDW$avgDryWtMgM32000YM))

zoopDW$adj2000DW <- log10(zoopDW$avgDryWtMgM32000YM)
hist(zoopDW$adj2000DW)

plot(adj2000DW ~ Date, data = zoopDW, type = "b",
     main = "2-5 mm zooplankton dry weight - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.2000 <- rollapply(zoopDW$adj2000DW, 12, function(x) mean(x, na.rm = T))

plot(adj2000DW ~ Date, data = zoopDW, type = "b",
     main = "2-5 mm zooplankton dry weight - BATS",
     ylim = c(-1.5, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.2000 ~ zoopDW$Date[6:345], type = "l",
     col = "red", lwd = 3,
     ylim = c(-1.5, 0.5), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.2000)


# >5 mm zoop dry weight

min(zoopDW$avgDryWtMgM35000YM)
min(na.omit(zoopDW$avgDryWtMgM35000YM))

zoopDW$adj5000DW <- log10(zoopDW$avgDryWtMgM35000YM)
hist(zoopDW$adj5000DW)

plot(adj5000DW ~ Date, data = zoopDW, type = "b",
     main = ">5 mm zooplankton dry weight - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.5000 <- rollapply(zoopDW$adj5000DW, 12, function(x) mean(x, na.rm = T))

plot(adj5000DW ~ Date, data = zoopDW, type = "b",
     main = ">5 mm zooplankton dry weight - BATS",
     ylim = c(-2, 1), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.5000 ~ zoopDW$Date[6:345], type = "l",
     col = "red", lwd = 3,
     ylim = c(-2, 1), xlim = c(zoopDW$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.5000)


# Surface chlorophyll a

chl <- read.csv(file.path("raw_data", site, "BATSChlAMonthly.csv"))

min(chl$avgChlA)
min(na.omit(chl$avgChlA))

# sometimes the chl concentration is zero. that must be wrong. so replacing with NAs
chl$avgChlA[chl$avgChlA == 0] <- NA
min(na.omit(chl$avgChlA))

chl$adjChl <- log10(chl$avgChlA)
hist(chl$adjChl)

chl$Date <- as.yearmon(paste(chl$Year, chl$Month), "%Y %m")

plot(adjChl ~ Date, data = chl, type = "b",
     main = "Chlorophyll a - BATS")

# I'm not quite sure how this rolling mean is being aligned, because it is even number 
# Can't use rollmean b/c there are NAs
run.mean.chl <- rollapply(chl$adjChl, 12, function(x) mean(x, na.rm = T))

plot(adjChl ~ Date, data = chl, type = "b",
     main = "Chlorophyll a - BATS",
     ylim = c(1, 3), xlim = c(chl$Date[1], zoopDW$Date[length(zoopDW$Date)]))
par(new = T)
# manually coding the dates for the rolling mean b/c uncertain about above alignment
plot(run.mean.chl ~ chl$Date[6:391], type = "l",
     col = "red", lwd = 3,
     ylim = c(1, 3), xlim = c(chl$Date[1], zoopDW$Date[length(zoopDW$Date)]),
     ylab = "", xlab = "")

sd(run.mean.chl)
