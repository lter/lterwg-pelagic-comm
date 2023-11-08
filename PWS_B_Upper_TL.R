
#install.packages("librarian")
librarian::shelf(tidyverse, googledrive)
library(here)

# Set site
site <- "NGA"

dir.create(path = file.path("raw_data"), showWarnings = F)
dir.create(path = file.path("raw_data", site), showWarnings = F)

raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/1ekR1tuhy-9PXiPaop9gig7jH7RgceYdN",
type = "csv"))

for(k in 1:nrow(raw_NGA_ids)){
  
  # Download file (but silence how chatty this function is)
  googledrive::with_drive_quiet(
    googledrive::drive_download(file = raw_NGA_ids[k, ]$id, overwrite = T,
                                path = file.path("raw_data", site, raw_NGA_ids[k, ]$name)) )
  
  # Print success message
  message("Downloaded file ", k, " of ", nrow(raw_NGA_ids))
}

# Read in a csv
YBPWS <- read.csv(file = file.path("raw_data", site,"Biomass_uppertrophic_PWS_bins_TL.csv"))


BPWS <- as_tibble(BPWS)

BPWSlog <- BPWS %>% mutate(across(c("harborSeal_B":"Herring"), log10))


BPWSlog_long <- tidyr::gather(BPWSlog, key = "Variable", value = "Value", -Year)
cbPalette <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Now, let's create the plot with the improved aesthetics
ggplot(data=BPWSlog_long, aes(x=Year, y=Value)) + 
  geom_smooth(aes(colour=Variable, fill=Variable), method=lm)+
  geom_line(aes(colour= Variable), size=1, alpha=0.2)+
  #geom_point(aes(colour= Variable), size=1)+
  scale_fill_manual(values=cbPalette) + 
  scale_color_manual(values= cbPalette)+# This will add different colors
  theme_minimal() +  # Cleaner theme
  labs(
    title="Upper Trophic Level",
    x="Year",
    y="Log10 Biomass",
    color=""  # Legend title for the colors
  ) + 
  theme(legend.position="bottom") +  # Move legend to the bottom
  guides(size=FALSE, fill= FALSE)


# Function to calculate the simple moving average with proper padding

sma = function(parameter, k = 100) {
  sma = rep(NA, length(parameter)) ## need this vector to be the same length as parameter
  for (i in 1:length(parameter)) {
    if (i < k/2 | i > length(parameter) - k/2) {
      sma[i] = NA
    } else {
      series = c((i - floor(k/2) + 1):(i + floor(k/2) - 1))
      #message(k, ' = k\t', length(series), ' = series')
      sma[i] = mean(parameter[series], na.rm = TRUE)
    }
  }
  
  sma
}

#------------------------------

species <- colnames(BPWSlog[,2:10])

# Rolling standard deviation 
sma_sd <- function(parameter, k = 100) {
  # Calculate rolling standard deviation using the same window size
  roll_sd <- sqrt(stats::filter(parameter^2, rep(1/k, k), sides = 1) - (stats::filter(parameter, rep(1/k, k), sides = 1))^2)
  return(roll_sd)
}

# SMA and SD for each fg
sma_data <- lapply(BPWSlog[species], sma, k = 5)
sma_sd_data <- lapply(BPWSlog[species], sma_sd, k = 5)


# setting the plots in the sheet (row,column)
par(mfrow = c(2, 4))

# 
for (i in seq_along(species)) {
  # Extract the appropriate Year data to match the SMA and SD lengths
  valid_indices <- which(!is.na(sma_data[[i]]))
  adjusted_years <- BPWSlog$Year[valid_indices]
  
  # Plot the SMA with automatically determined y-axis limits
  if(length(valid_indices) > 0) {  # Check if there are valid indices to plot
    plot(adjusted_years, sma_data[[i]][valid_indices], type = 'l', main = species[i],
         xlab = "Year", ylab = "Value")
    
    # Adding a legend with the average standard deviation for the window
    avg_sd <- mean(sma_sd_data[[i]][valid_indices], na.rm = TRUE)
    legend("topleft", legend = paste("SD =", round(avg_sd, 2)), bty = "n")
  } else {
    # Plot empty plot if there are no valid indices
    plot(1, type = 'n', main = species[i], xlab = "Year", ylab = "Value")
    legend("topright", legend = "No data available", bty = "n")
  }
}

#---------------------------------------

# Assuming species is a character vector containing the names of the species columns
species <- c("Halibut_C", "Herring", "Salmon_C", "Gandids", "Groundfish", "Humpbacks_B", "Seabirds", "seaotter_B")

# Calculate SMA for each species and store it in a list
sma_list <- setNames(lapply(BPWSlog[species], sma, k = 5), species)

# Plotting all at once
par(mfrow = c(2, 4))  # Adjust the dimensions as needed for your number of species

# Loop over the species and plot
for (spec in species) {
  valid_indices <- which(!is.na(sma_list[[spec]]))
  adjusted_years <- BPWSlog$Year[valid_indices]
  
  # Plot if there are valid data points to plot
  if (length(valid_indices) > 0) {
    plot(adjusted_years, sma_list[[spec]][valid_indices], type = 'l', main = spec,
         xlab = "Year", ylab = "Value")
    

  } else {
    # Handle the case where there are no valid data points
    plot(1, type = 'n', main = spec, xlab = "Year", ylab = "Value")
    text(1, 0.5, "No data available")
  }
}













