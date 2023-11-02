## ------------------------------------------ ##
#        Pelagic Community Structure 
#         Chlorophyll Extraction         
## ------------------------------------------ ##
# Script author(s): Angel Chen

# Purpose: Extracts the chlorophyll values that are inside the boundaries of the GAK shapefile

## ------------------------------------------ ##
#            Housekeeping -----
## ------------------------------------------ ##

# Load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive, sf, terra)

# Set site
site <- "NGA"

# Create necessary sub-folder(s)
dir.create(path = file.path("raw_data"), showWarnings = F)
dir.create(path = file.path("raw_data", site), showWarnings = F)

## -------------------------------------------- ##
#             Data Acquisition ----
## -------------------------------------------- ##

# Identify raw chlorophyll data file
raw_NGA_ids1 <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M"),
                                    pattern = "Chlo_GAK") 

# Identify raw GAK spatial files
raw_NGA_ids2 <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M"),
                                    pattern = "GAK_conc") 

# Combine the file IDs together
raw_ids <- rbind(raw_NGA_ids1, raw_NGA_ids2)

# For each raw data file, download it into its own site folder
for(k in 1:nrow(raw_ids)){
  
  # Download file (but silence how chatty this function is)
  googledrive::with_drive_quiet(
    googledrive::drive_download(file = raw_ids[k, ]$id, overwrite = T,
                                path = file.path("raw_data", site, raw_ids[k, ]$name)) )
  
  # Print success message
  message("Downloaded file ", k, " of ", nrow(raw_ids))
}

## ------------------------------------------ ##
#            Data Extraction -----
## ------------------------------------------ ##

# Read in shapefile
GAK_shapefile <- sf::st_read(dsn = file.path("raw_data", site, "GAK_conc.shp"))

# Plot shapefile
plot(GAK_shapefile, axes = T)

# Read in chlorophyll csv
chloro <- read.csv(file = file.path("raw_data", site, "Chlo_GAK.csv")) %>%
  # Move the time column
  dplyr::relocate(time, .after = chlorophyll)

# Convert chlorophyll csv to raster
chloro_rast <- terra::rast(chloro, type="xyz", crs="4326")
# Error: [raster,matrix(xyz)] x cell sizes are not regular



