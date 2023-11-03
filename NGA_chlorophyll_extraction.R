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
  # Add an ID column
  dplyr::mutate(id = 1:nrow(.), .before = dplyr::everything()) %>% 
  # Move the time column
  dplyr::relocate(time, .after = chlorophyll) %>%
  # Remove missing chlorophyll values
  ## Makes a smaller/more manageable data object going forward
  dplyr::filter(!is.na(chlorophyll))

# Make a spatial variant
chloro_sf <- chloro %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), 
               crs = sf::st_crs(GAK_shapefile))

# Check whether that worked
plot(GAK_shapefile, axes = T, reset = F)
plot(chloro_sf, add = T, col = "black", pch = 20)

# Identify whether each point intersects the shapefile
chloro_extract <- chloro_sf %>%
  # Get a "1" if the point is in the shapefile or "NA" if not
  dplyr::mutate(ixn = as.integer(sf::st_intersects(x = geometry, y = GAK_shapefile))) %>%
  # Filter to only points that do interesect the shapefile
  dplyr::filter(!is.na(ixn))

# Check structure of that
dplyr::glimpse(chloro_extract)

# Should have far fewer rows than we started with
nrow(chloro); nrow(chloro_extract)

# Make another test graph
plot(GAK_shapefile, axes = T, reset = F)
plot(chloro_extract, add = T, col = "black", pch = 20)

# Do some final object processing
chloro_actual <- chloro_extract %>%
  # Drop sf geometry
  sf::st_drop_geometry(x = .) %>% 
  # Re-attach coordinates to the data
  dplyr::left_join(y = chloro, by = c("id", "chlorophyll", "time")) %>% 
  # Re-order columns (implicitly drop ixn and ID column)
  dplyr::select(latitude, longitude, chlorophyll, time)

# Final structure check
dplyr::glimpse(chloro_actual)
