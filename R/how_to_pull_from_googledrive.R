## ------------------------------------------ ##
#        Pelagic Community Structure 
# Example Data Acquisition from Google Drive
## ------------------------------------------ ##
# Script author(s): Angel Chen

# Purpose:
## An example script showing how you can pull in data from Google Drive to use in your local RStudio session

## ------------------------------------------ ##
#            Housekeeping -----
## ------------------------------------------ ##

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

# Set site
site <- "NGA"

# Create necessary sub-folder(s)
dir.create(path = file.path("raw_data"), showWarnings = F)
dir.create(path = file.path("raw_data", site), showWarnings = F)

## -------------------------------------------- ##
#             Data Acquisition ----
## -------------------------------------------- ##

# Identify raw data files
# For example, here I'm pulling all the NGA csv files from Google Drive
# A new window will pop up asking you to select the appropriate Google Drive account
# For more help, see: https://nceas.github.io/scicomp.github.io/tutorials.html#using-the-googledrive-r-package
raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M"),
                                     type = "csv") 

# Identify raw data files
# For example, here I'm pulling specific NGA files from Google Drive
# raw_NGA_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/14lYLzawcQUy0ruAG6norhJ4VcocfPk8M")) %>%
#   dplyr::filter(name %in% c("NGA_signature_microzoop_carbon_biomass_2011_2022.csv ",
#                             "NGA_signature_microzoop_abundance_2011_2022.csv"))

# For each raw data file, download it into its own site folder
for(k in 1:nrow(raw_NGA_ids)){
  
  # Download file (but silence how chatty this function is)
  googledrive::with_drive_quiet(
    googledrive::drive_download(file = raw_NGA_ids[k, ]$id, overwrite = T,
                                path = file.path("raw_data", site, raw_NGA_ids[k, ]$name)) )
  
  # Print success message
  message("Downloaded file ", k, " of ", nrow(raw_NGA_ids))
}

# Read in a csv
NGA_microzoop_abundance <- read.csv(file = file.path("raw_data",
                                                     site,
                                                     "NGA_signature_microzoop_abundance_2011_2022.csv"))
