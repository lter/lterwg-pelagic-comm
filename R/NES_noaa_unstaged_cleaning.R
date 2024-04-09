## ------------------------------------------ ##
#  Pelagic Community Structure -- Data Cleaning (Unstaged)
## ------------------------------------------ ##
# Script author(s): Angel Chen

# Purpose: Wrangling the NES NOAA unstaged data into a tidier format

## ------------------------------------------ ##
#            Housekeeping -----
## ------------------------------------------ ##

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

# Set site
site <- "NES"

# Create necessary sub-folder(s)
dir.create(path = file.path("raw_data"), showWarnings = F)
dir.create(path = file.path("raw_data", site), showWarnings = F)

## -------------------------------------------- ##
#             Data Acquisition ----
## -------------------------------------------- ##

# Identify raw data files
raw_NES_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/1TsZuPSjzz2mwyekoDCZgXpbKJkkzGRpv")) %>%
  dplyr::filter(name != "~$ZIW-EN627-6B3-LTER.xlsx")

# For each raw data file, download it into its own site folder
for(k in 1:nrow(raw_NES_ids)){
  
  # Download file (but silence how chatty this function is)
  googledrive::with_drive_quiet(
    googledrive::drive_download(file = raw_NES_ids[k, ]$id, overwrite = T,
                                path = file.path("raw_data", site, raw_NES_ids[k, ]$name)) )
  
  # Print success message
  message("Downloaded file ", k, " of ", nrow(raw_NES_ids))
}

## -------------------------------------------- ##
#             Data Cleaning ----
## -------------------------------------------- ##
# Identify all downloaded files
( raw_files <- dir(path = file.path("raw_data", site)) )

# For each file...
for (j in 1:length(raw_files)){
  
  # Message procesing start
  message("Cleaning '", raw_files[j], "' (file ",  j, " of ", length(raw_files), ")")
  
  # Read in file
  tidy_v0 <- readxl::read_excel(path = file.path("raw_data", site, raw_files[j]), skip = 5) %>%
    # Select only relevant columns
    dplyr::select(CRUCOD:TOTCNT)
  
  # Grab the indices of one "section"
  indices <- which(tidy_v0$CRUCOD=="\"NON-MAJOR COPEPODA TAXA\"")
  more_indices <- which(tidy_v0$CRUCOD=="4019")
  
  # Create empty list to store our tidied sections for later
  unstaged_list <- list()
  
  # For each section...
  for (i in 1:length(indices)){
    
    # Grab the rows for one section
    tidy_v1 <- tidy_v0 %>% dplyr::slice((indices[i]-1):(more_indices[i]+3)) 
    
    # Get the station and cruise metadata by grabbing relevant rows and columns
    station <- tidy_v1 %>%
      dplyr::slice(1) %>%
      dplyr::select(CRUNAM, STA, GERCOD)
    
    # Get the non-major copepoda data by grabbing relevant rows and columns
    # Rename columns as needed
    non_major <- tidy_v1 %>% 
      dplyr::slice(5:nrow(tidy_v1)) %>%
      dplyr::rename(PLKTAX_NUM = CRUCOD,
                    PLKTAX_NAM = ...2,
                    ZOOSTG_Vial_No = ...6,
                    ZOOCNT = BONNUM) %>%
      dplyr::select(PLKTAX_NUM, PLKTAX_NAM, ZOOSTG_Vial_No, ZOOCNT) %>%
      dplyr::filter(!is.na(PLKTAX_NUM))
    
    # Get the other taxa data by grabbing relevant rows and columns
    # Rename columns as needed
    other_taxa <- tidy_v1 %>% 
      dplyr::slice(5:nrow(tidy_v1)) %>%
      dplyr::rename(PLKTAX_NUM = 8,
                    PLKTAX_NAM = ...9,
                    ZOOSTG_Vial_No = ...13,
                    ZOOCNT = TOTCNT) %>%
      dplyr::select(PLKTAX_NUM, PLKTAX_NAM, ZOOSTG_Vial_No, ZOOCNT) %>%
      dplyr::filter(!is.na(PLKTAX_NUM))
    
    # Combine the station, non-major copepoda, and other taxa data together
    tidy_v2 <- bind_rows(non_major, other_taxa) %>%
      dplyr::bind_cols(station) %>%
      dplyr::relocate(CRUNAM:GERCOD, .before = PLKTAX_NUM)
    
    # Save our combined result into our list
    unstaged_list[[i]] <- tidy_v2
  }
  
  # Unlist to combine all the sections in one file together
  unstaged <- unstaged_list %>%
    purrr::list_rbind(x = .)
  
  
  ## -------------------------------------------- ##
  #                   Export ----
  ## -------------------------------------------- ##
  
  # Make a tidy name for our exported data
  tidy_filename_0 <- gsub(pattern = ".xlsx?", replacement = "", raw_files[j])
  tidy_filename <- gsub(pattern = "-", replacement = "_", paste0(tidy_filename_0, "_unstaged.csv"))
  
  # Create necessary sub-folder(s)
  dir.create(path = file.path("tidy"), showWarnings = F)
  
  # Export locally
  write.csv(x = unstaged, file = file.path("tidy", tidy_filename), na = '', row.names = F)
  
  # Export clean dataset to Drive
  googledrive::drive_upload(media = file.path("tidy", tidy_filename), overwrite = T,
                            path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1VRnPKGxWULZ4NyGZfgzTMh0H_Cr6WJsn"))
  
}
