## ------------------------------------------ ##
#  Pelagic Community Structure -- Data Cleaning (Staged)
## ------------------------------------------ ##
# Script author(s): Angel Chen

# Purpose: Wrangling the NES NOAA staged data into a tidier format

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
  indices <- which(tidy_v0$CRUCOD=="(Cruise Code)")
  more_indices <- which(tidy_v0$CRUCOD=="\"ORGANISMS  >  2.5 cm\"")
  
  # Create empty list to store our tidied sections for later
  staged_list <- list()
  
  # For each section...
  for (i in 1:length(indices)){
    
    # Grab the rows for one section
    tidy_v1 <- tidy_v0 %>% dplyr::slice(indices[i]:(more_indices[i]-1))
    
    # Get the station and cruise metadata by grabbing relevant rows and columns
    station <- tidy_v1 %>%
      dplyr::slice(2) %>%
      dplyr::select(CRUNAM, STA, GERCOD, BONNUM, VOLSML, ALQFCTR, TOTCNT)
    
    # Grab the indices of the copepoda rows
    cope_indices <- which(tidy_v1$CRUCOD=="\"MAJOR COPEPODA TAXA\"")
    more_cope_indices <- which(tidy_v1$CRUCOD=="0104")
    
    # Get the copepoda data by grabbing relevant rows and columns
    # Rename columns as needed
    cope <- tidy_v1 %>% 
      dplyr::slice((cope_indices+1):more_cope_indices) %>%
      dplyr::rename(PLKTAX_NUM = CRUCOD,
             PLKTAX_NAM = ...2,
             ZOOSTG_Vial_No = ...6,
             ZOOSTG_000 = BONNUM,
             ZOOSTG_024 = 8,
             ZOOSTG_023 = ...9,
             ZOOSTG_022 = VOLSML,
             ZOOSTG_021 = ...11,
             ZOOSTG_020 = ALQFCTR,
             ZOOSTG_999 = ...13,
             ZOOCNT = TOTCNT) %>%
      dplyr::filter(PLKTAX_NUM != "PLKTAX" & PLKTAX_NUM != "NUM") %>%
      dplyr::select(-CRUNAM, -STA, -GERCOD)
      
    euph_indices <- which(tidy_v1$CRUCOD=="\"EUPHAUSIACEA TAXA\"")
    more_euph_indices <- which(tidy_v1$CRUCOD=="2003.0"|tidy_v1$CRUCOD=="2003")
    
    # Get the euphausiacea data by grabbing relevant rows and columns
    # Rename columns as needed
    euph <- tidy_v1 %>% 
      dplyr::slice((euph_indices+1):more_euph_indices) %>%
      dplyr::rename(PLKTAX_NUM = CRUCOD,
             PLKTAX_NAM = ...2,
             ZOOSTG_Vial_No = ...6,
             ZOOSTG_000 = BONNUM,
             ZOOSTG_030 = 8,
             ZOOSTG_029 = ...9,
             ZOOSTG_028 = VOLSML,
             ZOOSTG_013 = ...11,
             ZOOSTG_999 = ...13,
             ZOOCNT = TOTCNT) %>%
      dplyr::filter(PLKTAX_NUM != "PLKTAX" & PLKTAX_NUM != "NUM") %>%
      dplyr::select(-CRUNAM, -STA, -GERCOD, -ALQFCTR)
    
    fish_indices <- which(tidy_v1$CRUCOD=="\"FISH TAXA\"")
    more_fish_indices <- which(tidy_v1$CRUCOD=="3500.0"|tidy_v1$CRUCOD=="3500")
    
    # Get the fish data by grabbing relevant rows and columns
    # Rename columns as needed
    fish <- tidy_v1 %>% 
      dplyr::slice((fish_indices+1):more_fish_indices) %>%
      dplyr::rename(PLKTAX_NUM = CRUCOD,
             PLKTAX_NAM = ...2,
             ZOOSTG_Vial_No = ...6,
             ZOOSTG_000 = BONNUM,
             ZOOSTG_054 = 8,
             ZOOSTG_051 = ...9,
             ZOOSTG_050 = VOLSML,
             ZOOSTG_999 = ...13,
             ZOOCNT = TOTCNT) %>%
      dplyr::filter(PLKTAX_NUM != "PLKTAX" & PLKTAX_NUM != "NUM") %>%
      dplyr::select(-CRUNAM, -STA, -GERCOD, -ALQFCTR, -...11)
    
    # Combine the station, copepoda, euphausiacea, and fish data together
    tidy_v2 <- bind_rows(cope, euph, fish) %>%
      dplyr::bind_cols(station) %>%
      dplyr::relocate(CRUNAM:TOTCNT, .before = PLKTAX_NUM)
    
    # Save our combined result into our list
    staged_list[[i]] <- tidy_v2
  }
  
  # Unlist to combine all the sections in one file together
  staged <- staged_list %>%
    purrr::list_rbind(x = .)

  ## -------------------------------------------- ##
  #                   Export ----
  ## -------------------------------------------- ##
  
  # Make a tidy name for our exported data
  tidy_filename_0 <- gsub(pattern = ".xlsx?", replacement = "", raw_files[j])
  tidy_filename <- gsub(pattern = "-", replacement = "_", paste0(tidy_filename_0, "_staged.csv"))
  
  # Create necessary sub-folder(s)
  dir.create(path = file.path("tidy"), showWarnings = F)
  
  # Export locally
  write.csv(x = staged, file = file.path("tidy", tidy_filename), na = '', row.names = F)
  
  # Export clean dataset to Drive
  googledrive::drive_upload(media = file.path("tidy", tidy_filename), overwrite = T,
                            path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1OgDPxpq0zhFFXpCUce9j5Ez9P319k88M"))

}
