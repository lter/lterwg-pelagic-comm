## ------------------------------------------ ##
#        Pelagic Community Structure 
#          CCE ZooScan Synthesizing 
## ------------------------------------------ ##
# Script author(s): Angel Chen

# Purpose: Synthesize ZooScan variables (ugcarbon mean, ESD mean, Feret mean, abundance)
# across F1, F2, F3 datasets

## ------------------------------------------ ##
#            Housekeeping -----
## ------------------------------------------ ##

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

# Set site
site <- "CCE"

# Create necessary sub-folder(s)
dir.create(path = file.path("raw_data"), showWarnings = F)
dir.create(path = file.path("raw_data", site), showWarnings = F)

## -------------------------------------------- ##
#             Data Acquisition ----
## -------------------------------------------- ##

# Identify raw data files
# For example, here I'm pulling all the CCE csv files from Google Drive
# A new window will pop up asking you to select the appropriate Google Drive account
# For more help, see: https://nceas.github.io/scicomp.github.io/tutorials.html#using-the-googledrive-r-package
raw_CCE_F1 <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/18wxLfixf9RB8mpNXlujClBIBDk42Ebxg"),
                                     type = "csv") 

raw_CCE_F2 <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/18tGQWMsgeQDYQUKOsThBWslpjoo_CpPv"),
                                    type = "csv") 

raw_CCE_F3 <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/198JBjjhBwu4u7ZALPm-NKmATEST4zfQM"),
                                    type = "csv") 

raw_ids <- rbind(raw_CCE_F1, raw_CCE_F2, raw_CCE_F3)

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
#            Data Synthesizing -----
## ------------------------------------------ ##

# --------------UGCARBON MEAN---------------- # 

target_variable <- "ugcarbon_mean"

sizes <- c("F1", "F2", "F3")

# Make an empty list to store re-formatted raw data
df_list <- list()

# For each size type...
for (a_size_type in 1:length(sizes)){
  # Grab the size
  size_type <- sizes[a_size_type]
  
  # Subset the raw files to that size type
  raw_files_subset <- dir(path = file.path("raw_data", site), pattern = size_type) 

  # For each raw file of that size type...  
  for (i in 1:length(raw_files_subset)){
    # Grab the name
    raw_file_name <- raw_files_subset[i]
    
    # Read in the csv
    raw_df_v1 <- read.csv(file = file.path("raw_data", site, raw_file_name))
    
    raw_df_v2 <- raw_df_v1 %>%
      # Select the columns of interest
      dplyr::select(Cruise, Line, Station, Time_PST, Class, !!sym(target_variable)) %>%
      # Transpose the species names
      tidyr::pivot_wider(names_from = Class, 
                         values_from = !!sym(target_variable)) %>%
      # Rename columns with size type indicator
      dplyr::rename_with(.fn = ~paste0(size_type, "_", .x), .cols = 5:last_col())
    
    # Add to list
    df_list[[raw_file_name]] <- raw_df_v2
  }
  
}

# Unlist the list we just generated
tidy_v0 <- df_list %>%
  purrr::list_rbind(x = .) %>%
  dplyr::arrange(Cruise, Line, Station, Time_PST)

tidy_v1 <- tidy_v0 %>%
  # Group by the 4 identifier columns
  dplyr::group_by(Cruise, Line, Station, Time_PST) %>%
  # Coalesce across grouped rows by filling in missing values...
  tidyr::fill(everything(), .direction = "downup") %>%
  # Then filter just one row for each group 
  # For more info: https://stackoverflow.com/questions/45515218/combine-rows-in-data-frame-containing-na-to-make-complete-row
  dplyr::slice(1)

# Check that out
dplyr::glimpse(tidy_v1)

# Create the export file name
tidy_filename <- paste0("CCE_synthesized_", target_variable, ".csv")

# Create necessary sub-folder(s)
dir.create(path = file.path("tidy"), showWarnings = F)

# Export locally
write.csv(x = tidy_v1, file = file.path("tidy", tidy_filename), na = "NaN", row.names = F)

# Export to Drive
googledrive::drive_upload(media = file.path("tidy", tidy_filename), overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1ipqKV1kwIMXF5CwydyGGFDkY2noe_wSf"))

# Clear environment
rm(list = setdiff(ls(), "site"))

# --------------ESD MEAN---------------- # 

target_variable <- "ESD_mean"

sizes <- c("F1", "F2", "F3")

# Make an empty list to store re-formatted raw data
df_list <- list()

# For each size type...
for (a_size_type in 1:length(sizes)){
  # Grab the size
  size_type <- sizes[a_size_type]
  
  # Subset the raw files to that size type
  raw_files_subset <- dir(path = file.path("raw_data", site), pattern = size_type) 
  
  # For each raw file of that size type...  
  for (i in 1:length(raw_files_subset)){
    # Grab the name
    raw_file_name <- raw_files_subset[i]
    
    # Read in the csv
    raw_df_v1 <- read.csv(file = file.path("raw_data", site, raw_file_name))
    
    raw_df_v2 <- raw_df_v1 %>%
      # Select the columns of interest
      dplyr::select(Cruise, Line, Station, Time_PST, Class, !!sym(target_variable)) %>%
      # Transpose the species names
      tidyr::pivot_wider(names_from = Class, 
                         values_from = !!sym(target_variable)) %>%
      # Rename columns with size type indicator
      dplyr::rename_with(.fn = ~paste0(size_type, "_", .x), .cols = 5:last_col())
    
    # Add to list
    df_list[[raw_file_name]] <- raw_df_v2
  }
  
}

# Unlist the list we just generated
tidy_v0 <- df_list %>%
  purrr::list_rbind(x = .) %>%
  dplyr::arrange(Cruise, Line, Station, Time_PST)

tidy_v1 <- tidy_v0 %>%
  # Group by the 4 identifier columns
  dplyr::group_by(Cruise, Line, Station, Time_PST) %>%
  # Coalesce across grouped rows by filling in missing values...
  tidyr::fill(everything(), .direction = "downup") %>%
  # Then filter just one row for each group 
  # For more info: https://stackoverflow.com/questions/45515218/combine-rows-in-data-frame-containing-na-to-make-complete-row
  dplyr::slice(1)

# Check that out
dplyr::glimpse(tidy_v1)

# Create the export file name
tidy_filename <- paste0("CCE_synthesized_", target_variable, ".csv")

# Create necessary sub-folder(s)
dir.create(path = file.path("tidy"), showWarnings = F)

# Export locally
write.csv(x = tidy_v1, file = file.path("tidy", tidy_filename), na = "NaN", row.names = F)

# Export to Drive
googledrive::drive_upload(media = file.path("tidy", tidy_filename), overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1ipqKV1kwIMXF5CwydyGGFDkY2noe_wSf"))

# Clear environment
rm(list = setdiff(ls(), "site"))

# --------------FERET MEAN---------------- # 

target_variable <- "Feret_mean"

sizes <- c("F1", "F2", "F3")

# Make an empty list to store re-formatted raw data
df_list <- list()

# For each size type...
for (a_size_type in 1:length(sizes)){
  # Grab the size
  size_type <- sizes[a_size_type]
  
  # Subset the raw files to that size type
  raw_files_subset <- dir(path = file.path("raw_data", site), pattern = size_type) 
  
  # For each raw file of that size type...  
  for (i in 1:length(raw_files_subset)){
    # Grab the name
    raw_file_name <- raw_files_subset[i]
    
    # Read in the csv
    raw_df_v1 <- read.csv(file = file.path("raw_data", site, raw_file_name))
    
    raw_df_v2 <- raw_df_v1 %>%
      # Select the columns of interest
      dplyr::select(Cruise, Line, Station, Time_PST, Class, !!sym(target_variable)) %>%
      # Transpose the species names
      tidyr::pivot_wider(names_from = Class, 
                         values_from = !!sym(target_variable)) %>%
      # Rename columns with size type indicator
      dplyr::rename_with(.fn = ~paste0(size_type, "_", .x), .cols = 5:last_col())
    
    # Add to list
    df_list[[raw_file_name]] <- raw_df_v2
  }
  
}

# Unlist the list we just generated
tidy_v0 <- df_list %>%
  purrr::list_rbind(x = .) %>%
  dplyr::arrange(Cruise, Line, Station, Time_PST)

tidy_v1 <- tidy_v0 %>%
  # Group by the 4 identifier columns
  dplyr::group_by(Cruise, Line, Station, Time_PST) %>%
  # Coalesce across grouped rows by filling in missing values...
  tidyr::fill(everything(), .direction = "downup") %>%
  # Then filter just one row for each group 
  # For more info: https://stackoverflow.com/questions/45515218/combine-rows-in-data-frame-containing-na-to-make-complete-row
  dplyr::slice(1)

# Check that out
dplyr::glimpse(tidy_v1)

# Create the export file name
tidy_filename <- paste0("CCE_synthesized_", target_variable, ".csv")

# Create necessary sub-folder(s)
dir.create(path = file.path("tidy"), showWarnings = F)

# Export locally
write.csv(x = tidy_v1, file = file.path("tidy", tidy_filename), na = "NaN", row.names = F)

# Export to Drive
googledrive::drive_upload(media = file.path("tidy", tidy_filename), overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1ipqKV1kwIMXF5CwydyGGFDkY2noe_wSf"))

# Clear environment
rm(list = setdiff(ls(), "site"))

# --------------ABUNDANCE---------------- # 

target_variable <- "Abundance..no.m.2."

sizes <- c("F1", "F2", "F3")

# Make an empty list to store re-formatted raw data
df_list <- list()

# For each size type...
for (a_size_type in 1:length(sizes)){
  # Grab the size
  size_type <- sizes[a_size_type]
  
  # Subset the raw files to that size type
  raw_files_subset <- dir(path = file.path("raw_data", site), pattern = size_type) 
  
  # For each raw file of that size type...  
  for (i in 1:length(raw_files_subset)){
    # Grab the name
    raw_file_name <- raw_files_subset[i]
    
    # Read in the csv
    raw_df_v1 <- read.csv(file = file.path("raw_data", site, raw_file_name))
    
    raw_df_v2 <- raw_df_v1 %>%
      # Select the columns of interest
      dplyr::select(Cruise, Line, Station, Time_PST, Class, !!sym(target_variable)) %>%
      # Transpose the species names
      tidyr::pivot_wider(names_from = Class, 
                         values_from = !!sym(target_variable)) %>%
      # Rename columns with size type indicator
      dplyr::rename_with(.fn = ~paste0(size_type, "_", .x), .cols = 5:last_col())
    
    # Add to list
    df_list[[raw_file_name]] <- raw_df_v2
  }
  
}

# Unlist the list we just generated
tidy_v0 <- df_list %>%
  purrr::list_rbind(x = .) %>%
  dplyr::arrange(Cruise, Line, Station, Time_PST)

tidy_v1 <- tidy_v0 %>%
  # Group by the 4 identifier columns
  dplyr::group_by(Cruise, Line, Station, Time_PST) %>%
  # Coalesce across grouped rows by filling in missing values...
  tidyr::fill(everything(), .direction = "downup") %>%
  # Then filter just one row for each group 
  # For more info: https://stackoverflow.com/questions/45515218/combine-rows-in-data-frame-containing-na-to-make-complete-row
  dplyr::slice(1)

# Check that out
dplyr::glimpse(tidy_v1)

# Create the export file name
tidy_filename <- paste0("CCE_synthesized_", target_variable)
tidy_filename <- gsub(pattern = "\\.\\.|\\.", x = tidy_filename, replacement = "_")
tidy_filename <- gsub(pattern = "2_", x = tidy_filename, replacement = "2.csv")

# Create necessary sub-folder(s)
dir.create(path = file.path("tidy"), showWarnings = F)

# Export locally
write.csv(x = tidy_v1, file = file.path("tidy", tidy_filename), na = "NaN", row.names = F)

# Export to Drive
googledrive::drive_upload(media = file.path("tidy", tidy_filename), overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1ipqKV1kwIMXF5CwydyGGFDkY2noe_wSf"))

# # Read in a csv
# L80S70_F1 <- read.csv(file = file.path("raw_data", site, "L80S70_F1.csv"))
# 
# test <- L80S70_F1 %>%
#   dplyr::select(Cruise, Line, Station, Time_PST, Class, !!sym(target_variable)) %>%
#   tidyr::pivot_wider(names_from = Class, 
#                      values_from = !!sym(target_variable)) %>%
#   dplyr::rename_with(.fn = ~paste0(size_type, "_", .x), .cols = 5:last_col())
