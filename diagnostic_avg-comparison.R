## ---------------------------------------------------- ##
          # Transformation Average Comparisons
## ---------------------------------------------------- ##
# Script authors: Nick J Lyon, 

# PURPOSE
## Assess effects of various statistical transformations on the data

## ----------------------------- ##
        # Housekeeping ----
## ----------------------------- ##
# Load libraries
librarian::shelf(tidyverse)

# Make decimals not in scientific notation
options(scipen = 999)

# Clear environment
rm(list = ls())

## ----------------------------- ##
      # Data Wrangling ----
## ----------------------------- ##

# Read in data
cce_v1 <- readxl::read_excel(path = file.path("raw_data", "cce_sample-size-spectra_fixed.xlsx"))

# Check structure
dplyr::glimpse(cce_v1)

# Grab size metadata info
size_meta <- read.csv(file = file.path("raw_data", "cce_size-metadata.csv"))
dplyr::glimpse(size_meta)

# Perform needed wrangling
cce_v2 <- cce_v1 %>% 
  # Pivot to long format
  tidyr::pivot_longer(cols = dplyr::starts_with("size_"),
                      names_to = "size_category",
                      values_to = "biomass_mgCm2") %>% 
  # Tidy old column name column
  dplyr::mutate(size_category = gsub(pattern = "size_|_biomass_mgCm2", replacement = "", 
                                     x = size_category)) %>% 
  # Join on size metadata info
  dplyr::left_join(y = size_meta, by = c("size_category"))

# Re-check structure
dplyr::glimpse(cce_v2)

## ----------------------------- ##
# Perform Transformations ----
## ----------------------------- ##

# Duplicate the data object to work on it safely
focal_df <- cce_v2

# Create a list for storing graphs
graph_list <- list()

# Loop across the desired transformations
for(focal_trans in c("log", "nudge log")){
  
  # Calculate desired transformation
  ## Simple log10
  if(focal_trans == "log"){
    focal_df <- focal_df %>% 
      dplyr::mutate(trans_biomass = log10(biomass_mgCm2))
  }
  
  ## "Nudge" and take log
  if(focal_trans == "nudge log"){
    focal_df <- focal_df %>% 
      dplyr::mutate(min_biomass = min(biomass_mgCm2[biomass_mgCm2 > 0], na.rm = T),
                    trans_biomass = log10( ( biomass_mgCm2 + (min_biomass / 2) ) )) %>% 
      dplyr::select(-min_biomass)
  }
  
  # Take average of transformed biomass
  focal_avg <- focal_df %>% 
    dplyr::group_by(size_midpoint) %>% 
    dplyr::summarize(replicates = dplyr::n(),
                     mean = mean(trans_biomass, na.rm = T),
                     std_dev = sd(trans_biomass, na.rm = T),
                     std_err = std_dev / sqrt(replicates)) %>% 
    dplyr::ungroup() %>% 
    # And attach transformation type as a column
    dplyr::mutate(trans = stringr::str_to_title(focal_trans),
                  .before = dplyr::everything()) %>% 
    # Remove bad numeric values (NA or infinity)
    dplyr::filter(!is.na(mean) & mean < Inf & mean > -Inf)
  
  # Fit model and get summary
  mod_sumry <- summary(lm(mean ~ size_midpoint, data = focal_avg))
  
  # Extract coefficient table
  mod_table <- as.data.frame(mod_sumry$coefficients) %>%
    dplyr::mutate(terms = rownames(x = .), .before = dplyr::everything())
  
  # Extract slope
  mod_slope <- mod_table %>% 
    dplyr::filter(terms == "size_midpoint") %>% 
    dplyr::pull(Estimate) %>% 
    round(x = ., digits = 6)
  
  # Attach to data
  focal_avg <- focal_avg %>% 
    dplyr::mutate(slope = mod_slope)
  
  # Create graph
  p <- ggplot(focal_avg, aes(x = size_midpoint, y = mean, fill = trans)) +
    geom_text(label = paste0("Slope = ", unique(focal_avg$slope)), 
              x = 1750, y = 3.5) +
    geom_smooth(method = "lm", formula = "y ~ x", se = F, color = "black") +
    geom_errorbar(aes(ymax = mean + std_err, ymin = mean - std_err)) +
    geom_point(pch = 21, size = 2.5) +
    facet_grid(. ~ trans) +
    lims(y = c(0, 4)) +
    labs(x = "Size Midpoint", y = "Mean Biomass (mg C / m2)",
         fill = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  # Add to list
  graph_list[[focal_trans]] <- p
  
  # Success message
  message("Graph created for transformation: ", focal_trans)
  
} # Close loop

  


## ----------------------------- ##
# Graphing ----
## ----------------------------- ##


# Demo plot
ggplot(cce_actual, aes(x = size_midpoint, y = mean, fill = trans)) +
  geom_text(label = unique(cce_actual$slope), x = 2000, y = 2.5) +
  geom_smooth(method = "lm", formula = "y ~ x", se = F, color = "black") +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se)) +
  geom_point(pch = 21, size = 2.5) +
  facet_grid(. ~ trans) +
  labs(x = "Size Midpoint", y = "Mean Biomass (mg C / m2)",
       fill = "") +
  theme_bw() +
  theme(legend.position = "none")



# End ----
