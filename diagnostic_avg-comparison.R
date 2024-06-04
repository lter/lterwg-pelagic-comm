
# Load libraries
librarian::shelf(tidyverse)

# Make decimals not in scientific notation
options(scipen = 999)

# Clear environment
rm(list = ls())

# Read in data
cce_df <- readxl::read_excel(path = file.path("raw_data", "cce_sample-size-spectra_fixed.xlsx"))

# Check structure
dplyr::glimpse(cce_df)

# Grab size metadata info
size_meta <- read.csv(file = file.path("raw_data", "cce_size-metadata.csv"))
dplyr::glimpse(size_meta)

# Perform needed wrangling
cce_v2 <- cce_df %>% 
  # Pivot to long format
  tidyr::pivot_longer(cols = dplyr::starts_with("size_"),
                      names_to = "size_category",
                      values_to = "biomass_mgCm2") %>% 
  # Tidy old column name column
  dplyr::mutate(size_category = gsub(pattern = "size_|_biomass_mgCm2", replacement = "", 
                                     x = size_category)) %>% 
  # Join on size metadata info
  dplyr::left_join(y = size_meta, by = c("size_category")) %>% 
  # Log transform in diff ways
  dplyr::mutate(log_biomass = log10(biomass_mgCm2),
                min_biomass = min(biomass_mgCm2[biomass_mgCm2 > 0], na.rm = T),
                nudge_log_biomass = log10( ( biomass_mgCm2 + (min_biomass / 2) ) ))

# Re-check structure
dplyr::glimpse(cce_v2)

# Perform different averaging methods
cce_avg <- cce_v2 %>% 
  dplyr::group_by(size_midpoint) %>% 
  dplyr::summarize(replicates = dplyr::n(),
                   log_mean = mean(log_biomass, na.rm = T),
                   log_sd = sd(log_biomass, na.rm = T),
                   log_se = log_sd / sqrt(replicates),
                   nudge_mean = mean(nudge_log_biomass, na.rm = T),
                   nudge_sd = sd(nudge_log_biomass, na.rm = T),
                   nudge_se = nudge_sd / sqrt(replicates)) %>% 
  dplyr::ungroup() %>% 
  # Pivot longer
  tidyr::pivot_longer(cols = dplyr::starts_with(c("log_", "nudge_"))) %>% 
  tidyr::separate_wider_delim(cols = name, delim = "_", names = c("trans", "metric")) %>% 
  # Pivot wider
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  # Tweak transformation names
  dplyr::mutate(trans = dplyr::case_when(
    trans == "log" ~ "Log",
    trans == "nudge" ~ "Log Nudge",
    T ~ trans))
  
# Structure check
dplyr::glimpse(cce_avg)
# view(cce_avg)

# Duplicate object to work on it safely
cce_actual <- cce_avg %>% 
  dplyr::mutate(slope = NA)

# Identify slopes
for(focal_trans in unique(cce_avg$trans)){
  
  # Subset the data to only that transformation
  data_sub <- cce_avg %>% 
    dplyr::filter(trans == focal_trans) %>% 
    dplyr::filter(!is.na(mean) & mean < Inf & mean > -Inf)
  
  # Fit model and get summary
  mod_sumry <- summary(lm(mean ~ size_midpoint, data = data_sub))
  
  # Extract coefficient table
  mod_table <- as.data.frame(mod_sumry$coefficients) %>%
    dplyr::mutate(terms = rownames(x = .), .before = dplyr::everything())
  
  # Extract slope
  mod_slope <- mod_table %>% 
    dplyr::filter(terms == "size_midpoint") %>% 
    dplyr::pull(Estimate) %>% 
    round(x = ., digits = 6)
  
  # Attach to data
  cce_actual <- cce_actual %>% 
    dplyr::mutate(slope = ifelse(trans == focal_trans,
                                 yes = mod_slope,
                                 no = slope))
}

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




