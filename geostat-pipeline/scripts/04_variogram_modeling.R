# ==============================================================================
# Script: 04_variogram_modeling.R
# Purpose: Extract spatial residuals from the trend model, compute empirical
#          residual variography, fit a stable locked-range spatial model,
#          and export diagnostic figures and seralized variogram objects.
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(gstat)

# Ensure output directories exist
dir.create("geostat-pipeline/data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("geostat-pipeline/output/figures", recursive = TRUE, showWarnings = FALSE)


# 1. Load Processed Datasets ---------------------------------------------------
message("Step 1: Loading clean datasets...")
stations_cleaned <- read_csv("geostat-pipeline/data/processed/stations_cleaned.csv", show_col_types = FALSE)
weather_cleaned <- read_csv("geostat-pipeline/data/processed/weather_cleaned.csv", show_col_types = FALSE)


# 2. Extract Spatial Residuals from Golden Year (1961) ------------------------
message("Step 2: Isolating 1961 spatial residuals...")

# Recompile the 1961 baseline (8 active stations)
station_means_1961 <- weather_cleaned %>%
    filter(date >= as.Date("1961-01-01"), date <= as.Date("1961-12-31")) %>%
    group_by(station_id) %>%
    summarise(mean_temp_1961 = mean(temperature), .groups = "drop") %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Fit the first-order trend model (Temperature ~ Altitude + Latitude)
trend_model <- lm(mean_temp_1961 ~ Altitude + Latitude, data = station_means_1961)

# Extract and append the residuals
station_means_1961$residuals <- resid(trend_model)

# Convert to a spatial sf vector object (WGS 84 / EPSG:4326)
stations_res_sf <- st_as_sf(
    station_means_1961,
    coords = c("Longitude", "Latitude"),
    crs = 4326
)


# 3. Compute and Fit Variogram -------------------------------------------------
message("Step 3: Computing and fitting spatial variogram...")

# Compute empirical variogram of residuals
emp_var <- variogram(residuals ~ 1, stations_res_sf)

# Define initial guess and fit model while locking the range at 100 km
init_model <- vgm(psill = 3, model = "Exp", range = 100, nugget = 1.5)
fitted_locked_range <- fit.variogram(
    emp_var,
    model = init_model,
    fit.method = 1,
    fit.ranges = FALSE
)

# Export the fitted variogram model as an RDS file
processed_var_path <- "geostat-pipeline/data/processed/fitted_variogram.rds"
saveRDS(fitted_locked_range, processed_var_path)
message("Fitted variogram model saved to: ", processed_var_path)


# 4. Generate and Export Diagnostic Plot ---------------------------------------
message("Step 4: Exporting diagnostic plot...")

# Save the diagnostic plot using R's native png device
png_path <- "geostat-pipeline/output/figures/04_variogram_fit.png"
png(filename = png_path, width = 800, height = 500, res = 120)

# Capture the lattice/trellis plot object
p <- plot(
    emp_var,
    model = fitted_locked_range,
    main = "Residual Spatial Variogram Fit (Locked Range = 100 km)"
)

# Explicitly print the trellis object to the active png device
print(p)

dev.off()

message("Diagnostic plot saved to: ", png_path)
message("================================================================")
