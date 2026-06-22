# ==============================================================================
# Script: 03_covariates_integration.R
# Purpose: Ingest topographical DEM data, perform coordinate cross-validation,
#          fit the multivariate lapse rate model, and build the 10-km prediction grid.
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(geodata)
library(terra)

# Ensure the output directory exists
dir.create("geostat-pipeline/data/processed", recursive = TRUE, showWarnings = FALSE)


# 1. Load Processed Phase 1 Datasets -------------------------------------------
message("Step 1: Loading clean datasets...")
stations_cleaned <- read_csv("geostat-pipeline/data/processed/stations_cleaned.csv", show_col_types = FALSE)
weather_cleaned <- read_csv("geostat-pipeline/data/processed/weather_cleaned.csv", show_col_types = FALSE)

# Convert stations to spatial sf points
stations_sf <- st_as_sf(stations_cleaned, coords = c("Longitude", "Latitude"), crs = 4326)


# 2. Acquire and Mask Morocco Terrain Raster (DEM) ----------------------------
message("\nStep 2: Ingesting high-resolution satellite DEM (~1-km)...")
morocco_spat <- gadm(country = "MAR", level = 0, path = tempdir())

# Download raw terrain and clip to borders
dem_raw <- elevation_30s(country = "MAR", path = tempdir())
dem_cropped <- crop(dem_raw, morocco_spat)
dem_masked <- mask(dem_cropped, morocco_spat)


# 3. Spatial Coordinate Cross-Validation ---------------------------------------
message("\nStep 3: Running topographical coordinate validation...")
extracted_vals <- terra::extract(dem_masked, stations_sf)

stations_validation <- stations_cleaned %>%
    mutate(DEM_Altitude = round(extracted_vals[[2]])) %>%
    mutate(Difference_m = abs(Altitude - DEM_Altitude)) %>%
    select(ID, Nom, Altitude, DEM_Altitude, Difference_m)

message("=== Validation Results ===")
print(as.data.frame(stations_validation))


# 4. Fit Multivariate Climatological Model -------------------------------------
message("\nStep 4: Fitting multivariate lapse rate model...")

all_time_means <- weather_cleaned %>%
    group_by(station_id) %>%
    summarise(mean_temp_all_time = mean(temperature), .groups = "drop") %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Regression: Temp as a function of Altitude and Latitude
lapse_rate_multiple <- lm(mean_temp_all_time ~ Altitude + Latitude, data = all_time_means)

message("=== Model Summary ===")
print(summary(lapse_rate_multiple))


# 5. Formulate and Export the Downsampled Prediction Grid ----------------------
message("\nStep 5: Formulating downsampled 10-km prediction grid...")

# Aggregate DEM by a factor of 10 to reduce size
dem_downsampled <- terra::aggregate(dem_masked, fact = 10, fun = "mean")

# Convert raster cells to regular spatial sf points
prediction_grid_sf <- terra::as.points(dem_downsampled) %>%
    st_as_sf() %>%
    rename(Altitude = MAR_elv_msk)

# Export the prediction grid as a serialized spatial RDS file
processed_grid_path <- "geostat-pipeline/data/processed/prediction_grid.rds"
saveRDS(prediction_grid_sf, processed_grid_path)

message("Spatial prediction grid (", nrow(prediction_grid_sf), " points) saved to: ", processed_grid_path)
message("================================================================")
