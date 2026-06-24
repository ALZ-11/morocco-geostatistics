# ==============================================================================
# Script: 03_lapse_rate_mlr.R
# Purpose: Ingest topographical DEM, perform geographical validation,
#          decouple the true physical lapse rate, and build the 10-km grid.
# ==============================================================================

library(tidyverse)
library(sf)
library(geodata)
library(terra)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Ingest and preprocess Morocco terrain (DEM)
# ------------------------------------------------------------------------------
message("downloading and masking GADM borders and high-res DEM...")

# Re-acquire GADM administrative borders
morocco_spat <- gadm(country = "MAR", level = 0, path = tempdir())

# Download the ~1 km (30 arc-second) resolution digital elevation dataset
dem_raw <- elevation_30s(country = "MAR", path = tempdir())

# Crop and mask the raw terrain model to Morocco's borders
dem_cropped <- crop(dem_raw, morocco_spat)
dem_masked <- mask(dem_cropped, morocco_spat)

# ------------------------------------------------------------------------------
# Spatial quality control (elevation cross-validation)
# ------------------------------------------------------------------------------
message("\nconducting spatial elevation cross-validation...")

stations_cleaned <- read_csv("data/processed/stations_cleaned.csv", show_col_types = FALSE)

# Convert to spatial points for extraction
stations_sf <- st_as_sf(stations_cleaned, coords = c("Longitude", "Latitude"), crs = 4326)

# Extract DEM elevation values at station coordinate locations
# terra::extract returns a dataframe; column 2 contains the elevation values.
extracted_vals <- terra::extract(dem_masked, stations_sf)

# Assemble and calculate differences
stations_validation <- stations_cleaned %>%
    mutate(
        DEM_Altitude = round(extracted_vals[[2]]),
        Difference_m = abs(Altitude - DEM_Altitude)
    ) %>%
    select(ID, Nom, Altitude, DEM_Altitude, Difference_m)

message("=== station coordinate validation results ===")
print(as.data.frame(stations_validation))

# ------------------------------------------------------------------------------
# Decouple the physical environmental lapse rate
# ------------------------------------------------------------------------------
message("\nfitting the global climatological Multiple Linear Regression (MLR)...")

weather_cleaned <- read_csv("data/processed/weather_cleaned.csv", show_col_types = FALSE)

# Calculate all-time mean temperature for each station across all years
all_time_means <- weather_cleaned %>%
    group_by(station_id) %>%
    summarise(mean_temp_all_time = mean(temperature), .groups = "drop") %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# MLR: controlling altitude against latitude confounding
lapse_rate_model <- lm(mean_temp_all_time ~ Altitude + Latitude, data = all_time_means)

message("=== MLR model summary ===")
print(summary(lapse_rate_model))

# Save the stable global model to data/processed for interpolation
saveRDS(lapse_rate_model, "data/processed/lapse_rate_model.rds")
message("Global climatological model saved to: data/processed/lapse_rate_model.rds")

# ------------------------------------------------------------------------------
# Formulate downsampled 10-km prediction grid
# ------------------------------------------------------------------------------
message("\nformulating downsampled 10-km spatial prediction grid...")

# Aggregate the 1-km DEM by a factor of 10 to create a regular ~10 km grid
# Using average cell aggregation.
dem_downsampled <- terra::aggregate(dem_masked, fact = 10, fun = "mean")

# rename raster layer before converting to points
# (to guarantee that the output sf column is named "Altitude")
names(dem_downsampled) <- "Altitude"

# Convert raster cells to a regular spatial sf point vector collection
prediction_grid_sf <- terra::as.points(dem_downsampled) %>%
    st_as_sf()

# Export
processed_grid_path <- "data/processed/prediction_grid.rds"
saveRDS(prediction_grid_sf, processed_grid_path)

message("Prediction grid (", nrow(prediction_grid_sf), " points) saved to: ", processed_grid_path)
message("================================================================")
message("Execution successful.")
