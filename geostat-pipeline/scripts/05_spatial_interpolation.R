# ==============================================================================
# Script: 05_spatial_interpolation.R
# Purpose: Execute Ordinary Kriging, run the stable climatological model,
#          and export high-resolution continuous temperature maps of Morocco.
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

prediction_grid_sf <- readRDS("geostat-pipeline/data/processed/prediction_grid.rds")
fitted_variogram <- readRDS("geostat-pipeline/data/processed/fitted_variogram.rds")


# 2. Slice Target Date (June 30, 1961) ----------------------------------------
message("Step 2: Slicing weather database for June 30, 1961...")
target_date <- as.Date("1961-06-30")

oneday_data <- weather_cleaned %>%
    filter(date == target_date) %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Convert to spatial sf points (WGS 84) and extract coordinates
oneday_sf <- st_as_sf(oneday_data, coords = c("Longitude", "Latitude"), crs = 4326)
oneday_sf$Latitude <- st_coordinates(oneday_sf)[, 2]

# Extract Latitude for the prediction grid
prediction_grid_sf$Latitude <- st_coordinates(prediction_grid_sf)[, 2]


# 3. Execute Ordinary Kriging (OK) ---------------------------------------------
message("Step 3: Executing Ordinary Kriging...")
ok_results <- krige(
    formula = temperature ~ 1,
    locations = oneday_sf,
    newdata = prediction_grid_sf,
    model = fitted_variogram
)

# Export Ordinary Kriging results
processed_ok_path <- "geostat-pipeline/data/processed/ok_predictions.rds"
saveRDS(ok_results, processed_ok_path)


# 4. Execute Stable Climatological Model ---------------------------------------
message("Step 4: Executing stable climatological model on grid...")
prediction_grid_sf <- prediction_grid_sf %>%
    mutate(climatology_pred = 67.37 - 0.00213 * Altitude - 1.271 * Latitude)

# Export Climatological predictions
processed_clim_path <- "geostat-pipeline/data/processed/climatology_predictions.rds"
saveRDS(prediction_grid_sf, processed_clim_path)


# 5. Download Morocco Boundary for Plotting Map Base --------------------------
morocco_spat <- geodata::gadm(country = "MAR", level = 0, path = tempdir())
morocco_sf <- st_as_sf(morocco_spat)


# 6. Generate and Export Comparative Figures ----------------------------------
message("Step 5: Exporting final prediction maps...")

# 6.1. Save the flat Ordinary Kriging Map
map_ok <- ggplot() +
    geom_sf(data = morocco_sf, fill = "#f5f5f5", color = "#8c8c8c", linewidth = 0.4) +
    geom_sf(data = ok_results, aes(color = var1.pred), size = 1.5) +
    scale_color_viridis_c(option = "plasma", name = "Temp (°C)", limits = c(10, 35)) +
    labs(
        title = "Ordinary Kriging Prediction Surface (June 30, 1961)",
        subtitle = "Pure Nugget Model: Predicts the global mean due to spatial uncertainty",
        x = "Longitude", y = "Latitude",
        caption = "Source: GADM & ECA&D"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8))

ggsave(
    filename = "geostat-pipeline/output/figures/05_ordinary_kriging_map.png",
    plot = map_ok, width = 8, height = 6, dpi = 300, bg = "white"
)

# 6.2. Save the Stable Climatological Map
map_stable_climate <- ggplot() +
    geom_sf(data = morocco_sf, fill = "#f5f5f5", color = "#8c8c8c", linewidth = 0.4) +
    geom_sf(data = prediction_grid_sf, aes(color = climatology_pred), size = 1.5) +
    scale_color_viridis_c(option = "plasma", name = "Temp (°C)", limits = c(10, 35)) +
    labs(
        title = "Morocco Climatological Temperature Surface",
        subtitle = "Model: Temp ~ Altitude + Latitude (Stable Longitudinal & Topographic Gradients)",
        x = "Longitude", y = "Latitude",
        caption = "Source: GADM & ECA&D"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "#555555", size = 8.5)
    )

ggsave(
    filename = "geostat-pipeline/output/figures/06_climatological_surface_map.png",
    plot = map_stable_climate, width = 8, height = 6, dpi = 300, bg = "white"
)

message("================================================================")
message("All Phase 5 predictions and maps generated successfully")
