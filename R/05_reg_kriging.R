# ==============================================================================
# Script: 05_reg_kriging.R
# Purpose: Execute Ordinary Kriging on daily anomalies, reconstruct the
#          terrain-adjusted daily temperature surface, and export maps.
# ==============================================================================

library(tidyverse)
library(sf)
library(gstat)
library(geodata)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load spatial grids, residuals, and global trend model
# ------------------------------------------------------------------------------
message("ingesting spatial grids, residuals, and global trend model...")

# Load spatial prediction grid, target day residuals, and trend model
prediction_grid_sf <- readRDS("data/processed/prediction_grid.rds")
oneday_sf <- readRDS("data/processed/oneday_residuals_sf.rds")
lapse_rate_model <- readRDS("data/processed/lapse_rate_model.rds")
prescribed_vgm <- readRDS("data/processed/prescribed_variogram.rds")

# Extract Latitude and Longitude columns from prediction grid geometries.
# Since prediction_grid_sf is compiled from raw points, it lacks explicit numeric columns
# for coordinates, which causes predict() to fail. Extracting them dynamically solves this.
grid_coords <- st_coordinates(prediction_grid_sf)
prediction_grid_sf$Longitude <- grid_coords[, 1]
prediction_grid_sf$Latitude <- grid_coords[, 2]

# ------------------------------------------------------------------------------
# Execute Ordinary Kriging (Anomaly vs. Baseline)
# ------------------------------------------------------------------------------
message("executing ordinary kriging interpolations...")

# Interpolation A: Krige raw temperatures to create our naive baseline
kriged_raw <- krige(
    formula = temperature ~ 1,
    locations = oneday_sf,
    newdata = prediction_grid_sf,
    model = prescribed_vgm
)

# Interpolation B: Krige daily residuals (anomalies) for Regression-Kriging
kriged_residuals <- krige(
    formula = residual ~ 1,
    locations = oneday_sf,
    newdata = prediction_grid_sf,
    model = prescribed_vgm
)

# ------------------------------------------------------------------------------
# Recombine trend & anomaly fields (RK Reconstruction)
# ------------------------------------------------------------------------------
message("rebuilding daily temperature surface (Trend + Kriged Anomaly)...")

# Predict the stable climatological trend across the entire 10-km prediction grid
prediction_grid_sf$expected_temp <- predict(lapse_rate_model, newdata = prediction_grid_sf)

# Add the kriged daily residual anomaly back to the trend grid to reconstruct RK
prediction_grid_sf$rk_pred <- prediction_grid_sf$expected_temp + kriged_residuals$var1.pred

# Save the final continuous predictions for later validation
saveRDS(prediction_grid_sf, "data/processed/climatology_predictions.rds")
message("Final Regression-Kriging predictions saved to: data/processed/climatology_predictions.rds")

# ------------------------------------------------------------------------------
# Map and plot spatial surfaces
# ------------------------------------------------------------------------------
message("compiling comparative maps...")

# Re-acquire GADM national boundary vector and "Western Sahara" for map backing
morocco_spat <- geodata::gadm(country = "MAR", level = 0, path = tempdir())
ws_spat <- geodata::gadm(country = "ESH", level = 0, path = tempdir())

# Convert to sf and union them to dissolve the internal border
morocco_sf <- st_union(st_as_sf(morocco_spat), st_as_sf(ws_spat)) %>%
    st_as_sf()

temp_scale <- scale_color_viridis_c(
    option = "plasma",
    name = "Temp (°C)",
    limits = c(10, 48)
)

# Map 4.1: Standard Ordinary Kriging Surface (Baseline)
map_ok <- ggplot() +
    geom_sf(data = morocco_sf, fill = "#f5f5f5", color = "#8c8c8c", linewidth = 0.4) +
    geom_sf(data = kriged_raw, aes(color = var1.pred), size = 1.2) +
    temp_scale +
    labs(
        title = "Ordinary Kriging Prediction Surface (June 30, 1961)",
        subtitle = "Naïve Distance Model: Severe spatial smoothing due to station sparsity",
        x = "Longitude", y = "Latitude",
        caption = "Source: GADM national borders & ECA&D weather station observations"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "#555555", size = 8.5),
        plot.caption = element_text(size = 7, color = "#777777")
    )

ggsave(
    filename = "output/figures/05_ordinary_kriging_map.png",
    plot = map_ok, width = 8, height = 6, dpi = 300, bg = "white"
)
message("Naive Ordinary Kriging map saved to: output/figures/05_ordinary_kriging_map.png")

# Map 4.2: Regression-Kriging Surface (Continuous & Terrain-Carved)
map_rk <- ggplot() +
    geom_sf(data = morocco_sf, fill = "#f5f5f5", color = "#8c8c8c", linewidth = 0.4) +
    geom_sf(data = prediction_grid_sf, aes(color = rk_pred), size = 1.2) +
    temp_scale +
    labs(
        title = "Morocco Regression-Kriging Daily Surface (June 30, 1961)",
        subtitle = "Reconstructed Model: Expected Climatology (MLR) + Interpolated Local Weather Anomaly",
        x = "Longitude", y = "Latitude",
        caption = "Source: GADM national borders, satellite DEM, and ECA&D weather observations"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "#555555", size = 8.5),
        plot.caption = element_text(size = 7, color = "#777777")
    )

ggsave(
    filename = "output/figures/06_climatological_surface_map.png",
    plot = map_rk, width = 8, height = 6, dpi = 300, bg = "white"
)
message("Regression-Kriging map saved to: output/figures/06_climatological_surface_map.png")
message("================================================================")
message("Execution successful.")
