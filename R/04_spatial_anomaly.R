# ==============================================================================
# Script: 04_spatial_anomaly.R
# Purpose: Extract daily weather slice, compute spatial residual anomalies,
#          define the prescribed geostatistical variogram, and export diagnostics.
# ==============================================================================

library(tidyverse)
library(sf)
library(gstat)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Slice daily weather and load spatial trend model
# ------------------------------------------------------------------------------
message("ingesting datasets and isolating target date (June 30, 1961)...")

stations_cleaned <- read_csv("data/processed/stations_cleaned.csv", show_col_types = FALSE)
weather_cleaned <- read_csv("data/processed/weather_cleaned.csv", show_col_types = FALSE)

lapse_rate_model <- readRDS("data/processed/lapse_rate_model.rds")

target_date <- as.Date("1961-06-30")

# Extract weather records strictly for our target date
oneday_weather <- weather_cleaned %>%
    filter(date == target_date) %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Convert daily data to a spatial sf point collection (WGS 84 / EPSG:4326)
# remove = FALSE used to preserve 'Latitude' and 'Longitude' as numeric columns
# which ensures predict() can resolve variables correctly without geometry parsing.
oneday_sf <- st_as_sf(
    oneday_weather,
    coords = c("Longitude", "Latitude"),
    crs = 4326,
    remove = FALSE
)

# ------------------------------------------------------------------------------
# Calculate daily spatial residual anomalies
# ------------------------------------------------------------------------------
message("conducting daily spatial detrending...")

# Use the global MLR trend to calculate expected climatological temperatures
oneday_sf <- oneday_sf %>%
    mutate(
        expected_temp = predict(lapse_rate_model, newdata = .),
        # Residual represents the daily anomaly unaccounted for by geography
        residual = temperature - expected_temp
    )

message("=== daily anomaly (residual) calculations ===")
print(oneday_sf %>% select(Nom, temperature, expected_temp, residual) %>% as.data.frame() %>% select(-geometry))

# Save the target day spatial residuals
saveRDS(oneday_sf, "data/processed/oneday_residuals_sf.rds")
message("Daily spatial residuals saved to: data/processed/oneday_residuals_sf.rds")

# ------------------------------------------------------------------------------
# Define the prescribed spatial variogram model
# ------------------------------------------------------------------------------
message("\nformulating prescribed climatological variogram model...")

# Compute the empirical variogram of the residuals
# Since the trend has been subtracted, residuals are stationary (residuals ~ 1)
emp_var <- variogram(residual ~ 1, oneday_sf)

# Define a prescribed Exponential spatial covariance model
# Because the input sf uses geographic coordinates (degrees), gstat computes
# spatial distance in kilometers. so range = 150 represents 150 km.
# Setting a realistic correlation range prevents spatial model collapsing.
prescribed_vgm <- vgm(
    psill = 5.0,
    model = "Exp",
    range = 150,
    nugget = 1.5
)

# Save the prescribed variogram model object for the Kriging step
saveRDS(prescribed_vgm, "data/processed/prescribed_variogram.rds")
message("Prescribed variogram model saved to: data/processed/prescribed_variogram.rds")

# ------------------------------------------------------------------------------
# Export Diagnostic Variogram Plot
# ------------------------------------------------------------------------------
message("exporting diagnostic variogram plot...")

png_path <- "output/figures/04_variogram_fit.png"
png(filename = png_path, width = 800, height = 500, res = 120)

v_plot <- plot(
    emp_var,
    model = prescribed_vgm,
    main = "Prescribed Residual Variogram Model (Range = 150 km)",
    xlab = "Distance (km)",
    ylab = "Semivariance"
)

print(v_plot)

dev.off()

message("Diagnostic variogram plot saved to: ", png_path)
message("================================================================")
message("Execution successful.")
