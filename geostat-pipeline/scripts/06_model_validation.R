# ==============================================================================
# Script: 06_model_validation.R
# Purpose: Execute Leave-One-Out Cross-Validation (LOOCV) across Ordinary and
#          Universal Kriging models, compute standardized error metrics,
#          export benchmarking datasets, and save residual figures.
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
fitted_variogram <- readRDS("geostat-pipeline/data/processed/fitted_variogram.rds")


# 2. Slice Target Date (June 30, 1961) ----------------------------------------
message("Step 2: Preparing single-day spatial slice...")
target_date <- as.Date("1961-06-30")

oneday_data <- weather_cleaned %>%
    filter(date == target_date) %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Convert to spatial sf points (WGS 84) and extract coordinates
oneday_sf <- st_as_sf(oneday_data, coords = c("Longitude", "Latitude"), crs = 4326)
oneday_sf$Latitude <- st_coordinates(oneday_sf)[, 2]


# 3. Execute Cross-Validation (LOOCV) ------------------------------------------
message("Step 3: Executing LOOCV loop across models...")

cv_ok <- krige.cv(temperature ~ 1, locations = oneday_sf, model = fitted_variogram, verbose = FALSE)
cv_uk1 <- krige.cv(temperature ~ Altitude, locations = oneday_sf, model = fitted_variogram, verbose = FALSE)
cv_uk2 <- krige.cv(temperature ~ Altitude + Latitude, locations = oneday_sf, model = fitted_variogram, verbose = FALSE)


# 4. Compute Metrics and Export Benchmarking Table -----------------------------
message("Step 4: Compiling statistical benchmarking metrics...")

calculate_metrics <- function(cv_obj) {
    rmse <- sqrt(mean(cv_obj$residual^2))
    mae <- mean(abs(cv_obj$residual))
    return(c(RMSE = rmse, MAE = mae))
}

ok_m <- calculate_metrics(cv_ok)
uk1_m <- calculate_metrics(cv_uk1)
uk2_m <- calculate_metrics(cv_uk2)

benchmarking_table <- data.frame(
    Model = c(
        "Ordinary Kriging (OK)",
        "Universal Kriging (UK-1: Elevation)",
        "Multivariate Universal Kriging (UK-2: Elev + Lat)"
    ),
    RMSE = c(ok_m["RMSE"], uk1_m["RMSE"], uk2_m["RMSE"]),
    MAE = c(ok_m["MAE"], uk1_m["MAE"], uk2_m["MAE"])
)

# Export the benchmarking table to a clean CSV
processed_bench_path <- "geostat-pipeline/data/processed/model_benchmarking.csv"
write_csv(benchmarking_table, processed_bench_path)
message("Benchmarking metrics CSV saved to: ", processed_bench_path)


# 5. Export Residual Error Diagnostics Bar Plot --------------------------------
message("Step 5: Exporting residual diagnostic figures...")

cv_uk2_df <- data.frame(
    Nom = oneday_sf$Nom,
    observed = cv_uk2$observed,
    predicted = cv_uk2$var1.pred,
    residual = cv_uk2$residual
)

error_plot <- ggplot(cv_uk2_df, aes(x = reorder(Nom, residual), y = residual, fill = residual > 0)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(
        values = c("TRUE" = "#2b8cbe", "FALSE" = "#de2d26"),
        labels = c("TRUE" = "Underpredicted (Actual was Hotter)", "FALSE" = "Overpredicted (Actual was Cooler)"),
        name = "Error Type"
    ) +
    coord_flip() +
    labs(
        title = "Cross-Validation Residuals by Station (June 30, 1961)",
        subtitle = "Model: Temp ~ Altitude + Latitude (LOOCV Error Diagnostics)",
        x = "Station Name",
        y = "Residual Error (°C) [Observed - Predicted]"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "#555555", size = 9),
        legend.position = "top"
    )

processed_error_fig_path <- "geostat-pipeline/output/figures/07_residual_diagnostics.png"
ggsave(
    filename = processed_error_fig_path,
    plot = error_plot, width = 8, height = 5, dpi = 300, bg = "white"
)

message("Diagnostic bar plot saved to: ", processed_error_fig_path)
message("================================================================")
message("Phase 6 Validation executed successfully.")
