# ==============================================================================
# Script: 02_exploratory_analysis.R
# Purpose: Build spatial geometries, download Morocco boundaries, and generate
#          exploratory maps, climatological time series, and spatial diagnostics.
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(geodata)

# Ensure the output directory exists
dir.create("geostat-pipeline/output/figures", recursive = TRUE, showWarnings = FALSE)


# 1. Load Cleaned Datasets -----------------------------------------------------
message("Loading processed datasets...")
stations_cleaned <- read_csv("geostat-pipeline/data/processed/stations_cleaned.csv", show_col_types = FALSE)
weather_cleaned <- read_csv("geostat-pipeline/data/processed/weather_cleaned.csv", show_col_types = FALSE)


# 2. Geospatial Setup (Morocco Boundary and Station sf Objects) ----------------
message("Setting up spatial data frameworks...")

# Convert station coordinates to formal sf spatial points
stations_sf <- st_as_sf(
    stations_cleaned,
    coords = c("Longitude", "Latitude"),
    crs = 4326
)

# Download and convert Morocco national boundary
morocco_spat <- gadm(country = "MAR", level = 0, path = tempdir())
morocco_sf <- st_as_sf(morocco_spat)


# 3. Generate and Export Station Cartography ----------------------------------
message("Generating station cartography map...")

station_map <- ggplot() +
    geom_sf(data = morocco_sf, fill = "#f9f6f0", color = "#8c8c8c", linewidth = 0.4) +
    geom_sf(data = stations_sf, aes(color = Altitude), size = 3.5) +
    scale_color_viridis_c(option = "viridis", name = "Altitude (m)") +
    labs(
        title = "Meteorological Station Networks in Morocco",
        subtitle = "Geographic distribution mapped by elevation (m)",
        x = "Longitude",
        y = "Latitude",
        caption = "Source: GADM administrative borders & ECA&D station metadata"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "#555555", size = 10),
        panel.background = element_rect(fill = "#fbfbfa", color = NA),
        panel.grid.major = element_line(color = "#e8e8e8", linewidth = 0.2)
    )

ggsave(
    filename = "geostat-pipeline/output/figures/01_station_map.png",
    plot = station_map, width = 8, height = 6, dpi = 300, bg = "white"
)


# 4. Generate and Export Climatological Time Series (1966 - 1976) -------------
message("Generating historical climatology time series...")

aligned_stations <- c(2156, 2160, 2144) # Midelt, Ouarzazate, Tangier Airport

weather_aligned <- weather_cleaned %>%
    filter(
        station_id %in% aligned_stations,
        date >= as.Date("1966-01-01"),
        date <= as.Date("1966-12-31") | date >= as.Date("1967-01-01") & date <= as.Date("1976-12-31")
    ) %>%
    mutate(year_month = floor_date(date, "month")) %>%
    group_by(station_id, year_month) %>%
    summarise(mean_temp = mean(temperature), .groups = "drop") %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

aligned_ts_plot <- ggplot(weather_aligned, aes(x = year_month, y = mean_temp, color = Nom)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    scale_color_brewer(palette = "Set1", name = "Climate Zone") +
    labs(
        title = "Historical Monthly Mean Maximum Temperature (1966 - 1976)",
        subtitle = "Comparing Mountain (Midelt), Desert (Ouarzazate), and Maritime (Tangier) Microclimates",
        x = "Year",
        y = "Monthly Mean Max Temp (°C)",
        caption = "Data Source: ECA&D"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "#555555", size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebeb", linewidth = 0.2),
        legend.position = "top"
    )

ggsave(
    filename = "geostat-pipeline/output/figures/02_climate_time_series.png",
    plot = aligned_ts_plot, width = 10, height = 6, dpi = 300, bg = "white"
)


# 5. Generate and Export Latitudinal Temperature Correlation -----------------
message("Generating latitudinal cooling correlation plot...")

station_means_1961 <- weather_cleaned %>%
    filter(date >= as.Date("1961-01-01"), date <= as.Date("1961-12-31")) %>%
    group_by(station_id) %>%
    summarise(mean_temp_1961 = mean(temperature), .groups = "drop") %>%
    left_join(stations_cleaned, by = c("station_id" = "ID"))

latitude_plot_1961 <- ggplot(station_means_1961, aes(x = Latitude, y = mean_temp_1961)) +
    geom_point(color = "#31a354", size = 4, alpha = 0.8) +
    geom_smooth(method = "lm", color = "#de2d26", se = FALSE, linetype = "dashed", formula = y ~ x) +
    geom_text(aes(label = Nom), vjust = -1, size = 3, check_overlap = TRUE) +
    labs(
        title = "Temperature vs. Latitude in Morocco (1961)",
        subtitle = "Evaluating the latitudinal cooling gradient across 8 stations",
        x = "Latitude (Degrees North)",
        y = "Annual Mean Temperature (°C)"
    ) +
    theme_minimal() +
    ylim(min(station_means_1961$mean_temp_1961) - 2, max(station_means_1961$mean_temp_1961) + 2)

ggsave(
    filename = "geostat-pipeline/output/figures/03_latitude_temp_correlation.png",
    plot = latitude_plot_1961, width = 8, height = 5, dpi = 300, bg = "white"
)

message("================================================================")
message("All figures regenerated and exported successfully")
