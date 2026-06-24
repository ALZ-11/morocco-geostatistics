# ==============================================================================
# Script: 02_eda_spatial.R
# Purpose: Build spatial objects, acquire national administrative boundaries,
#          and generate baseline exploratory visualizations.
# ==============================================================================

library(tidyverse)
library(sf)
library(geodata)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Spatial ingestion & administrative boundaries
# ------------------------------------------------------------------------------
message("setting up spatial frameworks and downloading GADM boundaries...")

stations_cleaned <- read_csv("data/processed/stations_cleaned.csv", show_col_types = FALSE)

# Convert tabular coordinates to formal spatial sf point collection (WGS 84 / EPSG:4326)
# x = longitude, y = latitude
stations_sf <- st_as_sf(
    stations_cleaned,
    coords = c("Longitude", "Latitude"),
    crs = 4326
)

# download Morocco national boundary from GADM (level 0)
morocco_spat <- gadm(country = "MAR", level = 0, path = tempdir())
morocco_sf <- st_as_sf(morocco_spat)

# ------------------------------------------------------------------------------
# Spatial cartography plot
# ------------------------------------------------------------------------------
message("generating station network cartography...")

station_map <- ggplot() +
    # Layer 1: Morocco administrative boundary
    geom_sf(
        data = morocco_sf,
        fill = "#f9f6f0",
        color = "#8c8c8c",
        linewidth = 0.4
    ) +
    # Layer 2: Station points colored by Altitude
    geom_sf(
        data = stations_sf,
        aes(color = Altitude),
        size = 3.5,
        alpha = 0.9
    ) +
    # Apply continuous color scale for elevation
    scale_color_viridis_c(
        option = "viridis",
        name = "Elevation (m)"
    ) +
    labs(
        title = "Meteorological Station Networks in Morocco",
        subtitle = "Geographic distribution mapped by station elevation (m)",
        x = "Longitude",
        y = "Latitude",
        caption = "Source: GADM national borders & ECA&D station metadata"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "#555555", size = 9),
        panel.background = element_rect(fill = "#fbfbfa", color = NA),
        panel.grid.major = element_line(color = "#e8e8e8", linewidth = 0.2),
        plot.caption = element_text(size = 7, color = "#777777")
    )

# Export
ggsave(
    filename = "output/figures/01_station_map.png",
    plot = station_map,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
)
message("Station map saved to: output/figures/01_station_map.png")

# ------------------------------------------------------------------------------
# Climatological time series (Microclimate comparison)
# ------------------------------------------------------------------------------
message("compiling monthly climatological time series...")

weather_cleaned <- read_csv("data/processed/weather_cleaned.csv", show_col_types = FALSE)

# Select representative stations: Midelt (Mountain), Ouarzazate (Desert), Tangier Airport (Coast)
target_stations <- c(2156, 2160, 2144)

weather_aggregated <- weather_cleaned %>%
    # Filter for targets and representative decade (1966-1976)
    filter(
        station_id %in% target_stations,
        date >= as.Date("1966-01-01"),
        date <= as.Date("1976-12-31")
    ) %>%
    # Group by month to compute monthly mean maximum temperatures
    mutate(year_month = floor_date(date, "month")) %>%
    group_by(station_id, year_month) %>%
    summarise(mean_temp = mean(temperature), .groups = "drop") %>%
    # Merge with station names
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Generate seasonal time series plot
time_series_plot <- ggplot(weather_aggregated, aes(x = year_month, y = mean_temp, color = Nom)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    scale_color_brewer(palette = "Set1", name = "Climate Zone") +
    labs(
        title = "Historical Monthly Mean Maximum Temperature (1966 - 1976)",
        subtitle = "Comparing Mountain (Midelt), Desert (Ouarzazate), and Maritime (Tangier) Microclimates",
        x = "Year",
        y = "Monthly Mean Max Temp (°C)",
        caption = "Source: ECA&D Daily Maximum Temperature"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "#555555", size = 8.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebeb", linewidth = 0.2),
        legend.position = "top",
        plot.caption = element_text(size = 7, color = "#777777")
    )

# Export
ggsave(
    filename = "output/figures/02_climate_time_series.png",
    plot = time_series_plot,
    width = 10,
    height = 5.5,
    dpi = 300,
    bg = "white"
)
message("Historical time series saved to: output/figures/02_climate_time_series.png")

# ------------------------------------------------------------------------------
# Latitudinal temperature correlation (1961)
# ------------------------------------------------------------------------------
message("computing latitudinal cooling correlation...")

# Isolate the "golden year" 1961 (maximizing active station overlap), the golden year was retrieved by trial on console
station_means_1961 <- weather_cleaned %>%
    filter(
        date >= as.Date("1961-01-01"),
        date <= as.Date("1961-12-31")
    ) %>%
    group_by(station_id) %>%
    summarise(mean_temp_1961 = mean(temperature), .groups = "drop") %>%
    # Join with station metadata
    left_join(stations_cleaned, by = c("station_id" = "ID"))

# Generate latitudinal scatter plot with linear trend
latitude_plot <- ggplot(station_means_1961, aes(x = Latitude, y = mean_temp_1961)) +
    geom_point(color = "#31a354", size = 3.5, alpha = 0.8) +
    geom_smooth(method = "lm", color = "#de2d26", se = FALSE, linetype = "dashed", formula = y ~ x) +
    geom_text(aes(label = Nom), vjust = -1, size = 2.5, check_overlap = TRUE) +
    labs(
        title = "Temperature vs. Latitude in Morocco (1961)",
        subtitle = "Evaluating the regional latitudinal cooling gradient across 8 active stations",
        x = "Latitude (Degrees North)",
        y = "Annual Mean Temperature (°C)",
        caption = "Source: ECA&D Climatological Baseline"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "#555555", size = 8.5),
        plot.caption = element_text(size = 7, color = "#777777")
    ) +
    ylim(min(station_means_1961$mean_temp_1961) - 1.5, max(station_means_1961$mean_temp_1961) + 1.5)

# Export
ggsave(
    filename = "output/figures/03_latitude_temp_correlation.png",
    plot = latitude_plot,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
)
message("Latitudinal correlation plot saved to: output/figures/03_latitude_temp_correlation.png")
message("================================================================")
message("EDA Execution successful.")
