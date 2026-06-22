# ==============================================================================
# Script: 01_data_prep.R
# Purpose: Clean station coordinates (DMS to DD) and assemble daily weather
#          time series across 11 meteorological stations in Morocco.
# ==============================================================================

# Load required libraries
library(tidyverse)

# ------------------------------------------------------------------------------
# Part 1: Define Coordinate Conversion Helpers
# ------------------------------------------------------------------------------

convert_dms_to_dd <- function(dms_str) {
    if (is.na(dms_str)) {
        return(NA_real_)
    }

    # Trim padding whitespace
    dms_str <- trimws(dms_str)

    # Extract negative direction sign
    is_negative <- grepl("^-", dms_str)

    # Strip leading sign character
    clean_str <- gsub("^[+-]", "", dms_str)

    # Split into Degrees, Minutes, Seconds
    parts <- as.numeric(strsplit(clean_str, ":")[[1]])

    # Calculate decimal degrees
    dd <- parts[1] + (parts[2] / 60) + (parts[3] / 3600)

    # Reapply negative direction
    if (is_negative) {
        dd <- -dd
    }

    return(dd)
}

# Vectorize the function for data frame operations
convert_dms_to_dd_vec <- Vectorize(convert_dms_to_dd)


# ------------------------------------------------------------------------------
# Part 2: Process Station Metadata Table
# ------------------------------------------------------------------------------
message("Step 1: Processing station metadata...")

raw_stations_path <- "geostat-pipeline/data/raw/stations.txt"
processed_stations_path <- "geostat-pipeline/data/processed/stations_cleaned.csv"

stations_raw <- read_csv(raw_stations_path, skip = 17, show_col_types = FALSE)

stations_clean <- stations_raw %>%
    # Filter out the empty spacer line
    filter(!is.na(STAID)) %>%
    # Trim whitespace from headers
    rename_all(trimws) %>%
    # Clean data types and convert coordinates
    mutate(
        ID = as.integer(STAID),
        Nom = trimws(STANAME),
        Country = trimws(CN),
        Altitude = as.numeric(HGHT),
        Latitude = convert_dms_to_dd_vec(LAT),
        Longitude = convert_dms_to_dd_vec(LON)
    ) %>%
    select(ID, Nom, Country, Latitude, Longitude, Altitude)

# Write output
write_csv(stations_clean, processed_stations_path)
message("Station coordinates converted and saved to: ", processed_stations_path)


# ------------------------------------------------------------------------------
# Part 3: Assemble Weather Time Series Datasets
# ------------------------------------------------------------------------------
message("\nStep 2: Scanning raw directory for climate datasets...")

# Discover files matching the ECA&D station naming convention
weather_files <- list.files(
    path = "geostat-pipeline/data/raw",
    pattern = "^TX_STAID.*\\.txt$",
    full.names = TRUE
)

# Processing helper function for a single weather file
process_single_weather_file <- function(file_path) {
    # Extract the integer station ID using a regular expression
    sta_id <- as.integer(stringr::str_extract(basename(file_path), "\\d+"))

    # Skip 19 metadata comment rows to parse the column header at row 20
    raw <- read_csv(file_path, skip = 19, show_col_types = FALSE)

    clean <- raw %>%
        rename_all(trimws) %>%
        # Filter for mathematically valid records and drop missing data codes
        filter(Q_TX == 0, TX != -9999) %>%
        mutate(
            station_id = sta_id,
            date = as.Date(as.character(DATE), format = "%Y%m%d"),
            temperature = TX * 0.1
        ) %>%
        select(station_id, date, temperature)

    return(clean)
}

# Process all files and merge them together
message("Processing all weather files. This may take a moment...")
combined_weather_data <- weather_files %>%
    map(process_single_weather_file) %>%
    bind_rows()

# Write output
processed_weather_path <- "geostat-pipeline/data/processed/weather_cleaned.csv"
write_csv(combined_weather_data, processed_weather_path)

message("Combined weather database saved to: ", processed_weather_path)
message("================================================================")
message("Phase 1 Pipeline executed successfully. Total rows: ", nrow(combined_weather_data))
