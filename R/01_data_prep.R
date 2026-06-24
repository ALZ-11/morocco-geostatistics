# ==============================================================================
# Script: 01_data_prep.R
# Purpose: Extract, clean, and standardize weather station metadata and daily
#          climatological records from ECA&D source files.
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# Coordinate conversion helper (DMS to DD)
# ------------------------------------------------------------------------------

#' Convert DMS string coordinates to decimal degrees
#'
#' @param dms_str Character. DMS coordinate in "+/-DD:MM:SS" format.
#' @return Numeric. Coordinate converted to decimal degrees (DD).
convert_dms_to_dd <- function(dms_str) {
    if (is.na(dms_str) || is.null(dms_str)) {
        return(NA_real_)
    }

    # Standardize white space and trim padding
    clean_str <- trimws(dms_str)

    if (clean_str == "" || clean_str == "NA") {
        return(NA_real_)
    }

    # Detect negative sign (south lat. or west long.)
    is_negative <- grepl("^-", clean_str)

    # Strip leading + or -
    clean_str <- gsub("^[+-]", "", clean_str)

    parts <- strsplit(clean_str, ":")[[1]]

    if (length(parts) != 3) {
        warning(paste("Malformed DMS string ignored:", dms_str))
        return(NA_real_)
    }

    numeric_parts <- as.numeric(parts)

    if (any(is.na(numeric_parts))) {
        warning(paste("Non-numeric components found in DMS string:", dms_str))
        return(NA_real_)
    }

    # Compute decimal degrees: deg + min/60 + sec/3600
    dd <- numeric_parts[1] + (numeric_parts[2] / 60) + (numeric_parts[3] / 3600)

    # Reapply negative direction if applicable
    if (is_negative) {
        dd <- -dd
    }

    return(dd)
}

# vectorize function
convert_dms_to_dd_vec <- Vectorize(convert_dms_to_dd, USE.NAMES = FALSE)


# ------------------------------------------------------------------------------
# Process weather station metadata
# ------------------------------------------------------------------------------
message("parsing and standardizing station metadata...")

raw_stations_path <- "data/raw/stations.txt"
processed_stations_path <- "data/processed/stations_cleaned.csv"

if (!file.exists(raw_stations_path)) {
    stop(paste("Station metadata not found at", raw_stations_path))
}

# Read raw dataset bypassing the 17 metadata comment headers
raw_stations <- read_csv(raw_stations_path, skip = 17, show_col_types = FALSE)

# structure and convert stations spatial table
stations_cleaned <- raw_stations %>%
    # Filter out blank spacer rows or comment anomalies
    filter(!is.na(STAID)) %>%
    # remove leading/trailing space anomalies
    rename_all(trimws) %>%
    # Cast columns to standardized types and process coordinates
    mutate(
        ID = as.integer(STAID),
        Nom = trimws(STANAME),
        Country = trimws(CN),
        Altitude = as.numeric(HGHT),
        Latitude = convert_dms_to_dd_vec(LAT),
        Longitude = convert_dms_to_dd_vec(LON)
    ) %>%
    # Keep strictly standard spatial columns
    select(ID, Nom, Country, Latitude, Longitude, Altitude)

# Write metadata file to processed data directory
write_csv(stations_cleaned, processed_stations_path)
message("Station metadata successfully parsed and written to: ", processed_stations_path)


# ------------------------------------------------------------------------------
# Ingest, clean, and merge daily climatological records
# ------------------------------------------------------------------------------
message("\nprocessing and merging raw climatological records...")

# Individual file processor helper
process_single_weather_file <- function(file_path) {
    # extract integer station ID using regex from filename
    sta_id <- as.integer(stringr::str_extract(basename(file_path), "\\d+"))

    if (!file.exists(file_path)) {
        warning(paste("Target weather file not found:", file_path))
        return(NULL)
    }

    # Read raw CSV skipping exactly 19 metadata lines (header is on line 20)
    raw <- read_csv(file_path, skip = 19, show_col_types = FALSE)

    # Process and format table
    clean <- raw %>%
        # trim whitespace from headers to resolve accidental spacing anomalies
        rename_all(trimws) %>%
        # drop invalid quality flags and filter out missing data codes (-9999)
        filter(Q_TX == 0, TX != -9999) %>%
        # Format dates and convert temperature from decidegrees (0.1 °C) to standard Celsius
        mutate(
            station_id = sta_id,
            date = as.Date(as.character(DATE), format = "%Y%m%d"),
            temperature = TX * 0.1
        ) %>%
        # Select standardized variables
        select(station_id, date, temperature)

    return(clean)
}

# Discover all files matching the ECA&D station file standard
weather_files <- list.files(
    path = "data/raw",
    pattern = "^TX_STAID.*\\.txt$",
    full.names = TRUE
)

if (length(weather_files) == 0) {
    stop("No raw weather files (TX_STAID*.txt) detected in data/raw/")
}

# Map over all files and bind them together into a unified data frame
combined_weather_data <- weather_files %>%
    map(process_single_weather_file) %>%
    bind_rows()

# Write clean dataset to "processed" directory
processed_weather_path <- "data/processed/weather_cleaned.csv"
write_csv(combined_weather_data, processed_weather_path)

message("Unified daily weather database saved to: ", processed_weather_path)
message("================================================================")
message("Pipeline complete. Total rows processed: ", nrow(combined_weather_data))
