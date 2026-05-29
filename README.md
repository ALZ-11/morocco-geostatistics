# Spatial and Temporal Analysis of Temperature in Morocco

This project performs a geostatistical analysis of temperature data in Morocco using R and Python. It leverages data from meteorological stations to model spatial and temporal trends, providing insights for climate resilience and resource management.

## Features
- **Data Ingestion Pipeline**: Automated parsing of ECA&D fixed-width climatic records.
- **Coordinate Conversion**: Python-based DMS to Decimal Degree transformation.
- **Spatial Modeling**: Empirical variogram computation and visualization using `geoR`.
- **Reproducible Analysis**: Complete workflow documented in an R Markdown notebook.

## Tech Stack
- **R**: `geoR`, `data.table`, `geodata`, `base R`
- **Python**: `pandas`
- **Notebooks**: R Markdown

## Setup
1. **Requirements**: 
   - R (>= 4.0)
   - Python (>= 3.8) with `pandas`
2. **Library Installation**:
   ```R
   install.packages(c("geoR", "data.table", "geodata"))
   ```
3. **Execution**:
   Open `Geostat_data and notebook/geostat.Rmd` in RStudio or any R environment and run the chunks sequentially.

## Data Source
Temperature data sourced from the **ECA&D (European Climate Assessment & Dataset)**.
