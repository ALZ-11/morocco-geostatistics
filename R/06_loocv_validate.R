# ==============================================================================
# Script: 06_loocv_validate.R
# Purpose: Execute Leave-One-Out Cross-Validation (LOOCV) comparing Ordinary
#          Kriging vs. Regression-Kriging, and compile diagnostic figures.
# ==============================================================================

library(tidyverse)
library(sf)
library(gstat)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Part 1: Load spatial residuals and execute OK LOOCV
# ------------------------------------------------------------------------------
message("loading spatial residuals and executing OK LOOCV...")

oneday_sf <- readRDS("data/processed/oneday_residuals_sf.rds")
prescribed_vgm <- readRDS("data/processed/prescribed_variogram.rds")

# Run native Leave-One-Out Cross-Validation on the raw temperatures (OK Baseline)
# nfold = nrow() explicitly guarantees leave-one-out cross-validation.
cv_ok <- krige.cv(
    formula = temperature ~ 1,
    locations = oneday_sf,
    model = prescribed_vgm,
    nfold = nrow(oneday_sf)
)

# ------------------------------------------------------------------------------
# Custom Regression-Kriging (RK) LOOCV
# ------------------------------------------------------------------------------
message("executing custom Regression-Kriging LOOCV loop...")

rk_predicted <- numeric(nrow(oneday_sf))
rk_residuals <- numeric(nrow(oneday_sf))

# Custom LOOCV Loop
for (i in 1:nrow(oneday_sf)) {
    # Leave out station 'i' as the testing set
    cv_train <- oneday_sf[-i, ]
    cv_test <- oneday_sf[i, ]

    # Krige-interpolate training station residuals to predict at testing station 'i' (Ordinary Kriging step strictly on the anomalies)
    kriged_residual_test <- krige(
        formula = residual ~ 1,
        locations = cv_train,
        newdata = cv_test,
        model = prescribed_vgm
    )

    # Reconstruct final prediction: Expected climatological trend + predicted anomaly
    rk_predicted[i] <- cv_test$expected_temp + kriged_residual_test$var1.pred

    # Compute error (observed - predicted)
    rk_residuals[i] <- cv_test$temperature - rk_predicted[i]
}

# ------------------------------------------------------------------------------
# Calculate accuracy metrics and export benchmarks
# ------------------------------------------------------------------------------
message("calculating statistical validation benchmarks...")

# Calculate OK Metrics
rmse_ok <- sqrt(mean(cv_ok$residual^2))
mae_ok <- mean(abs(cv_ok$residual))

# Calculate RK Metrics
rmse_rk <- sqrt(mean(rk_residuals^2))
mae_rk <- mean(abs(rk_residuals))

# Compile comparative benchmarking table
benchmarking_table <- data.frame(
    Model = c("Ordinary Kriging (OK)", "Regression-Kriging (RK)"),
    RMSE = c(rmse_ok, rmse_rk),
    MAE = c(mae_ok, mae_rk)
)

processed_bench_path <- "data/processed/model_benchmarking.csv"
write_csv(benchmarking_table, processed_bench_path)

message("\n=================== model benchmarking table ===================")
print(benchmarking_table)
message("================================================================")

# ------------------------------------------------------------------------------
# Export residual error diagnostics plot
# ------------------------------------------------------------------------------
message("generating station-by-station residual error diagnostics plot...")

# Compile cross-validation dataframe for RK
cv_rk_df <- data.frame(
    Nom = oneday_sf$Nom,
    observed = oneday_sf$temperature,
    predicted = rk_predicted,
    residual = rk_residuals
)

error_plot <- ggplot(cv_rk_df, aes(x = reorder(Nom, residual), y = residual, fill = residual > 0)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(
        values = c("TRUE" = "#2b8cbe", "FALSE" = "#de2d26"),
        labels = c("TRUE" = "Underpredicted (Actual was Hotter)", "FALSE" = "Overpredicted (Actual was Cooler)"),
        name = "Error Type"
    ) +
    coord_flip() +
    labs(
        title = "Cross-Validation Residuals by Station (June 30, 1961)",
        subtitle = "Model: Regression-Kriging (LOOCV Error Diagnostics)",
        x = "Station Name",
        y = "Residual Error (°C) [Observed - Predicted]",
        caption = "Source: GADM national borders & ECA&D daily weather station observations"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "#555555", size = 8.5),
        legend.position = "top",
        plot.caption = element_text(size = 7, color = "#777777")
    )

processed_error_fig_path <- "output/figures/07_residual_diagnostics.png"
ggsave(
    filename = processed_error_fig_path,
    plot = error_plot, width = 8, height = 5.5, dpi = 300, bg = "white"
)

message("bar plot successfully saved to: ", processed_error_fig_path)
message("================================================================")
message("Execution successful.")
