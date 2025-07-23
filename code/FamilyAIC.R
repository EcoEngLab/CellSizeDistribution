getwd()
setwd("C:/Users/catbb/OneDrive/Desktop/bacdive-client")
library(tidyverse)

input_dir <- "nopool_family_level_output_bacteria"
output_dir <- "FamilyAICplots.Bacteria"
dir.create(output_dir, showWarnings = FALSE)

# Provided bin edges
bin_edges <- c(-2.7810013, -2.6515427, -2.52208411, -2.39262551, -2.26316691, -2.13370831,
               -2.00424972, -1.87479112, -1.74533252, -1.61587392, -1.48641533, -1.35695673,
               -1.22749813, -1.09803953, -0.96858094, -0.83912234, -0.70966374, -0.58020514,
               -0.45074655, -0.32128795, -0.19182935, -0.06237076, 0.06708784, 0.19654644,
               0.32600504, 0.45546363, 0.58492223, 0.71438083, 0.84383943, 0.97329802,
               1.10275662, 1.23221522, 1.36167382, 1.49113241, 1.62059101, 1.75004961)
bin_centers <- (head(bin_edges, -1) + tail(bin_edges, -1)) / 2

# Helper function to extract host species from filename
extract_species <- function(filename) {
  name <- sub("_family_with_cellsize\\.csv$", "", filename)
  name <- gsub("_", " ", name)
  return(name)
}

csv_files <- list.files(input_dir, pattern = "_family_with_cellsize\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  message(sprintf("No input files found in directory: %s", input_dir))
} else {
  message(sprintf("Found %d input files in directory: %s", length(csv_files), input_dir))
  print(csv_files)
}

# Create output directories for model types and slopes
# dir.create(file.path(output_base, "linear", "positive_slope"), recursive = TRUE, showWarnings = FALSE)
# dir.create(file.path(output_base, "linear", "negative_slope"), recursive = TRUE, showWarnings = FALSE)
# dir.create(file.path(output_base, "linear", "uncertain_slope"), recursive = TRUE, showWarnings = FALSE)
# dir.create(file.path(output_base, "quadratic"), recursive = TRUE, showWarnings = FALSE)
# dir.create(file.path(output_base, "no_fit"), recursive = TRUE, showWarnings = FALSE)

summary_list <- list()

for (file in csv_files) {
  # Extract host species
  filename <- basename(file)
  host_species <- extract_species(filename)
  message(sprintf("Processing file: %s (host species: %s)", filename, host_species))
  
  # Read data
  df <- read_csv(file, show_col_types = FALSE)
  
  # Gather sample columns (ending with _rel_abund)
  rel_abund_cols <- grep("_rel_abund$", names(df), value = TRUE)
  sample_ids <- sub("_rel_abund$", "", rel_abund_cols)
  
  # Reshape to long format
  long_df <- df %>%
    select(family, LogGeoMeanVolume, all_of(rel_abund_cols)) %>%
    pivot_longer(
      cols = all_of(rel_abund_cols),
      names_to = "sample_id",
      values_to = "rel_abund"
    ) %>%
    mutate(
      sample_id = sub("_rel_abund$", "", sample_id),
      host_species = host_species
    )
  
  # For each sample, bin LogGeoMeanVolume and sum rel_abund per bin
  for (sid in unique(long_df$sample_id)) {
    message(sprintf("  Plotting sample: %s", sid))
    sample_data <- long_df %>%
      filter(sample_id == sid) %>%
      mutate(LogGeoMeanVolume = as.numeric(LogGeoMeanVolume)) %>%
      filter(!is.na(LogGeoMeanVolume)) %>%  # Drop rows without LogGeoMeanVolume
      mutate(bin = cut(LogGeoMeanVolume, breaks = bin_edges, include.lowest = TRUE, labels = FALSE))
    
    # Precompute bin labels
    bin_labels <- paste0(round(head(bin_edges, -1), 2), "–", round(tail(bin_edges, -1), 2))
    # Summed rel_abund per bin
    binned <- sample_data %>%
      filter(!is.na(bin)) %>%
      group_by(bin) %>%
      summarise(rel_abund = sum(rel_abund, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(
        bin_center = bin_centers[bin],
        bin_label = bin_labels[bin]
      )
    
    # Standard bar plot style with white background, grid lines, and improved fit line
    subtitle_text <- NULL
    # Only fit models if there are at least 3 non-zero, non-NA rel_abund values
    valid_rows <- binned %>% filter(!is.na(rel_abund) & rel_abund > 0)
    # --- Model fitting code commented out below ---
    # if (nrow(valid_rows) >= 3) {
    #   lm1 <- tryCatch(lm(rel_abund ~ bin_center, data = valid_rows), error = function(e) NULL)
    #   lm2 <- tryCatch(lm(rel_abund ~ poly(bin_center, 2, raw = TRUE), data = valid_rows), error = function(e) NULL)
    #   aic1 <- if (!is.null(lm1)) AIC(lm1) else NA
    #   aic2 <- if (!is.null(lm2)) AIC(lm2) else NA
    #
    #   # Only proceed if both AICs are not NA and are finite
    #   if (!is.na(aic1) && !is.na(aic2) && is.finite(aic1) && is.finite(aic2)) {
    #     if ((aic2 - aic1) > 2) {
    #       best_model <- lm1
    #       model_type <- "Linear"
    #     } else if ((aic1 - aic2) > 2) {
    #       best_model <- lm2
    #       model_type <- "Quadratic"
    #     } else {
    #       best_model <- lm1
    #       model_type <- "Linear"
    #     }
    #     # Overlay best fit
    #     if (model_type == "Quadratic") {
    #       fit_x <- seq(min(bin_centers), max(bin_centers), length.out = 300)
    #     } else {
    #       fit_x <- bin_centers
    #     }
    #     fit_df <- data.frame(bin_center = fit_x)
    #     fit_df$fit <- predict(best_model, newdata = fit_df)
    #     # Calculate slope if linear
    #     if (model_type == "Linear" && !is.null(lm1)) {
    #       slope_val <- coef(lm1)[2]
    #       subtitle_text <- paste0("Best fit: ", model_type, ", Slope: ", signif(slope_val, 3))
    #     } else {
    #       subtitle_text <- paste0("Best fit: ", model_type)
    #     }
    #     p <- ggplot(binned, aes(x = bin_center, y = rel_abund)) +
    #       geom_col(color = "black", position = "identity") +
    #       labs(
    #         title = paste0(host_species, " (", sid, ")"),
    #         subtitle = subtitle_text,
    #         x = expression(Log[10]*" Geometric Mean Volume (μm³)"),
    #         y = "Total Relative Abundance (%)"
    #       ) +
    #       ylim(0, 1) +
    #       geom_line(data = fit_df, aes(x = bin_center, y = fit), color = "red", size = 0.7) +
    #       theme_bw()
    #   } else {
    #     p <- ggplot(binned, aes(x = bin_center, y = rel_abund)) +
    #       geom_col(color = "black", position = "identity") +
    #       labs(
    #         title = paste0(host_species, " (", sid, ")"),
    #         subtitle = "No suitable model fit",
    #         x = expression(Log[10]*" Geometric Mean Volume (μm³)"),
    #         y = "Total Relative Abundance (%)"
    #       ) +
    #       ylim(0, 1) +
    #       theme_bw()
    #   }
    # } else {
    p <- ggplot(binned, aes(x = bin_center, y = rel_abund)) +
      geom_col(color = "black", position = "identity") +
      labs(
        title = paste0(host_species, " (", sid, ")"),
        subtitle = "Model fitting disabled",
        x = expression(Log[10]*" Geometric Mean Volume (μm³)"),
        y = "Total Relative Abundance (%)"
      ) +
      # ylim(0, 1) +  # y-axis released, comment out or remove this line
      theme_bw()
    # }
    
    # --- Model and slope classification ---
    # model_class <- NA
    # slope_class <- NA
    # slope_val <- NA
    # if (!is.null(subtitle_text)) {
    #   if (grepl("Linear", subtitle_text)) {
    #     model_class <- "linear"
    #     if (exists("lm1") && !is.null(lm1)) {
    #       slope_val <- coef(lm1)[2]
    #       if (is.na(slope_val) || abs(slope_val) < 1e-8) {
    #         slope_class <- "uncertain_slope"
    #       } else if (slope_val > 0) {
    #         slope_class <- "positive_slope"
    #       } else {
    #         slope_class <- "negative_slope"
    #       }
    #     } else {
    #       slope_class <- "uncertain_slope"
    #     }
    #   } else if (grepl("Quadratic", subtitle_text)) {
    #     model_class <- "quadratic"
    #   } else {
    #     model_class <- "no_fit"
    #   }
    # } else {
    #   model_class <- "no_fit"
    # }

    # --- Determine output folder ---
    out_dir <- file.path(output_dir, "barplots_only")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    pdf_file <- file.path(out_dir, paste0(gsub(" ", "_", host_species), "_", sid, "_binned_plot.pdf"))
    message(sprintf("    Saving plot to: %s", pdf_file))
    if (nrow(binned) > 0) {
      ggsave(pdf_file, plot = p, width = 7, height = 5)
    } else {
      message(sprintf("    No data to plot for %s (%s), skipping PDF.", host_species, sid))
    }

     --- Save summary info ---
    summary_list[[length(summary_list) + 1]] <- data.frame(
      host_species = host_species,
      sample_id = sid,
      model = "barplots_only", # Changed to reflect new output folder
      slope_class = "", # No slope classification for barplots
      slope = NA # No slope value for barplots
    )
  }
}

summary_df <- do.call(rbind, summary_list)

# Model fit proportions
model_counts <- table(summary_df$model)
model_props <- model_counts / sum(model_counts)

# Linear slope proportions
linear_df <- summary_df[summary_df$model == "linear", ]
slope_counts <- table(linear_df$slope_class)
slope_props <- slope_counts / sum(slope_counts)

# Write summary CSV
#write.csv(summary_df, file = file.path(output_base, "model_fit_summary.csv"), row.names = FALSE)

# Write proportions to a separate CSV
model_prop_df <- data.frame(model = names(model_props), proportion = as.numeric(model_props))
slope_prop_df <- data.frame(slope_class = names(slope_props), proportion = as.numeric(slope_props))
#write.csv(model_prop_df, file = file.path(output_base, "model_fit_proportions.csv"), row.names = FALSE)
#write.csv(slope_prop_df, file = file.path(output_base, "linear_slope_proportions.csv"), row.names = FALSE)
