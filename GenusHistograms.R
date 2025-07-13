# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)

# Set working directory
getwd()
setwd("C:/Users/catbb/OneDrive/Desktop/bacdive-client")

# Configuration
input_dir <- "genus_level_output"
output_dir <- "plots_output"

dir.create(output_dir, showWarnings = FALSE)

# Function to calculate cell size coverage percentage (genus level)
calculate_coverage <- function(data) {
  total_genera <- length(unique(data$genus))
  matched_genera <- length(unique(data$genus[!is.na(data$GeoMeanVolume_genus_median)]))
  coverage_percentage <- (matched_genera / total_genera) * 100
  return(coverage_percentage)
}

# Function to create histogram (genus level)
create_histogram <- function(data, study_id, file_name) {
  matched_data <- data %>% 
    filter(!is.na(GeoMeanVolume_genus_median))
  if(nrow(matched_data) == 0) {
    return(NULL)
  }
  p <- ggplot(matched_data, aes(x = log10(GeoMeanVolume_genus_median))) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "black") +
    labs(
      title = paste(file_name, study_id),
      x = "Log10(Geometric Mean Volume)",
      y = "Number of Genera"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  return(p)
}

# Function to create scatter plot (genus level)
create_scatter_plot <- function(data, study_id, file_name) {
  matched_data <- data %>% 
    filter(!is.na(GeoMeanVolume_genus_median),
           !is.na(relative_abundance),
           relative_abundance > 0,
           GeoMeanVolume_genus_median > 0,
           !is.infinite(GeoMeanVolume_genus_median),
           !is.infinite(relative_abundance))
  if(nrow(matched_data) < 2) {
    return(NULL)
  }
  matched_data <- matched_data %>%
    mutate(
      log_volume = log10(GeoMeanVolume_genus_median),
      log_abundance = log10(relative_abundance)
    ) %>%
    filter(!is.na(log_volume), 
           !is.na(log_abundance),
           !is.infinite(log_volume),
           !is.infinite(log_abundance))
  if(nrow(matched_data) < 2) {
    return(NULL)
  }
  tryCatch({
    model <- lm(log_abundance ~ log_volume, data = matched_data)
    slope <- coef(model)[2]
    r_squared <- summary(model)$r.squared
    p <- ggplot(matched_data, aes(x = log_volume, y = log_abundance)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(
        title = paste(file_name, study_id),
        subtitle = paste("Slope:", round(slope, 3), "| R²:", round(r_squared, 3)),
        x = "Log10(Geometric Mean Volume)",
        y = "Log10(Relative Abundance)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      )
    return(p)
  }, error = function(e) {
    cat("Error creating scatter plot for", study_id, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to process a single file (genus level)
process_file <- function(file_path) {
  data <- read_csv(file_path, show_col_types = FALSE, na = c("", "NA", "NaN", "Inf", "-Inf"))
  file_name <- tools::file_path_sans_ext(basename(file_path))
  required_cols <- c("genus", "study_accession", "read_count", "GeoMeanVolume_genus_median")
  missing_cols <- setdiff(required_cols, colnames(data))
  if(length(missing_cols) > 0) {
    cat("Skipping", file_name, "- Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }
  data <- data %>%
    filter(!is.na(genus), 
           !is.na(study_accession), 
           !is.na(read_count),
           read_count > 0,
           !is.infinite(read_count))
  if(nrow(data) == 0) {
    cat("Skipping", file_name, "- No valid data after cleaning\n")
    return(NULL)
  }
  studies <- unique(data$study_accession)
  file_coverage_info <- list()
  for(study_id in studies) {
    study_data <- data %>% filter(study_accession == study_id)
    coverage_percentage <- calculate_coverage(study_data)
    file_coverage_info[[study_id]] <- coverage_percentage
    clean_study_id <- gsub("[^A-Za-z0-9]", "_", study_id)
    hist_plot <- create_histogram(study_data, study_id, file_name)
    if(!is.null(hist_plot)) {
      hist_filename <- paste0(file_name, "_", clean_study_id, "_histogram.png")
      hist_path <- file.path(output_dir, hist_filename)
      ggsave(hist_path, hist_plot, width = 10, height = 6, dpi = 300)
      cat("Saved histogram:", hist_filename, "\n")
    }
    scatter_plot <- create_scatter_plot(study_data, study_id, file_name)
    if(!is.null(scatter_plot)) {
      scatter_filename <- paste0(file_name, "_", clean_study_id, "_scatter.png")
      scatter_path <- file.path(output_dir, scatter_filename)
      ggsave(scatter_path, scatter_plot, width = 10, height = 6, dpi = 300)
      cat("Saved scatter plot:", scatter_filename, "\n")
    }
  }
  return(file_coverage_info)
}

cat("Starting plot generation...\n")
csv_files <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)
if(length(csv_files) == 0) {
  cat("No CSV files found in", input_dir, "\n")
} else {
  cat("Found", length(csv_files), "CSV files to process\n")
  all_coverage_info <- list()
  for(file_path in csv_files) {
    cat("\nProcessing:", basename(file_path), "\n")
    coverage_info <- process_file(file_path)
    if(!is.null(coverage_info)) {
      all_coverage_info[[basename(file_path)]] <- coverage_info
    }
  }
  coverage_summary <- do.call(rbind, lapply(names(all_coverage_info), function(file_name) {
    file_info <- all_coverage_info[[file_name]]
    data <- read_csv(file.path(input_dir, file_name), show_col_types = FALSE, na = c("", "NA", "NaN", "Inf", "-Inf"))
    studies <- names(file_info)
    do.call(rbind, lapply(studies, function(study_id) {
      study_data <- data[data$study_accession == study_id, ]
      filtered <- study_data[!is.na(study_data$GeoMeanVolume_genus_median), ]
      n_genera_filtered <- length(unique(filtered$genus))
      data.frame(
        file = file_name,
        study = study_id,
        coverage_percentage = file_info[[study_id]],
        n_genera_filtered = n_genera_filtered,
        stringsAsFactors = FALSE
      )
    }))
  }))
  write_csv(coverage_summary, file.path(output_dir, "coverage_summary.csv"))
  coverage_bins <- list(
    "coverage_0_20" = c(0, 20),
    "coverage_20_40" = c(20, 40),
    "coverage_40_60" = c(40, 60),
    "coverage_60_80" = c(60, 80),
    "coverage_80_100" = c(80, 100.0001)
  )
  for (folder in names(coverage_bins)) {
    dir.create(file.path(output_dir, folder), showWarnings = FALSE)
  }
  for (i in seq_len(nrow(coverage_summary))) {
    cov <- coverage_summary$coverage_percentage[i]
    file_base <- tools::file_path_sans_ext(coverage_summary$file[i])
    study <- coverage_summary$study[i]
    clean_study_id <- gsub("[^A-Za-z0-9]", "_", study)
    for (j in seq_along(coverage_bins)) {
      folder <- names(coverage_bins)[j]
      bin <- coverage_bins[[j]]
      is_last_bin <- (j == length(coverage_bins))
      if ((cov >= bin[1] && cov < bin[2]) || (is_last_bin && cov == 100)) {
        hist_file <- paste0(file_base, "_", clean_study_id, "_histogram.png")
        scatter_file <- paste0(file_base, "_", clean_study_id, "_scatter.png")
        if (file.exists(file.path(output_dir, hist_file))) {
          file.copy(
            file.path(output_dir, hist_file),
            file.path(output_dir, folder, hist_file),
            overwrite = TRUE
          )
        }
        if (file.exists(file.path(output_dir, scatter_file))) {
          file.copy(
            file.path(output_dir, scatter_file),
            file.path(output_dir, folder, scatter_file),
            overwrite = TRUE
          )
        }
        break
      }
    }
  }
}

# Load required libraries
library(png)
library(grid)
library(patchwork)
library(dplyr)
library(readr)

output_dir <- "plots_output"
coverage_folders <- c("coverage_0_20", "coverage_20_40", "coverage_40_60", "coverage_60_80", "coverage_80_100")

coverage_summary <- read_csv(file.path(output_dir, "coverage_summary.csv"), show_col_types = FALSE)

for (folder in coverage_folders) {
  folder_path <- file.path(output_dir, folder)
  bin_range <- as.numeric(unlist(regmatches(folder, gregexpr("[0-9]+", folder))))
  lower <- bin_range[1]
  upper <- bin_range[2]
  is_last_bin <- folder == "coverage_80_100"
  if (is_last_bin) {
    bin_summary <- coverage_summary %>%
      filter(coverage_percentage >= lower & coverage_percentage <= upper)
  } else {
    bin_summary <- coverage_summary %>%
      filter(coverage_percentage >= lower & coverage_percentage < upper)
  }
  top6 <- bin_summary %>%
    arrange(desc(n_genera_filtered)) %>%
    head(6)
  if (nrow(top6) == 0) next
  plot_files <- mapply(function(file, study) {
    file_base <- tools::file_path_sans_ext(file)
    clean_study_id <- gsub("[^A-Za-z0-9]", "_", study)
    file.path(folder_path, paste0(file_base, "_", clean_study_id, "_scatter.png"))
  }, top6$file, top6$study, SIMPLIFY = TRUE)
  # Only keep genus-level scatter plots
  plot_files <- plot_files[grepl("with_metadata_pooled_by_study_genus_with_cellsize_", plot_files)]
  plots <- lapply(plot_files, function(f) {
    tryCatch({
      img <- png::readPNG(f)
      grid::rasterGrob(img, interpolate = TRUE)
    }, error = function(e) {
      message(sprintf("Skipping file (not a valid PNG): %s", f))
      NULL
    })
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = 3, nrow = 2)
    ggsave(
      filename = file.path(folder_path, paste0("top6_scatter_array_genus_", folder, ".png")),
      plot = combined_plot,
      width = 18, height = 8, dpi = 300
    )
    print(paste("Saved array for", folder))
  }
}