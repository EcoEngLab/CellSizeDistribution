# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)

# Set working directory
getwd()
setwd("C:/Users/catbb/OneDrive/Desktop/bacdive-client")

# Configuration
input_dir <- "family_level_output"
output_dir <- "plots_output"

# Create output directories
dir.create(output_dir, showWarnings = FALSE)

# Function to calculate cell size coverage percentage
calculate_coverage <- function(data) {
  total_families <- length(unique(data$family))
  matched_families <- length(unique(data$family[!is.na(data$GeoMeanVolume_family_median)]))
  coverage_percentage <- (matched_families / total_families) * 100
  return(coverage_percentage)
}

# Function to create histogram
create_histogram <- function(data, study_id, file_name) {
  # Filter for matched families only
  matched_data <- data %>% 
    filter(!is.na(GeoMeanVolume_family_median))
  
  if(nrow(matched_data) == 0) {
    return(NULL)
  }
  
  # Use only host species name and study accession for title
  p <- ggplot(matched_data, aes(x = log10(GeoMeanVolume_family_median))) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "black") +
    labs(
      title = paste(file_name, study_id),
      x = "Log10(Geometric Mean Volume)",
      y = "Number of Families"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

# Function to create scatter plot
create_scatter_plot <- function(data, study_id, file_name) {
  # Filter for matched families only and clean data
  matched_data <- data %>% 
    filter(!is.na(GeoMeanVolume_family_median),
           !is.na(relative_abundance),
           relative_abundance > 0,
           GeoMeanVolume_family_median > 0,
           !is.infinite(GeoMeanVolume_family_median),
           !is.infinite(relative_abundance))
  
  if(nrow(matched_data) < 2) {
    return(NULL)
  }
  
  # Check for valid log values
  matched_data <- matched_data %>%
    mutate(
      log_volume = log10(GeoMeanVolume_family_median),
      log_abundance = log10(relative_abundance)
    ) %>%
    filter(!is.na(log_volume), 
           !is.na(log_abundance),
           !is.infinite(log_volume),
           !is.infinite(log_abundance))
  
  if(nrow(matched_data) < 2) {
    return(NULL)
  }
  
  # Fit regression line
  tryCatch({
    model <- lm(log_abundance ~ log_volume, data = matched_data)
    slope <- coef(model)[2]
    r_squared <- summary(model)$r.squared
    
    # Use only host species name and study accession for title, slope and R2 in subtitle
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

# Function to process a single file
process_file <- function(file_path) {
  # Read the CSV file with better error handling
  data <- read_csv(file_path, show_col_types = FALSE, na = c("", "NA", "NaN", "Inf", "-Inf"))
  
  # Get file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Check if required columns exist
  required_cols <- c("family", "study_accession", "read_count", "GeoMeanVolume_family_median")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if(length(missing_cols) > 0) {
    cat("Skipping", file_name, "- Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Clean the data - remove rows with NA or infinite values in key columns
  data <- data %>%
    filter(!is.na(family), 
           !is.na(study_accession), 
           !is.na(read_count),
           read_count > 0,
           !is.infinite(read_count))
  
  if(nrow(data) == 0) {
    cat("Skipping", file_name, "- No valid data after cleaning\n")
    return(NULL)
  }
  
  # Get unique studies in this file
  studies <- unique(data$study_accession)
  
  # Store coverage information for this file
  file_coverage_info <- list()
  
  # Process each study
  for(study_id in studies) {
    # Filter data for this study
    study_data <- data %>% filter(study_accession == study_id)
    
    # Calculate coverage for this study
    coverage_percentage <- calculate_coverage(study_data)
    file_coverage_info[[study_id]] <- coverage_percentage
    
    # Create clean study ID for filename
    clean_study_id <- gsub("[^A-Za-z0-9]", "_", study_id)
    
    # Create histogram
    hist_plot <- create_histogram(study_data, study_id, file_name)
    if(!is.null(hist_plot)) {
      hist_filename <- paste0(file_name, "_", clean_study_id, "_histogram.png")
      hist_path <- file.path(output_dir, hist_filename)
      ggsave(hist_path, hist_plot, width = 10, height = 6, dpi = 300)
      cat("Saved histogram:", hist_filename, "\n")
    }
    
    # Create scatter plot
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

# Main execution
cat("Starting plot generation...\n")

# Get all CSV files from input directory
csv_files <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)

if(length(csv_files) == 0) {
  cat("No CSV files found in", input_dir, "\n")
} else {
  cat("Found", length(csv_files), "CSV files to process\n")
  
  # Process each file
  all_coverage_info <- list()
  
  for(file_path in csv_files) {
    cat("\nProcessing:", basename(file_path), "\n")
    coverage_info <- process_file(file_path)
    if(!is.null(coverage_info)) {
      all_coverage_info[[basename(file_path)]] <- coverage_info
    }
  }
  
  # Save coverage summary
  coverage_summary <- do.call(rbind, lapply(names(all_coverage_info), function(file_name) {
    file_info <- all_coverage_info[[file_name]]
    # Read the file again to get per-study stats
    data <- read_csv(file.path(input_dir, file_name), show_col_types = FALSE, na = c("", "NA", "NaN", "Inf", "-Inf"))
    studies <- names(file_info)
    do.call(rbind, lapply(studies, function(study_id) {
      study_data <- data[data$study_accession == study_id, ]
      # Filter for families with non-NA GeoMeanVolume_family_median
      filtered <- study_data[!is.na(study_data$GeoMeanVolume_family_median), ]
      n_families_filtered <- length(unique(filtered$family))
      data.frame(
        file = file_name,
        study = study_id,
        coverage_percentage = file_info[[study_id]],
        n_families_filtered = n_families_filtered,
        stringsAsFactors = FALSE
      )
    }))
  }))

write_csv(coverage_summary, file.path(output_dir, "coverage_summary.csv"))

# Define coverage bins and output folders
coverage_bins <- list(
  "coverage_0_20" = c(0, 20),
  "coverage_20_40" = c(20, 40),
  "coverage_40_60" = c(40, 60),
  "coverage_60_80" = c(60, 80),
  "coverage_80_100" = c(80, 100.0001) # 100 inclusive
)

# Create folders if they don't exist
for (folder in names(coverage_bins)) {
  dir.create(file.path(output_dir, folder), showWarnings = FALSE)
}

# Copy plots to appropriate folders based on coverage
for (i in seq_len(nrow(coverage_summary))) {
  cov <- coverage_summary$coverage_percentage[i]
  file_base <- tools::file_path_sans_ext(coverage_summary$file[i])
  study <- coverage_summary$study[i]
  clean_study_id <- gsub("[^A-Za-z0-9]", "_", study)
  
  # Find which bin this coverage falls into (upper bound exclusive except last bin)
  for (j in seq_along(coverage_bins)) {
    folder <- names(coverage_bins)[j]
    bin <- coverage_bins[[j]]
    is_last_bin <- (j == length(coverage_bins))
    if ((cov >= bin[1] && cov < bin[2]) || (is_last_bin && cov == 100)) {
      # Histogram and scatter plot filenames
      hist_file <- paste0(file_base, "_", clean_study_id, "_histogram.png")
      scatter_file <- paste0(file_base, "_", clean_study_id, "_scatter.png")
      # Copy if file exists
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
      break # Only copy to one folder
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

# Read the coverage summary
coverage_summary <- read_csv(file.path(output_dir, "coverage_summary.csv"), show_col_types = FALSE)

for (folder in coverage_folders) {
  folder_path <- file.path(output_dir, folder)
  
  # Get the relevant rows for this coverage bin
  bin_range <- as.numeric(unlist(regmatches(folder, gregexpr("[0-9]+", folder))))
  lower <- bin_range[1]
  upper <- bin_range[2]
  is_last_bin <- folder == "coverage_80_100"
  
  # Filter coverage_summary for this bin
  if (is_last_bin) {
    bin_summary <- coverage_summary %>%
      filter(coverage_percentage >= lower & coverage_percentage <= upper)
  } else {
    bin_summary <- coverage_summary %>%
      filter(coverage_percentage >= lower & coverage_percentage < upper)
  }
  
  # Pick top 6 by n_families_filtered
  top6 <- bin_summary %>%
    arrange(desc(n_families_filtered)) %>%
    head(6)
  
  if (nrow(top6) == 0) next
  
  # Get the corresponding scatter plot file names
  plot_files <- mapply(function(file, study) {
    file_base <- tools::file_path_sans_ext(file)
    clean_study_id <- gsub("[^A-Za-z0-9]", "_", study)
    file.path(folder_path, paste0(file_base, "_", clean_study_id, "_scatter.png"))
  }, top6$file, top6$study, SIMPLIFY = TRUE)
  
  # Read the images as grobs
  plots <- lapply(plot_files, function(f) {
    if (file.exists(f)) {
      img <- png::readPNG(f)
      grid::rasterGrob(img, interpolate = TRUE)
    } else {
      NULL
    }
  })
  plots <- Filter(Negate(is.null), plots)
  
  # Combine into a 2x3 grid
  if (length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = 3, nrow = 2)
    # Save the combined plot
    ggsave(
      filename = file.path(folder_path, paste0("top6_scatter_array_", folder, ".png")),
      plot = combined_plot,
      width = 18, height = 8, dpi = 300
    )
    print(paste("Saved array for", folder))
  }
}