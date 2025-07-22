# Load required libraries
library(tidyverse)
library(fs)

# Set input directory
input_dir <- "nopool_family_level_output_bacteria"

# List all CSV files in the directory
csv_files <- dir_ls(input_dir, regexp = "\\.csv$")

# Function to read and reshape a single file
df_list <- lapply(csv_files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  # Extract host species from file name (remove directory, extension, and suffix, replace underscores with spaces)
  host_species <- fs::path_file(file) %>%
    fs::path_ext_remove() %>%
    str_replace("_L[0-9]+_family_with_cellsize$", "") %>%
    str_replace("_family_with_cellsize$", "") %>%
    str_replace_all("_", " ")
  # Identify sample columns (ending with _rel_abund and starting with SRR/ERR/DRR)
  sample_cols <- grep("^(SRR|ERR|DRR).+_rel_abund$", names(df), value = TRUE)
  if (length(sample_cols) == 0) return(NULL)
  # Reshape to long format
  df_long <- df %>%
    select(family, LogGeoMeanVolume, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample_id",
      values_to = "rel_abund"
    ) %>%
    mutate(
      sample_id = str_replace(sample_id, "_rel_abund$", ""),
      host_species = host_species
    )
  return(df_long)
})

# Combine all data frames, remove NULLs
all_data <- bind_rows(df_list)

# Remove rows with missing LogGeoMeanVolume or rel_abund
all_data <- all_data %>% filter(!is.na(LogGeoMeanVolume), !is.na(rel_abund))

# Define global histogram bins for LogGeoMeanVolume
logvol_min <- min(all_data$LogGeoMeanVolume, na.rm = TRUE)
logvol_max <- max(all_data$LogGeoMeanVolume, na.rm = TRUE)
bin_breaks <- seq(logvol_min, logvol_max, length.out = 41) # 40 bins

# For each host species and sample, plot weighted histogram
unique_hosts <- unique(all_data$host_species)
unique_samples <- unique(all_data$sample_id)

for (host in unique_hosts) {
  for (sid in unique_samples) {
    df_sample <- all_data %>% filter(host_species == host, sample_id == sid)
    if (nrow(df_sample) == 0) next
    p <- ggplot(df_sample, aes(x = LogGeoMeanVolume, weight = rel_abund)) +
      geom_histogram(breaks = bin_breaks, color = "black", fill = "skyblue", closed = "left") +
      labs(
        title = paste("Cell Size Histogram (Weighted by Relative Abundance):", sid, "in", host),
        x = "Log10(GeoMeanVolume)",
        y = "Summed Relative Abundance"
      ) +
      xlim(logvol_min, logvol_max) +
      theme_minimal()
    # Clean up host and sample names for file output
    host_clean <- str_replace_all(host, "[ /]", "_")
    sid_clean <- str_replace_all(sid, "[ /]", "_")
    ggsave(
      filename = file.path(input_dir, paste0("family_bacteria_cell_size_hist_rel_abund_", host_clean, "_", sid_clean, ".pdf")),
      plot = p,
      width = 7, height = 5
    )
  }
}

# (Optional) Overlay KDEs for selected samples (example: first 3 samples of first host)
selected_host <- unique_hosts[1]
selected_samples <- head(unique(all_data %>% filter(host_species == selected_host) %>% pull(sample_id)), 3)
if (length(selected_samples) > 1) {
  df_selected <- all_data %>% filter(host_species == selected_host, sample_id %in% selected_samples)
  p_kde <- ggplot(df_selected, aes(x = LogGeoMeanVolume, color = sample_id, weight = rel_abund)) +
    geom_density(adjust = 1.2, size = 1) +
    labs(
      title = paste0("Cell Size KDE (Weighted by Relative Abundance): ", selected_host),
      x = "Log10(GeoMeanVolume)",
      y = "Weighted Density"
    ) +
    xlim(logvol_min, logvol_max) +
    theme_minimal()
  host_clean <- str_replace_all(selected_host, "[ /]", "_")
  ggsave(
    filename = file.path(input_dir, paste0("family_bacteria_cell_size_kde_overlay_selected_samples_", host_clean, ".pdf")),
    plot = p_kde,
    width = 7, height = 5
  )
}
