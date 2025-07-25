setwd("C:/Users/charmaine/OneDrive/Desktop/bacdive-client")

suppressPackageStartupMessages({
  library(data.table)  
  library(moments)      
})

# Set parameters
input_dir <- "C:/Users/charmaine/OneDrive/Desktop/bacdive-client/nopool_genus_level_output_bacteria"
output_dir <- "C:/Users/charmaine/OneDrive/Desktop/bacdive-client/shuffling_results_genus_bacteria"   
n_shuffles <- 1000
set.seed(12345)  # For reproducibility

if (!dir.exists(output_dir)) dir.create(output_dir)

# Helper function: Weighted mean
weighted_mean <- function(x, w) {
  sum(x * w) / sum(w)
}

# Helper function: Weighted variance
weighted_var <- function(x, w) {
  wm <- weighted_mean(x, w)
  sum(w * (x - wm)^2) / sum(w)
}

# Helper function: Weighted kurtosis
weighted_kurtosis <- function(x, w) {
  wm <- weighted_mean(x, w)
  wv <- weighted_var(x, w)
  n <- length(x)
  if (is.na(wv) || wv == 0) return(NA)
  m4 <- sum(w * (x - wm)^4) / sum(w)
  kurt <- m4 / (wv^2)
  return(kurt)
}

# Helper function: Generate unique shuffles
unique_shuffles <- function(x, n) {
  shuffles <- list()
  seen <- new.env(hash=TRUE, parent=emptyenv())
  count <- 0
  while (length(shuffles) < n) {
    s <- sample(x, length(x), replace=FALSE)
    key <- paste(s, collapse="|")
    if (!exists(key, envir=seen, inherits=FALSE)) {
      shuffles[[length(shuffles)+1]] <- s
      assign(key, TRUE, envir=seen)
    }
    count <- count + 1
    if (count > n * 10) break  # Prevent infinite loop if n > factorial(length(x))
  }
  return(shuffles)
}

# List all input files
input_files <- list.files(input_dir, pattern="_genus_with_cellsize\\.csv$", full.names=TRUE)

# Main loop over files
for (file in input_files) {
  # Extract host species name
  fname <- basename(file)
  host_species <- sub("_L7_genus_with_cellsize\\.csv$", "", fname)
  host_species_clean <- gsub("_", " ", host_species)

  # Read data, skipping comment lines
  dt <- as.data.table(read.csv(file, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE))

  # Identify sample columns and rel_abund columns
  sample_cols <- grep("^(SRR|ERR|DRR)", names(dt), value=TRUE)
  rel_abund_cols <- grep("^(SRR|ERR|DRR).+_rel_abund$", names(dt), value=TRUE)
  # For each sample, find the corresponding rel_abund column
  for (sample_col in sample_cols) {
    sample_id <- sample_col
    rel_abund_col <- paste0(sample_id, "_rel_abund")
    if (!(rel_abund_col %in% rel_abund_cols)) next  # Skip if rel_abund column missing

    # Subset data for this sample, remove NAs in GeoMeanVolume_genus_median
    sub_dt <- dt[!is.na(GeoMeanVolume_genus_median) & !is.na(get(sample_col)) & !is.na(get(rel_abund_col)),
                 .(genus, GeoMeanVolume_genus_median, rel_abund=get(rel_abund_col))]
    if (nrow(sub_dt) < 2) next  # Need at least 2 genera for shuffling

    # Check if enough unique shuffles are possible
    n_entries <- nrow(sub_dt)
    if (factorial(n_entries) < n_shuffles) {
      warning(sprintf("Not enough unique shuffles possible for %s %s (n=%d, n!=%d, required=%d). Skipping.", 
                      host_species, sample_id, n_entries, factorial(n_entries), n_shuffles))
      next
    }

    # Normalize relative abundances (should sum to 1, but just in case)
    sub_dt[, rel_abund := rel_abund / sum(rel_abund)]

    # Original statistics
    orig_mean <- weighted_mean(sub_dt$GeoMeanVolume_genus_median, sub_dt$rel_abund)
    orig_var <- weighted_var(sub_dt$GeoMeanVolume_genus_median, sub_dt$rel_abund)
    orig_kurt <- weighted_kurtosis(sub_dt$GeoMeanVolume_genus_median, sub_dt$rel_abund)

    # Pre-allocate result vectors
    mean_shuffled <- numeric(n_shuffles)
    var_shuffled <- numeric(n_shuffles)
    kurt_shuffled <- numeric(n_shuffles)
    t_value_mean <- numeric(n_shuffles)
    p_value_mean <- numeric(n_shuffles)

    # Generate unique shuffles
    shuffles <- unique_shuffles(sub_dt$GeoMeanVolume_genus_median, n_shuffles)
    actual_shuffles <- length(shuffles)
    if (actual_shuffles < n_shuffles) {
      warning(sprintf("Only %d unique shuffles possible for %s %s", actual_shuffles, host_species, sample_id))
    }

    # For t-test, create a vector of the original values (weighted by rel_abund)
    times_vec <- round(sub_dt$rel_abund * 1000)
    times_vec[is.na(times_vec) | times_vec < 0] <- 0
    orig_weighted <- rep(sub_dt$GeoMeanVolume_genus_median, times=times_vec)

    for (i in seq_len(actual_shuffles)) {
      shuffled <- shuffles[[i]]
      mean_shuffled[i] <- weighted_mean(shuffled, sub_dt$rel_abund)
      var_shuffled[i] <- weighted_var(shuffled, sub_dt$rel_abund)
      kurt_shuffled[i] <- weighted_kurtosis(shuffled, sub_dt$rel_abund)
      # For t-test, create weighted vector for shuffled
      shuffled_weighted <- rep(shuffled, times=times_vec)
      if (sum(times_vec) == 0) {
        t_value_mean[i] <- NA
        p_value_mean[i] <- NA
      } else {
        ttest <- tryCatch(t.test(orig_weighted, shuffled_weighted, paired=FALSE), error=function(e) list(statistic=NA, p.value=NA))
        t_value_mean[i] <- ifelse(is.null(ttest$statistic), NA, as.numeric(ttest$statistic))
        p_value_mean[i] <- ifelse(is.null(ttest$p.value), NA, as.numeric(ttest$p.value))
      }
    }

    # Prepare output data.frame (only up to actual_shuffles)
    out_df <- data.frame(
      iteration = 1:actual_shuffles,
      t_value_mean = t_value_mean[1:actual_shuffles],
      p_value_mean = p_value_mean[1:actual_shuffles],
      mean_original = orig_mean,
      mean_shuffled = mean_shuffled[1:actual_shuffles],
      mean_diff = orig_mean - mean_shuffled[1:actual_shuffles],
      variance_original = orig_var,
      variance_shuffled = var_shuffled[1:actual_shuffles],
      variance_diff = orig_var - var_shuffled[1:actual_shuffles],
      kurtosis_original = orig_kurt,
      kurtosis_shuffled = kurt_shuffled[1:actual_shuffles],
      kurtosis_diff = orig_kurt - kurt_shuffled[1:actual_shuffles]
    )

    # Output file name
    out_file <- file.path(output_dir, paste0(host_species, "_", sample_id, "_shuffling_results.csv"))
    fwrite(out_df, out_file)
  }
}

cat("Shuffling analysis complete. Results saved to:", output_dir, "\n")

# --- Empirical P-value Categorization and Output Summary ---

suppressPackageStartupMessages({
  library(ggplot2)
})

# List all result CSVs in the output directory
result_files <- list.files(output_dir, pattern="_shuffling_results\\.csv$", full.names=TRUE)

# Prepare summary data.frame
summary_list <- list()

for (f in result_files) {
  df <- tryCatch(read.csv(f, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(df) || !"p_value_mean" %in% names(df)) next
  sample_id <- sub("_shuffling_results.csv$", "", basename(f))
  # Empirical p-value: proportion of shuffles with p < 0.05
  emp_p <- mean(df$p_value_mean < 0.05, na.rm=TRUE)
  summary_list[[length(summary_list)+1]] <- data.frame(
    sample_id = sample_id,
    empirical_p = emp_p
  )
}

summary_df <- do.call(rbind, summary_list)

# Remove rows with empirical_p == 0 or is.na(empirical_p)
summary_df <- summary_df[!is.na(summary_df$empirical_p) & summary_df$empirical_p != 0, ]

# Categorize empirical p-values into bins
breaks <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.995, 1)
labels <- c(
  "0-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.995", ">0.995"
)
summary_df$category <- cut(summary_df$empirical_p, breaks=breaks, labels=labels, include.lowest=TRUE, right=FALSE)

# Output summary CSV
summary_csv <- file.path(output_dir, "empirical_pvalue_summary.csv")
write.csv(summary_df, summary_csv, row.names=FALSE)

# Prepare and append the proportions to the CSV
prop_table <- prop.table(table(summary_df$category))
prop_df <- data.frame(
  category = names(prop_table),
  proportion = as.numeric(round(prop_table, 3))
)

# Append a blank row and then the proportions
write.table(
  rbind(rep(NA, ncol(summary_df))), 
  summary_csv, 
  append=TRUE, 
  sep = ",", 
  col.names=FALSE, 
  row.names=FALSE, 
  na = ""
)
write.table(
  data.frame(category_proportion = "category_proportion"), 
  summary_csv, 
  append=TRUE, 
  sep = ",", 
  col.names=FALSE, 
  row.names=FALSE
)
write.table(
  prop_df, 
  summary_csv, 
  append=TRUE, 
  sep = ",", 
  col.names=TRUE, 
  row.names=FALSE
)

# Plot histogram
hist_pdf <- file.path(output_dir, "empirical_pvalue_histogram.pdf")
p <- ggplot(summary_df, aes(x=empirical_p)) +
  geom_histogram(binwidth=0.025, fill="skyblue", color="black", boundary=0) +
  geom_vline(xintercept=breaks, linetype="dashed", color="red") +
  scale_x_continuous(breaks=seq(0,1,0.1), limits=c(0,1)) +
  labs(title="Empirical P-value Distribution Across Samples",
       x="Empirical P-value (proportion p < 0.05)",
       y="Number of Samples") +
  theme_minimal()
ggsave(hist_pdf, p, width=8, height=5)

cat("Summary CSV and histogram saved to:", summary_csv, "and", hist_pdf, "\n")

# Print proportion of samples in each empirical p-value bin
cat("\nProportion of samples in each empirical p-value bin (category):\n")
prop_table <- prop.table(table(summary_df$category))
print(round(prop_table, 3))
