# Load packages
check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

required_packages <- c("dplyr", "tools", "readr", "stringr")
invisible(lapply(required_packages, check_install_load))

# Define all run folders
run_paths <- file.path("./data/", paste0("run", 1:6))

# Rename files in all runs (ABC_3.csv = ABC_03.csv)
for (folder in run_paths) {
  
  file_list <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  
  for (file in file_list) {
    file_name <- basename(file)
    
    # Match 3 letters + "_" + single digit
    if (grepl("^[A-Za-z]{3}_\\d\\.csv$", file_name)) {
      parts <- strsplit(file_name, "_|\\.csv")[[1]]     # -> c("ABC", "3")
      new_name <- paste0(parts[1], "_0", parts[2], ".csv")  # -> ABC_03.csv
      file.rename(file, file.path(folder, new_name))
    }
  }
}

# Collect all files across all runs
all_files <- unlist(lapply(run_paths, function(folder) {
  list.files(folder, pattern = "\\.csv$", full.names = TRUE)
}), use.names = FALSE)

# Load each file into a list and tag with run name
data_list <- lapply(all_files, function(file) {
  run_name <- basename(dirname(file))  # e.g. "run3"
  
  df <- tryCatch(
    read.csv(file),
    error = function(e) {
      warning(paste("Failed to read:", file))
      return(NULL)
    }
  )
  
  if (!is.null(df)) {
    df$file_source <- tools::file_path_sans_ext(basename(file))
    df$run <- run_name
  }
  
  df
})

# Drop failed reads
data_list <- Filter(Negate(is.null), data_list)

# Combine into one dataframe
combined_data <- dplyr::bind_rows(data_list)

# Rename columns and clean data
cleaned_data <- combined_data %>%
  # Remove rows where Time is NA or empty
  filter(!is.na(Time) & Time != "") %>%
  # Rename columns
  rename(
    Time = Time,               # keep the same name
    temperature = X1,          # X1 -> temperature
    treatment = file_source,   # file_source -> treatment
    run = run                  # keep run
  ) %>%
mutate(
  # Extract digits after the underscore (keeps leading zeros, e.g. "03")
  depth = stringr::str_extract(treatment, "(?<=_)\\d+"),
  # Remove the underscore and trailing digits from treatment
  treatment = stringr::str_remove(treatment, "_\\d+$")
) %>%
select(-X)                   # remove X completely

# Define time windows for each run 
run_time_limits <- list(
  run1 = list(start = "2025-10-13 05:00", end = "2025-10-19 05:00"),
  run2 = list(start = "2025-10-22 05:00", end = "2025-10-28 05:00"),
  run3 = list(start = "2025-11-02 05:00", end = "2025-11-08 05:00"),
  run4 = list(start = "2025-11-11 05:00", end = "2025-11-17 05:00"),
  run5 = list(start = "2025-11-20 05:00", end = "2025-11-26 05:00"),
  run6 = list(start = "2025-11-29 05:00", end = "2025-12-05 05:00")
)

# Convert Time column to POSIXct 
cleaned_data$Time <- as.POSIXct(cleaned_data$Time, tz = "UTC")

# Apply per-run time filtering
filtered_data <- cleaned_data %>%
  rowwise() %>%
  mutate(
    start_lim = as.POSIXct(run_time_limits[[run]]$start, tz = "UTC"),
    end_lim   = as.POSIXct(run_time_limits[[run]]$end,   tz = "UTC")
  ) %>%
  filter(Time >= start_lim & Time <= end_lim) %>%
  ungroup() %>%
  select(-start_lim, -end_lim)   # clean up helper columns


# Save the combined CSV, not recommended as the file is very big.
 # write.csv(combined_data, "combined_temperature_runs.csv", row.names = FALSE)