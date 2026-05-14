# Function to check, install, and load packages
check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load required packages
invisible(lapply(c("dplyr", "ggplot2", "tools", "readr", "stringr", "Cairo"), check_install_load))

# Set folder path
folder_path <- "data_files/pipes_temperature/run1"

# Rename files with pattern: three letters + underscore + single digit
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

for (file in file_list) {
  file_name <- basename(file)
  if (grepl("^[A-Za-z]{3}_\\d\\.csv$", file_name)) {
    parts <- strsplit(file_name, "_|\\.csv")[[1]]
    new_name <- paste0(parts[1], "_0", parts[2], ".csv")
    file.rename(file, file.path(folder_path, new_name))
  }
}

# Refresh file list
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read and tag each file
data_list <- lapply(file_list, function(file) {
  tryCatch({
    df <- read.csv(file)
    df$file_source <- tools::file_path_sans_ext(basename(file))
    return(df)
  }, error = function(e) {
    warning(paste("Failed to read:", file))
    return(NULL)
  })
})
data_list <- Filter(Negate(is.null), data_list)

# Combine all data
combined_data <- bind_rows(data_list)

# Define treatment-specific start times
start_times <- list(
  SCW = "2025-10-13 05:07",
  BUW = "2025-10-13 05:07",
  BCW = "2025-10-13 05:22",
  BCD = "2025-10-13 05:22",
  SUW = "2025-10-13 05:50",
  SCD = "2025-10-13 05:40",
  SUD = "2025-10-13 06:01",
  BUD = "2025-10-13 06:01"
)

# Define experiment end time
end_time <- as.POSIXct("2025-10-19 05:00", tz = "UTC")

# Clean and trim the data
cleaned_data <- combined_data %>%
  filter(grepl("^\\d+$", X)) %>%
  mutate(
    time = as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    temperature = readr::parse_number(X1),
    row = as.integer(X),
    source_file = file_source,
    prefix = str_sub(source_file, 1, 3),
    depth = as.numeric(str_sub(source_file, -2))
  ) %>%
  filter(!is.na(depth)) %>%
  mutate(
    start_time = as.POSIXct(sapply(prefix, function(p) start_times[[p]]), tz = "UTC")
  ) %>%
  filter(time >= start_time & time <= end_time) %>%
  select(-start_time)  # Optional: remove helper column

# Ensure depth is treated as a factor for consistent coloring
cleaned_data <- cleaned_data %>%
  mutate(depth = factor(depth))


CairoSVG("temperature_profiles.svg", width = 12, height = 8)

ggplot(cleaned_data, aes(x = time, y = temperature, color = depth)) +
  geom_line() +
  facet_wrap(~ prefix, ncol = 4) +
  labs(
    title = "Temperature profiles over time by treatment",
    x = "Time",
    y = "Temperature (C)",
    color = "Depth (cm)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

dev.off()

