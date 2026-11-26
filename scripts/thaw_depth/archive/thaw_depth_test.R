library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

# Step 1: Read and clean data
thaw_depth <- read_csv("./data_files/thaw_depth/test_4nov2025.csv")

# Clean and create factors
thaw_depth_clean <- thaw_depth %>%
  mutate(
    avg_thaw_depth_mid = suppressWarnings(as.numeric(avg_thaw_depth_mid)),
    treatment = str_trim(treatment),  # Remove whitespace
    plant = factor(str_sub(treatment, 1, 1), levels = c("B", "S")),
    cut = factor(str_sub(treatment, 2, 2), levels = c("U", "C")),
    moisture = factor(str_sub(treatment, 3, 3), levels = c("W", "D"))
  ) %>%
  filter(!is.na(avg_thaw_depth_mid))

# Check unique values (optional)
print(unique(thaw_depth_clean$treatment))
print(unique(thaw_depth_clean$plant))
print(unique(thaw_depth_clean$cut))
print(unique(thaw_depth_clean$moisture))

# Step 2: Save as SVG
svg("thaw_depth_plot_4nov2025.svg", width = 10, height = 6)

# Step 3: Plot
ggplot(thaw_depth_clean, aes(
  x = hrs_after_incubation,
  y = avg_thaw_depth_mid,
  group = treatment
)) +
  geom_line(aes(color = plant, linetype = cut), stat = "summary", fun = mean, linewidth = 1.2) +
  geom_point(aes(color = plant, shape = moisture), stat = "summary", fun = mean, size = 4) +
  scale_y_reverse() +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Thaw Depth Over Time by Treatment",
    x = "Hours After Incubation",
    y = "Thaw Depth (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal(base_size = 14)

# Close SVG device
dev.off()

# Read data
thaw_depth <- read_csv("./data_files/thaw_depth/test_4nov2025.csv")

# Clean and create factors
thaw_depth_clean <- thaw_depth %>%
  mutate(
    thaw_depth_center = suppressWarnings(as.numeric(thaw_depth_center)),
    thaw_depth_periphery = suppressWarnings(as.numeric(thaw_depth_periphery)),
    treatment = str_trim(treatment),
    plant = factor(str_sub(treatment, 1, 1), levels = c("B", "S")),
    cut = factor(str_sub(treatment, 2, 2), levels = c("U", "C")),
    moisture = factor(str_sub(treatment, 3, 3), levels = c("W", "D"))
  )

# ---- Plot 1: Thaw Depth Center ----
center_data <- thaw_depth_clean %>% filter(!is.na(thaw_depth_center))

svg("thaw_depth_center.svg", width = 10, height = 6)
ggplot(center_data, aes(
  x = hrs_after_incubation,
  y = thaw_depth_center,
  group = treatment
)) +
  geom_line(aes(color = plant, linetype = cut), stat = "summary", fun = mean, linewidth = 1.2) +
  geom_point(aes(color = plant, shape = moisture), stat = "summary", fun = mean, size = 4) +
  scale_y_reverse() +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Thaw Depth (Center) Over Time by Treatment",
    x = "Hours After Incubation",
    y = "Thaw Depth Center (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal(base_size = 14)
dev.off()

# ---- Plot 2: Thaw Depth Periphery ----
periphery_data <- thaw_depth_clean %>% filter(!is.na(thaw_depth_periphery))

svg("thaw_depth_periphery.svg", width = 10, height = 6)
ggplot(periphery_data, aes(
  x = hrs_after_incubation,
  y = thaw_depth_periphery,
  group = treatment
)) +
  geom_line(aes(color = plant, linetype = cut), stat = "summary", fun = mean, linewidth = 1.2) +
  geom_point(aes(color = plant, shape = moisture), stat = "summary", fun = mean, size = 4) +
  scale_y_reverse() +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Thaw Depth (Periphery) Over Time by Treatment",
    x = "Hours After Incubation",
    y = "Thaw Depth Periphery (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal(base_size = 14)
dev.off()









# --- Packages ---
check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(c("dplyr", "ggplot2", "tools", "readr", "stringr", "Cairo"), check_install_load))

# --- Load data ---
data_dir   <- "./data_files/thaw_depth/"
file_pattern <- "study_2_thaw_depth_run[1-6]\\.csv"
files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)
thaw_depth_list <- lapply(files, read.table, sep = ",", header = TRUE)
all_thaw_depth <- dplyr::bind_rows(thaw_depth_list)

# --- Clean data ---
thaw_depth_clean <- all_thaw_depth %>%
  dplyr::mutate(
    treatment = stringr::str_trim(treatment),
    thaw_depth_center = dplyr::coalesce(
      suppressWarnings(as.numeric(thaw_depth_center)),
      suppressWarnings(as.numeric(dplyr::na_if(thaw_depth_center, "-")))
    ),
    hrs_after_incubation = suppressWarnings(as.numeric(hrs_after_incubation)),
    vegetation_type = factor(stringr::str_sub(treatment, 1, 1), levels = c("B", "S")),
    cut            = factor(stringr::str_sub(treatment, 2, 2), levels = c("U", "C")),
    moisture       = factor(stringr::str_sub(treatment, 3, 3), levels = c("W", "D"))
  ) %>%
  dplyr::filter(!is.na(thaw_depth_center), !is.na(hrs_after_incubation))

# --- Summarize to handle repeated measurements ---
summary_td <- thaw_depth_clean %>%
  dplyr::group_by(treatment, hrs_after_incubation, vegetation_type, cut, moisture) %>%
  dplyr::summarise(
    n  = dplyr::n(),
    mean_td = mean(thaw_depth_center, na.rm = TRUE),
    sd_td   = sd(thaw_depth_center, na.rm = TRUE),
    se_td   = sd_td / sqrt(n),
    .groups = "drop"
  )

# --- Plot: mean with SE error bars ---
svg("thaw_depth_1.svg", width = 10, height = 6)

ggplot(summary_td, aes(
  x = hrs_after_incubation,
  y = mean_td,
  group = treatment
)) +
  geom_line(aes(color = vegetation_type, linetype = cut), linewidth = 1.2) +
  geom_point(aes(color = vegetation_type, shape = moisture), size = 3.8) +
  geom_errorbar(aes(ymin = mean_td - se_td, ymax = mean_td + se_td, color = vegetation_type),
                width = 0, linewidth = 0.8) +
  scale_y_reverse() +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Thaw depth over time by treatment (mean ± SE)",
    x = "Hours After Incubation",
    y = "Thaw Depth (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal()

dev.off()












# --- Packages ---
check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(c("dplyr", "ggplot2", "tools", "readr", "stringr", "Cairo"), check_install_load))

# --- Load data ---
data_dir     <- "./data_files/thaw_depth/"
file_pattern <- "study_2_thaw_depth_run[1-6]\\.csv"
files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)
thaw_depth_list <- lapply(files, read.table, sep = ",", header = TRUE)
all_thaw_depth <- dplyr::bind_rows(thaw_depth_list)

# --- Clean & derive paired average ---
thaw_depth_clean <- all_thaw_depth %>%
  dplyr::mutate(
    # Trim treatment
    treatment = stringr::str_trim(treatment),
    
    # Coerce depth columns to numeric (treat "-" as NA if present)
    thaw_depth_center = dplyr::coalesce(
      suppressWarnings(as.numeric(thaw_depth_center)),
      suppressWarnings(as.numeric(dplyr::na_if(thaw_depth_center, "-")))
    ),
    thaw_depth_periphery = dplyr::coalesce(
      suppressWarnings(as.numeric(thaw_depth_periphery)),
      suppressWarnings(as.numeric(dplyr::na_if(thaw_depth_periphery, "-")))
    ),
    
    # Time to numeric
    hrs_after_incubation = suppressWarnings(as.numeric(hrs_after_incubation)),
    
    # Factors from treatment string
    vegetation_type = factor(stringr::str_sub(treatment, 1, 1), levels = c("B", "S")),
    cut            = factor(stringr::str_sub(treatment, 2, 2), levels = c("U", "C")),
    moisture       = factor(stringr::str_sub(treatment, 3, 3), levels = c("W", "D"))
  ) %>%
  # Drop rows missing both center and periphery OR missing time
  dplyr::filter(!is.na(hrs_after_incubation)) %>%
  # Compute paired average per observation
  dplyr::rowwise() %>%
  dplyr::mutate(
    thaw_depth_avg = mean(c(thaw_depth_center, thaw_depth_periphery), na.rm = TRUE),
    both_missing   = all(is.na(c(thaw_depth_center, thaw_depth_periphery)))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!both_missing) %>%
  dplyr::select(-both_missing)

# Optional quick checks
print(unique(thaw_depth_clean$treatment))
print(unique(thaw_depth_clean$vegetation_type))
print(unique(thaw_depth_clean$cut))
print(unique(thaw_depth_clean$moisture))

# --- Summarise across replicates at each timepoint ---
summary_td <- thaw_depth_clean %>%
  dplyr::group_by(treatment, hrs_after_incubation, vegetation_type, cut, moisture) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    mean_td = mean(thaw_depth_avg, na.rm = TRUE),
    sd_td   = sd(thaw_depth_avg, na.rm = TRUE),
    se_td   = sd_td / sqrt(n),
    .groups = "drop"
  )

# --- Plot mean ± SE of the averaged (center+periphery) depth ---
svg("thaw_depth_avg_1.svg", width = 10, height = 6)

ggplot(summary_td, aes(
  x = hrs_after_incubation,
  y = mean_td,
  group = treatment
)) +
  geom_line(aes(color = vegetation_type, linetype = cut), linewidth = 1.2) +
  geom_point(aes(color = vegetation_type, shape = moisture), size = 3.8) +
  geom_errorbar(aes(ymin = mean_td - se_td, ymax = mean_td + se_td, color = vegetation_type),
                width = 0, linewidth = 0.8) +
  facet_grid(.~vegetation_type) +
  scale_y_reverse() +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Thaw depth (center & periphery averaged) over time by treatment (mean ± SE)",
    x = "Hours After Incubation",
    y = "Thaw Depth (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal()

dev.off()
