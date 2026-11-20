# Function to check, install, and load packages
load_or_install(c("dplyr", "ggplot2", "tools", "readr", "stringr", "Cairo"))

# define the directory and file pattern for thaw depth data
data_dir <- "./data/raw/thaw_depth/"
file_pattern <- "study_2_thaw_depth_run[1-6]\\.csv"

# get all matching file paths
files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)

# read all files into a list of data frames
thaw_depth_list <- lapply(files, read.table, sep = ",", header = TRUE)

# combine all into one data frame called all_EGM
all_thaw_depth <- bind_rows(thaw_depth_list)

# ---- Clean data (single, consistent step) ----
thaw_depth_clean <- all_thaw_depth %>%
  # Trim treatment spaces
  dplyr::mutate(treatment = stringr::str_trim(treatment)) %>%
  # Ensure thaw_depth_cm is numeric (handle "-" as missing if present)
  dplyr::mutate(
    thaw_depth_center = dplyr::coalesce(
      suppressWarnings(as.numeric(thaw_depth_center)),          # works if already numeric-like
      suppressWarnings(as.numeric(dplyr::na_if(thaw_depth_center, "-"))) # or replace "-" with NA then numeric
    ),
    thaw_depth_periphery = dplyr::coalesce(
      suppressWarnings(as.numeric(thaw_depth_periphery)),          # works if already numeric-like
      suppressWarnings(as.numeric(dplyr::na_if(thaw_depth_periphery, "-"))) # or replace "-" with NA then numeric
    ),
    # Ensure hrs_after_incubation is numeric too
    hrs_after_incubation = suppressWarnings(as.numeric(hrs_after_incubation)),
    # Derive factor columns from treatment
    vegetation_type   = factor(stringr::str_sub(treatment, 1, 1), levels = c("B", "S")),
    cut     = factor(stringr::str_sub(treatment, 2, 2), levels = c("U", "C")),
    moisture= factor(stringr::str_sub(treatment, 3, 3), levels = c("W", "D"))
  ) %>%
  # Keep only complete rows needed for plotting
  dplyr::filter(!is.na(thaw_depth_center), !is.na(hrs_after_incubation))


# Optional: Inspect unique values
print(unique(thaw_depth_clean$treatment))
print(unique(thaw_depth_clean$vegetation_type))
print(unique(thaw_depth_clean$cut))
print(unique(thaw_depth_clean$moisture))

# plot thw depth run 3-6 center values
ggplot(thaw_depth_clean, aes(
  x = hrs_after_incubation,
  y = thaw_depth_center,          # <-- use the numeric column
  group = treatment
)) +
  # Use stat="summary" to plot means at each time point per treatment
  geom_line(aes(color = vegetation_type, linetype = cut), stat = "summary", fun = mean, linewidth = 1) +
  scale_y_reverse() +
  geom_point(aes(color = vegetation_type, shape = moisture), stat = "summary", fun = mean, size = 4) +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Central thaw depth over time by treatment",
    x = "Hours After Incubation",
    y = "Thaw Depth (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal()

#plot thaw depth run 3-6 peripheral values
ggplot(thaw_depth_clean, aes(
  x = hrs_after_incubation,
  y = thaw_depth_periphery,          # <-- use the numeric column
  group = treatment
)) +
  # Use stat="summary" to plot means at each time point per treatment
  geom_line(aes(color = vegetation_type, linetype = cut), stat = "summary", fun = mean, linewidth = 1) +
  scale_y_reverse() +
  geom_point(aes(color = vegetation_type, shape = moisture), stat = "summary", fun = mean, size = 4) +
  scale_color_manual(values = c("B" = "#1b9e77", "S" = "#d95f02")) +
  scale_linetype_manual(values = c("U" = "solid", "C" = "dashed")) +
  scale_shape_manual(values = c("W" = 16, "D" = 17)) +
  labs(
    title = "Peripheral thaw depth over time by treatment",
    x = "Hours After Incubation",
    y = "Thaw Depth (cm)",
    color = "Plant Type",
    linetype = "Cutting",
    shape = "Moisture"
  ) +
  theme_minimal()

#create a calibration dataframe
calibration <- data.frame(
  hrs_after_incubation = thaw_depth_clean$hrs_after_incubation,
  thick = thaw_depth_clean$thaw_depth_cm,
  thin_center = thaw_depth_clean$thaw_depth_center,
  thin_periphery = thaw_depth_clean$thaw_depth_periphery)

calibration <- calibration %>% # remove unwanted rows
  filter(!is.na(thick), thick != "-", thick != "")
calibration <- calibration %>% 
  mutate(thick = as.numeric(thick))

# fit the regression (center points)
model <- lm(thin_center ~ thick, data = calibration)
summary(model)
# plot the regression
plot(calibration$thick,
     calibration$thin_center,
     xlab = "Thick probe",
     ylab = "Thin probe",
     pch = 16)

abline(model, col = "red", lwd = 2) # points are variably close to the fitted line....

# fit the regression (peripheral points)
#model <- lm(thin_periphery ~ thick, data = calibration)
#summary(model)
# plot the regression
#plot(calibration$thick,
#     calibration$thin_periphery,
#     xlab = "Thick probe",
#     ylab = "Thin probe",
#     pch = 16)

#abline(model, col = "red", lwd = 2) # points are NOT close to the fitted line



# maybe a simple offset is sufficient?
calibration$diff <- calibration$thick - calibration$thin_center

plot(calibration$thick,
     calibration$diff)#,
  
#   xlab = "Thick probe",
   #  ylab = "Difference (Thin - Thick)")
abline(h = mean(calibration$diff), col = "red")


