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

# ---- Clean data ----

thaw_depth_clean <- all_thaw_depth %>%
  mutate(
    treatment = str_trim(treatment), # Remove any extra spaces before or after the treatment code
    # Replace "-" with NA (missing value) and convert to numeric
    thin_center    = as.numeric(na_if(thaw_depth_center, "-")),
    thin_periphery = as.numeric(na_if(thaw_depth_periphery, "-")),
    thick          = as.numeric(na_if(thaw_depth_cm, "-")),
    hrs_after_incubation = as.numeric(hrs_after_incubation),# Convert incubation time to a numeric variable
    # Extract treatment components from the 3-letter treatment code
    vegetation_type = factor(str_sub(treatment, 1, 1)),# First letter = vegetation type (B = bryophyte, S = salix)
    cut = factor(str_sub(treatment, 2, 2)),# Second letter = cut treatment (U = uncut, C = cut)
    moisture = factor(str_sub(treatment, 3, 3))# Third letter = moisture treatment (W = wet, D = dry)
    ) %>%
      select(# Remove columns that are not needed
        -atmospheric_pressure,
        -soil_pipe_cm_1, -soil_pipe_cm_2, -soil_pipe_cm_3,
        -pipe_lid_cm_1, -pipe_lid_cm_2, -pipe_lid_cm_3,
        -C13_pipe_lid_cm_1, -C13_pipe_lid_cm_2, -C13_pipe_lid_cm_3,
        -X
      )  %>%
  # Keep only rows with valid thaw depth and incubation time values, the rest of the rows refer to lid measurements
  filter(
    if_any(c(thaw_depth_center, thaw_depth_periphery, thaw_depth_cm),
           ~ !is.na(.) & . != "-")
  )

# check unique values
print(unique(thaw_depth_clean$treatment))
print(unique(thaw_depth_clean$vegetation_type))
print(unique(thaw_depth_clean$cut))
print(unique(thaw_depth_clean$moisture))

# change hour 97 to 96
thaw_depth_clean$hrs_after_incubation[thaw_depth_clean$hrs_after_incubation == "97"] <- "96"

# remove hour 6, 30, 54


# change back to int

#create a calibration dataframe
calibration <- data.frame(
  hrs_after_incubation = thaw_depth_clean$hrs_after_incubation,
  thick = thaw_depth_clean$thick,
  thin_center = thaw_depth_clean$thin_center,
  thin_periphery = thaw_depth_clean$thin_periphery,
  run = thaw_depth_clean$run)

calibration <- calibration %>% # remove unwanted rows
  filter(run == 3) # keep only the overlapping run
calibration <- calibration %>% 
  mutate(thick = as.numeric(thick))

# fit the regression (center points)
model_center <- lm(thin_center ~ thick, data = calibration)
summary(model_center)
# plot the regression
plot(calibration$thick,
     calibration$thin_center,
     xlab = "Thick probe",
     ylab = "Thin probe",
     pch = 16)

abline(model_center, col = "red", lwd = 2) # thin_center = 1.107 + 0.271 x thick

# fit the regression (peripheral points)
model_periphery <- lm(thin_periphery ~ thick, data = calibration)
summary(model_periphery)
# plot the regression
plot(calibration$thick,
     calibration$thin_periphery,
     xlab = "Thick probe",
     ylab = "Thin probe",
     pch = 16)

abline(model_periphery, col = "red", lwd = 2) # thin_periphery = 1.232 + 0.298 x thick

# iterate calibration

missing_center <- is.na(thaw_depth_clean$thin_center)
missing_periphery <- is.na(thaw_depth_clean$thin_periphery)

thaw_depth_clean$thin_center[missing_center] <-
  predict(model_center, newdata = thaw_depth_clean[missing_center, ])
thaw_depth_clean$thin_periphery[missing_periphery] <-
  predict(model_periphery, newdata = thaw_depth_clean[missing_periphery, ])

thaw_depth_clean <- thaw_depth_clean %>% 
  select(-c(thaw_depth_cm, thaw_depth_center, thaw_depth_periphery, thick))

write.csv(thaw_depth_clean, "data/clean/thaw_depth/thaw_depth_clean.csv", row.names = FALSE)
