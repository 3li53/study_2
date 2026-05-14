# ---- Study 2 – CO₂ respiration & flux calculations ----
# Elise Blum, November 2025

# -------------------- Package handling --------------------

check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

pkgs <- c("dplyr", "ggplot2", "hms", "viridis", "stringr", 
          "plotly", "readr")
invisible(lapply(pkgs, check_install_load))


# -------------------- Load CO₂ data --------------------

data_dir <- "./data_files/EGM/"
file_pattern <- "study_2_EGM_run[1-6]\\.csv"

EGM_files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)
EGM_list  <- lapply(EGM_files, read.table, sep = ";", header = TRUE)
all_EGM   <- bind_rows(EGM_list)

# Clean 
EGM_data <- all_EGM %>%
  filter(minute != "ppm _room") %>% 
  transmute(
    run        = as.factor(run),
    date       = as.Date(date, format = "%d.%m.%Y"),
    hrs_after_incubation = as.factor(hrs_after_incubation),
    minute     = as.numeric(str_extract(minute, "\\d+")),
    treatment  = as.factor(treatment),
    ppm        = as.numeric(ppm)
  )

# Check problematic ppm drops between minute 0 and 1
ppm_check <- EGM_data %>%
  filter(minute %in% c(0, 1)) %>%
  arrange(run, treatment, hrs_after_incubation, minute) %>%
  group_by(run, treatment, hrs_after_incubation) %>%
  summarise(
    ppm_0 = ppm[minute == 0],
    ppm_1 = ppm[minute == 1],
    .groups = "drop"
  ) %>%
  filter(ppm_0 > ppm_1)

print(ppm_check) # what to do with these?

# Diagnostic plot ppm change over time
EGM_data <- EGM_data %>%
  group_by(run) %>%
  mutate(ppm_centered = ppm - first(ppm)) %>%
  ungroup()

p_diag <- ggplot(EGM_data, aes(minute, ppm_centered, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ hrs_after_incubation) +
  labs(
    title = "Diagnostic: CO₂ Accumulation Over Time (Centered)",
    x = "Minute",
    y = "Δ ppm from start"
  ) +
  theme_minimal()

ggplotly(p_diag)


# -------------------- Load chamber height data --------------------

chamber_dir <- "./data_files/thaw_depth"
chamber_pattern <- "study_2_thaw_depth_run[1-6]\\.csv"

chamber_files <- list.files(path = chamber_dir, 
                            pattern = chamber_pattern, full.names = TRUE)

chamber_list <- lapply(chamber_files, read.table, sep = ",", header = TRUE)
all_chamber  <- bind_rows(chamber_list)

# Convert height fields to numeric
height_cols <- c("soil_pipe_cm_1", "soil_pipe_cm_2", "soil_pipe_cm_3",
                 "pipe_lid_cm_1",  "pipe_lid_cm_2",  "pipe_lid_cm_3")

all_chamber <- all_chamber %>%
  mutate(across(all_of(height_cols), as.numeric))

# Compute chamber heights
chamber_heights <- all_chamber %>%
  group_by(treatment, run, hrs_after_incubation) %>%
  summarise(
    soil_to_pipe = mean(c(soil_pipe_cm_1, soil_pipe_cm_2, soil_pipe_cm_3)),
    pipe_to_lid  = mean(c(pipe_lid_cm_1,  pipe_lid_cm_2,  pipe_lid_cm_3)),
    chamber_height_m = (soil_to_pipe + pipe_to_lid) / 100,
    .groups = "drop"
  )

# Add chamber diameter (run-specific)
chamber_heights <- chamber_heights %>%
  mutate(
    diameter_m = ifelse(run %in% c(4, 5), 0.163, 0.17),
    area_m2    = pi * (diameter_m / 2)^2,
    volume_m3  = area_m2 * chamber_height_m
  )

# Ensure datatypes match flux data
chamber_heights <- chamber_heights %>%
  mutate(across(c(treatment, hrs_after_incubation, run), as.factor))

# -------------------- Flux calculations --------------------

# ppm increase per minute (linear slope)
flux_slopes <- EGM_data %>%
  mutate(minute = as.numeric(minute)) %>%
  group_by(treatment, hrs_after_incubation, run) %>%
  summarise(
    slope_ppm_min = coef(lm(ppm ~ minute))[2],
    .groups = "drop"
  )

# Convert ppm/min → mg CO₂/m³/hour (1 ppm = 1.96 mg/m³)
flux_slopes <- flux_slopes %>%
  mutate(delta_mg_m3_hr = slope_ppm_min * 60 * 1.96)

# Merge with chamber geometry
flux_data <- flux_slopes %>%
  left_join(chamber_heights, 
            by = c("treatment", "hrs_after_incubation", "run")) %>%
  mutate(
    flux_mg_m2_hr = (delta_mg_m3_hr * volume_m3) / area_m2
  )

str(flux_data)

# -------------------- Summaries & Plots --------------------

# Treatment-only summary
flux_summary <- flux_data %>%
  group_by(treatment) %>%
  summarise(
    mean_flux = mean(flux_mg_m2_hr, na.rm = TRUE),
    se_flux   = sd(flux_mg_m2_hr, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

print(flux_summary)


p1 <- ggplot(flux_summary, aes(treatment, mean_flux, fill = treatment)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_flux - se_flux,
                    ymax = mean_flux + se_flux), width = 0.2) +
  labs(
    title = "Average CO₂ Flux by Treatment",
    y = "Flux (mg CO₂ m⁻² hr⁻¹)"
  ) +
  theme_minimal()

ggplotly(p1)

# Treatment × time summary
flux_summary_time <- flux_data %>%
  group_by(treatment, hrs_after_incubation) %>%
  summarise(
    mean_flux = mean(flux_mg_m2_hr, na.rm = TRUE),
    se_flux   = sd(flux_mg_m2_hr, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

print(flux_summary_time)

p2 <- ggplot(flux_summary_time,
             aes(hrs_after_incubation, mean_flux, fill = treatment)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_flux - se_flux,
                    ymax = mean_flux + se_flux),
                width = 0.2, position = position_dodge(0.9)) +
  labs(
    title = "CO₂ Flux by Time and Treatment",
    x = "Hours After Incubation",
    y = "Flux (mg CO₂ m⁻² hr⁻¹)"
  ) +
  theme_minimal()

ggplotly(p2)
