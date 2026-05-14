# function to check whether package is loaded and installs it if not. 
check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    message(paste(pkg, "is already installed and loaded."))
  }
}

# loading packages - if necessary installing them
check_install_load('ggplot2') 
check_install_load('hms') 
check_install_load('stringr')
check_install_load('viridis')
check_install_load('plotly')
check_install_load('dplyr')

#Import data
EGM_test_import <- read.table(file = "./data_files/study2_EGMdata.csv", sep = ",", header = TRUE)

#see data types
str(EGM_test_import)

#new dataframe with relevant data in correct formats
EGM_test <- EGM_test_import %>%
  transmute(
    date = as.Date(Date, format = "%m/%d/%Y"),
    time = as_hms(strptime(Time, format = "%H:%M")), 
    measurement = as.factor(rep(rep(1:19, each = 8), length.out = n())),
    minute = as.factor(str_extract(minute, "\\d+|room")),  # Extract only digits (d+) or (|) "room"
    pipe = as.factor(Pipe),
    tinytag = as.factor(Tiny.Tag),
    atm = as.numeric(AP),
    ppm = as.numeric(ppm)
  )

str(EGM_test)
summary(EGM_test)

#subtract room ppm from pipe measurement group

EGM_test <- EGM_test %>%
  group_by(pipe, measurement) %>%
  mutate(
    room_ppm = ppm[minute == "room"], # Extract room ppm for each group
    ppm_adj = ppm - room_ppm # Subtract room ppm
  ) %>%
  ungroup()

#remove room row
EGM_test_filtered <- EGM_test %>%
  filter(minute != "room")


#extract chamber heights for volume estimates

EGM_test_import$measurement <- EGM_test$measurement

EGM_test_import <- EGM_test_import %>%
  mutate(across(starts_with("cm."), as.numeric))

chamber_heights <- EGM_test_import %>%
  group_by(Pipe, measurement) %>%
  summarize(
    soil_to_pipe = mean(c(cm.soil.pipe.1[1], cm.soil.pipe.2[1], cm.soil.pipe.3[1]), na.rm = TRUE),
    pipe_to_lid  = mean(c(cm.lid.pipe.1[1], cm.lid.pipe.2[1], cm.lid.pipe.3[1]), na.rm = TRUE),
    chamber_height_cm = soil_to_pipe + pipe_to_lid,
    .groups = "drop"
  )

#convert height to meters
chamber_heights <- chamber_heights %>%
  mutate(chamber_height_m = chamber_height_cm / 100)

#calculate chamber area
diameter_m <- 0.163
area_m2 <- pi * (diameter_m / 2)^2  # Cross-sectional area of the pipe

#calculate chamber volume
chamber_heights <- chamber_heights %>%
  mutate(volume_m3 = area_m2 * chamber_height_m)

#calculate fluxes!
# Flux =  (Δ𝐶⋅𝑉)/(A ⋅ Δt) Where
#:ΔC = change in CO₂ concentration (ppm → mg/m³)
#V = chamber volume (m³)
#A = chamber area (m²) ~0.0209
#Δt = time interval (hours)

#1 ppm CO₂ ≈ 1.96 mg/m³ at standard conditions

#calculate slope of ppm_adj per minute
flux_slopes <- EGM_test %>%
  filter(minute != "room") %>%
  mutate(minute_num = as.numeric(as.character(minute))) %>%
  group_by(pipe, measurement) %>%
  summarize(
    slope_ppm_per_min = coef(lm(ppm_adj ~ minute_num))[2],
    .groups = "drop"
  )

#convert slope to mg/m3 hr-1
flux_slopes <- flux_slopes %>%
  mutate(delta_mg_m3_per_hr = slope_ppm_per_min * 60 * 1.96)

#join chamber heights
flux_data <- flux_slopes %>%
  left_join(chamber_heights, by = c("pipe" = "Pipe", "measurement")) %>%
  mutate(area_m2 = 0.0209)  # pipe cross-sectional area

#cAlCuLaTe FlUx!!
flux_data <- flux_data %>%
  mutate(flux_mg_m2_hr = (delta_mg_m3_per_hr * volume_m3) / area_m2)

#Visulise
p <- ggplot(flux_data, aes(x = measurement, y = flux_mg_m2_hr, fill = pipe)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Interactive CO₂ Flux by Measurement and Pipe",
    x = "Measurement",
    y = "Flux (mg CO₂ m⁻² hr⁻¹)",
    fill = "Pipe"
  )

ggplotly(p)

