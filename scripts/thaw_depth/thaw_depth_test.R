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