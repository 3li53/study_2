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