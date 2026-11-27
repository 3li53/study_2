load_or_install(c("dplyr", "ggplot2", "tools", "readr", "stringr", "Cairo"))

thaw_depth_clean = read.csv("data/clean/thaw_depth/thaw_depth_clean.csv", header = TRUE)

ggplot(thaw_depth_clean, aes(x = hrs_after_incubation, y = thin_center)) +
  geom_point()

ggplot(thaw_depth_clean, aes(x = hrs_after_incubation, y = thin_periphery)) +
  geom_point()



