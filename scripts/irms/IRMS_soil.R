load_or_install(c("tidyverse", "tidyr", "dplyr"))

library(tidyverse)
library(tidyr)
library(ggthemes)
library(RColorBrewer)
library(dplyr)

# Load IRMS dataset
irms_soil <- read.csv("./data_files/IRMS/IRMS_soil.csv",
                      sep = ",", header = TRUE,
                      stringsAsFactors = FALSE)

# Quick inspection
summary(irms_soil)
head(irms_soil)
glimpse(irms_soil)

# define numeric and factor cols
num_cols <- c("pct_n", "pct_c", "c_to_n", "delta_15n", "delta_13c", "LOI", "mass_6_8_mg")
fact_cols <- c("treatment", "layer", "beriget", "veg", "cut", "wet") 

# ---- Rename columns to convenient names ----
irms_soil_raw <- irms_soil %>%
  rename(
    run                 = run,
    treatment           = treatment,
    layer               = layer,
    LOI                 = LOI,
    tray                = tray,
    well                = well,
    nr                  = nr,
    mass_6_8_mg         = `X6.8mg`,           # Excel "4-5mg"
    beriget             = beriget,
    
    pct_n               = `X.N`,              # "%N"
    pct_c               = `X.C`,              # "%C"
    c_to_n              = `C.N`,              # "C/N"
    delta_15n           = d15Nkorr,
    delta_13c           = d13Ckorr
  ) %>%
  mutate(
    # vegetation type from treatment prefix
    veg = case_when(
      grepl("S", treatment) ~ "S",
      grepl("B", treatment) ~ "B",
      TRUE ~ NA_character_
    ),
    cut = case_when(
      grepl("U", treatment) ~ "U",
      grepl("C", treatment) ~ "C",
      TRUE ~ NA_character_
    ),
    wet = case_when(
      grepl("W", treatment) ~ "W",
      grepl("D", treatment) ~ "D",
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(layer != "top_poop") %>%  #remove poop samples 
  mutate(across(all_of(intersect(num_cols, names(.))), as.numeric)) %>% 
  mutate(across(all_of(intersect(fact_cols, names(.))), as.factor))






irms_soil_raw <- irms_soil_raw %>% 
  mutate(across(all_of(intersect(fact_cols, names(.))), as.factor))
 

# =================== correct atom pct ================================



irms_extracts<- read.csv("./data_files/IRMS/IRMS_extracts.csv",
                      sep = ",", header = TRUE,
                      stringsAsFactors = FALSE)

head(irms_extracts)


















#======================= plots =============================

long_soil_df <- irms_soil_raw %>%
  select(
    veg, beriget, layer,
    pct_n, pct_c, c_to_n, delta_15n, delta_13c
  ) %>%
  pivot_longer(
    cols = c(pct_n, pct_c, c_to_n, delta_15n, delta_13c),
    names_to = "metric_raw",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric_raw,
                    "pct_c" = "%C",
                    "pct_n" = "%N",
                    "c_to_n"  = "C/N",
                    "delta_15n"   = "δ15N (‰)",
                    "delta_13C"  = "δ13C (‰)"
    )
  )

se_fun <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

cmp_summ <- long_soil_df %>%
  group_by(metric, veg, layer, beriget) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    se   = se_fun(value),
    .groups = "drop"
  )

# before and after raw BELOW IS FROM ROOTS_GOOD.

p_soil_raw <- ggplot() +
  # raw values
  geom_jitter(
    data = long_soil_df,
    aes(x = beriget, y = value, color = beriget),
    width = 0.15, alpha = 0.5, size = 1.8
  ) +
  # mean ± SE
  geom_point(
    data = cmp_summ,
    aes(x = beriget, y = mean, fill = beriget),
    size = 3, shape = 21, color = "black"
  ) +
  geom_errorbar(
    data = cmp_summ,
    aes(x = beriget, ymin = mean - se, ymax = mean + se),
    width = 0.10, linewidth = 0.5
  ) +
  facet_grid(metric ~ veg + layer, scales = "free_y") +
  scale_color_manual(values = c("0" = "#7FB3D5", "1" = "#E59866")) +
  scale_fill_manual(values = c("0" = "#7FB3D5", "1" = "#E59866")) +
  labs(
    x = "labelled",
    y = "Value",
    title = "labelled / unlabbelled (raw data + mean ± SE)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA)
  )

p_soil_raw

# Optionally save
# ggsave("fig_irms_veg_beriget_raw_mean_se.png", p_soil_raw, width = 9, height = 7, dpi = 300)