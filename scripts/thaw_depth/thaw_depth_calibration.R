# ---- thaw depth linear calibration ---- copilot conv. "what is it called when you take some data...

# ---- Minimal Thaw Depth Calibration & Append ----

# 0) Packages
pkgs <- c("dplyr", "readr", "ggplot2", "purrr", "broom", "tidyr")
invisible(lapply(pkgs, function(p) if (!require(p, character.only = TRUE)) {
  install.packages(p, dependencies = TRUE); library(p, character.only = TRUE)
}))

# 1) Read & clean data
data_dir <- "./data_files/thaw_depth/"
file_pattern <- "study_2_thaw_depth_run[1-6]\\.csv"
files <- list.files(data_dir, pattern = file_pattern, full.names = TRUE)

# Treat "-" and "" as NA at import
all_thaw_depth <- files |>
  purrr::map(~ read.csv(.x, stringsAsFactors = FALSE, na.strings = c("", "-"))) |>
  dplyr::bind_rows()

# Ensure numeric columns
df <- all_thaw_depth |>
  dplyr::mutate(
    thaw_depth_cm        = as.numeric(thaw_depth_cm),
    thaw_depth_center    = as.numeric(thaw_depth_center),
    thaw_depth_periphery = as.numeric(thaw_depth_periphery),
    hrs_after_incubation = as.numeric(hrs_after_incubation)
  )

# 2) Fit global calibrations (A ~ B) using overlap rows
min_points <- 3
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

overlap_center <- df |> dplyr::filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_center))
overlap_periph <- df |> dplyr::filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_periphery))

fit_center <- if (nrow(overlap_center) >= min_points)
  lm(thaw_depth_cm ~ thaw_depth_center, data = overlap_center) else NULL

fit_periph <- if (nrow(overlap_periph) >= min_points)
  lm(thaw_depth_cm ~ thaw_depth_periphery, data = overlap_periph) else NULL

# Weights from RMSE (down‑weight noisier predictor)
w_center <- if (!is.null(fit_center)) {
  r <- rmse(overlap_center$thaw_depth_cm, predict(fit_center))
  if (is.finite(r) && r > 0) 1 / r^2 else 0
} else 0

w_periph <- if (!is.null(fit_periph)) {
  r <- rmse(overlap_periph$thaw_depth_cm, predict(fit_periph))
  if (is.finite(r) && r > 0) 1 / r^2 else 0
} else 0

# 3) Predict A‑scale from Center/Periphery, combine, and append to original df
df_aug <- df |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  dplyr::mutate(
    pred_from_center = if (!is.null(fit_center) & !is.na(thaw_depth_center)) {
      as.numeric(predict(fit_center, newdata = data.frame(thaw_depth_center = thaw_depth_center)))
    } else NA_real_,
    pred_from_periph = if (!is.null(fit_periph) & !is.na(thaw_depth_periphery)) {
      as.numeric(predict(fit_periph, newdata = data.frame(thaw_depth_periphery = thaw_depth_periphery)))
    } else NA_real_,
    pred_from_B = dplyr::case_when(
      !is.na(pred_from_center) & !is.na(pred_from_periph) & (w_center + w_periph) > 0 ~
        (w_center * pred_from_center + w_periph * pred_from_periph) / (w_center + w_periph),
      !is.na(pred_from_center) ~ pred_from_center,
      !is.na(pred_from_periph) ~ pred_from_periph,
      TRUE ~ NA_real_
    ),
    thaw_depth_harmonized = dplyr::coalesce(thaw_depth_cm, pred_from_B),
    source = dplyr::case_when(
      !is.na(thaw_depth_cm) ~ "Method A (raw)",
      is.na(thaw_depth_cm) & !is.na(pred_from_B) ~ "Method B (calibrated)",
      TRUE ~ "Missing"
    ),
    resid_center = dplyr::if_else(!is.na(thaw_depth_cm) & !is.na(pred_from_center),
                                  thaw_depth_cm - pred_from_center, NA_real_),
    resid_periph = dplyr::if_else(!is.na(thaw_depth_cm) & !is.na(pred_from_periph),
                                  thaw_depth_cm - pred_from_periph, NA_real_)
  )

# 4) Quick sanity plot & save (works even if pane is moody)
options(device = "RStudioGD")
p_quick <- ggplot2::ggplot(df_aug, ggplot2::aes(hrs_after_incubation, thaw_depth_harmonized, color = treatment)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::geom_line(ggplot2::aes(group = interaction(treatment, pipe_nr)), alpha = 0.5) +
  ggplot2::labs(title = "Harmonized thaw depth (global calibration)",
                x = "Hours after incubation", y = "Thaw depth (A-scale, cm)") +
  ggplot2::theme_bw()

print(p_quick)
dir.create("./figures", showWarnings = FALSE)
ggplot2::ggsave("./figures/harmonized_quick.png", p_quick, width = 10, height = 7, dpi = 300)

# 5) Export augmented data if you want to archive it
# write.csv(df_aug, "./data_files/thaw_depth/thaw_depth_augmented.csv", row.names = FALSE)

# Done: df_aug now contains pred_from_center, pred_from_periph, pred_from_B, thaw_depth_harmonized, source, and residuals.
































# 1) Load required packages
packages <- c("dplyr", "ggplot2", "broom", "readr", "tidyr", "purrr")
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# 2) Read and combine thaw depth files
data_dir <- "./data_files/thaw_depth/"
file_pattern <- "study_2_thaw_depth_run[1-6]\\.csv"
files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)

# Treat "-" and "" as NA during import
thaw_depth_list <- lapply(files, function(f) {
  read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "-"))
})
df <- bind_rows(thaw_depth_list)

# 3) Clean numeric columns
df <- df %>%
  mutate(
    thaw_depth_cm        = as.numeric(thaw_depth_cm),
    thaw_depth_center    = as.numeric(thaw_depth_center),
    thaw_depth_periphery = as.numeric(thaw_depth_periphery),
    hrs_after_incubation = as.numeric(hrs_after_incubation)
  )

# Quick check
str(df[, c("thaw_depth_cm", "thaw_depth_center", "thaw_depth_periphery")])

# 4) Define helper for RMSE
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

# 5) Global fallback models
min_points <- 3
global_overlap_center <- df %>% filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_center))
global_overlap_periph <- df %>% filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_periphery))

global_fit_center <- if (nrow(global_overlap_center) >= min_points)
  lm(thaw_depth_cm ~ thaw_depth_center, data = global_overlap_center) else NULL
global_fit_periph <- if (nrow(global_overlap_periph) >= min_points)
  lm(thaw_depth_cm ~ thaw_depth_periphery, data = global_overlap_periph) else NULL

global_rmse_center <- if (!is.null(global_fit_center)) rmse(global_overlap_center$thaw_depth_cm, predict(global_fit_center)) else NA_real_
global_rmse_periph <- if (!is.null(global_fit_periph)) rmse(global_overlap_periph$thaw_depth_cm, predict(global_fit_periph)) else NA_real_
w_center_global <- if (is.finite(global_rmse_center) && global_rmse_center > 0) 1 / global_rmse_center^2 else 0
w_periph_global <- if (is.finite(global_rmse_periph) && global_rmse_periph > 0) 1 / global_rmse_periph^2 else 0

# 6) Per-treatment fits
center_fits <- df %>%
  filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_center)) %>%
  group_by(treatment) %>%
  nest() %>%
  mutate(
    n   = map_int(data, nrow),
    fit  = map2(data, n, ~ if (.y >= min_points) lm(thaw_depth_cm ~ thaw_depth_center, data = .x) else NULL),
    rmse = map2_dbl(data, fit, ~ if (!is.null(.y)) rmse(.x$thaw_depth_cm, predict(.y)) else NA_real_),
    w    = if_else(is.finite(rmse) & rmse > 0, 1 / rmse^2, 0)
  ) %>%
  select(treatment, fit, rmse, w)

periph_fits <- df %>%
  filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_periphery)) %>%
  group_by(treatment) %>%
  nest() %>%
  mutate(
    n    = map_int(data, nrow),
    fit  = map2(data, n, ~ if (.y >= min_points) lm(thaw_depth_cm ~ thaw_depth_periphery, data = .x) else NULL),
    rmse = map2_dbl(data, fit, ~ if (!is.null(.y)) rmse(.x$thaw_depth_cm, predict(.y)) else NA_real_),
    w    = if_else(is.finite(rmse) & rmse > 0, 1 / rmse^2, 0)
  ) %>%
  select(treatment, fit, rmse, w)

# 7) Predict and harmonize
df_pred <- df %>%
  left_join(center_fits, by = "treatment", suffix = c("", "_center")) %>%
  rename(fit_center = fit, rmse_center = rmse, w_center = w) %>%
  left_join(periph_fits, by = "treatment", suffix = c("", "_periph")) %>%
  rename(fit_periph = fit, rmse_periph = rmse, w_periph = w)

df_h <- df_pred %>%
  rowwise() %>%
  mutate(
    pred_from_center = if (!is.na(thaw_depth_center)) {
      if (!is.null(fit_center)) as.numeric(predict(fit_center, newdata = data.frame(thaw_depth_center = thaw_depth_center)))
      else if (!is.null(global_fit_center)) as.numeric(predict(global_fit_center, newdata = data.frame(thaw_depth_center = thaw_depth_center)))
      else NA_real_
    } else NA_real_,
    pred_from_periph = if (!is.na(thaw_depth_periphery)) {
      if (!is.null(fit_periph)) as.numeric(predict(fit_periph, newdata = data.frame(thaw_depth_periphery = thaw_depth_periphery)))
      else if (!is.null(global_fit_periph)) as.numeric(predict(global_fit_periph, newdata = data.frame(thaw_depth_periphery = thaw_depth_periphery)))
      else NA_real_
    } else NA_real_,
    w_c = if (!is.null(fit_center) && is.finite(w_center) && w_center > 0) w_center else w_center_global,
    w_p = if (!is.null(fit_periph) && is.finite(w_periph) && w_periph > 0) w_periph else w_periph_global,
    pred_from_B = case_when(
      !is.na(pred_from_center) & !is.na(pred_from_periph) & (w_c + w_p) > 0 ~ (w_c * pred_from_center + w_p * pred_from_periph) / (w_c + w_p),
      !is.na(pred_from_center) ~ pred_from_center,
      !is.na(pred_from_periph) ~ pred_from_periph,
      TRUE ~ NA_real_
    ),
    thaw_depth_harmonized = coalesce(thaw_depth_cm, pred_from_B),
    source = case_when(
      !is.na(thaw_depth_cm) ~ "Method A (raw)",
      is.na(thaw_depth_cm) & !is.na(pred_from_B) ~ "Method B (calibrated)",
      TRUE ~ "Missing"
    )
  ) %>%
  ungroup()

# 8) Visualization
options(device = "RStudioGD") # Ensure plots show in RStudio
dir.create("./figures", showWarnings = FALSE)

# Calibration plots
p_center <- df %>%
  filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_center)) %>%
  ggplot(aes(thaw_depth_center, thaw_depth_cm, color = treatment)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ treatment, scales = "free") +
  labs(title = "Calibration by treatment: Center → A-scale",
       x = "Center (Method B, cm)", y = "Method A (cm)") +
  theme_bw()

p_periph <- df %>%
  filter(!is.na(thaw_depth_cm), !is.na(thaw_depth_periphery)) %>%
  ggplot(aes(thaw_depth_periphery, thaw_depth_cm, color = treatment)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ treatment, scales = "free") +
  labs(title = "Calibration by treatment: Periphery → A-scale",
       x = "Periphery (Method B, cm)", y = "Method A (cm)") +
  theme_bw()

p_ts <- df_h %>%
  ggplot(aes(hrs_after_incubation, thaw_depth_harmonized, group = pipe_nr)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(shape = source, color = source), size = 1.8) +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_color_manual(values = c("Method A (raw)" = "#1b9e77",
                                "Method B (calibrated)" = "#d95f02",
                                "Missing" = "grey70")) +
  scale_shape_manual(values = c("Method A (raw)" = 16,
                                "Method B (calibrated)" = 17,
                                "Missing" = 4)) +
  labs(title = "Harmonized thaw depth by treatment",
       x = "Hours after incubation", y = "Thaw depth (A-scale, cm)") +
  theme_bw()

# Print plots
print(p_center)
print(p_periph)
print(p_ts)

# Save plots
ggsave("./figures/calibration_center.png", p_center, width = 9, height = 7, dpi = 300)
ggsave("./figures/calibration_periphery.png", p_periph, width = 9, height = 7, dpi = 300)
ggsave("./figures/harmonized_timeseries.png", p_ts, width = 10, height = 7, dpi = 300)





# OPTION A Make sure df_h exists from the harmonization code
# df_h contains pred_from_center, pred_from_periph, thaw_depth_harmonized, source

# Add a stable row_id on both sides to align row-wise
df_aug <- df %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::left_join(
    df_h %>%
      dplyr::mutate(row_id = dplyr::row_number()) %>%
      dplyr::select(
        row_id,
        pred_from_center, pred_from_periph, pred_from_B,
        thaw_depth_harmonized, source
      ),
    by = "row_id"
  ) %>%
  # Add residuals on overlap (only where A and the corresponding B exist)
  dplyr::mutate(
    resid_center  = dplyr::if_else(!is.na(thaw_depth_cm) & !is.na(thaw_depth_center),
                                   thaw_depth_cm - pred_from_center, NA_real_),
    resid_periph  = dplyr::if_else(!is.na(thaw_depth_cm) & !is.na(thaw_depth_periphery),
                                   thaw_depth_cm - pred_from_periph, NA_real_)
  )

# Inspect
str(df_aug[, c("pred_from_center","pred_from_periph","pred_from_B",
               "thaw_depth_harmonized","source","resid_center","resid_periph")])

#OPTION B
all_thaw_depth_aug <- all_thaw_depth %>%
  dplyr::mutate(row_id = dplyr::row_number())

df_h_small <- df_h %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::select(
    row_id,
    pred_from_center, pred_from_periph, pred_from_B,
    thaw_depth_harmonized, source
  )

all_thaw_depth_aug <- all_thaw_depth_aug %>%
  dplyr::left_join(df_h_small, by = "row_id")

# Check new columns are present
names(all_thaw_depth_aug)
