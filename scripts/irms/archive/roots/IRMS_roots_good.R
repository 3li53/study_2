# ============================================================
# IRMS roots: before/after comparison (veg × diameter specific)
# Author: Elise

# ---- Packages ----
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)

# ---- 0) Read and pre-filter raw data ----
# NOTE: if your CSV has an extra header row, 'slice(-1)' removes it.
#       If you ever need to remove more lines, adjust accordingly.
irms_roots_raw <- read.csv("./data_files/IRMS/IRMS_roots_raw.csv", header = TRUE)

irms_d <- irms_roots_raw %>%
  filter(comment != "PEACH" | is.na(comment)) %>%  # drop PEACH rows
  slice(-1)                                         # remove duplicated header row (if present)

# ---- 1) Rename columns to convenient names ----
irms_da <- irms_d %>%
  rename(
    run                 = run,
    treatment           = treatment,
    layer               = layer,
    diameter            = diameter,
    tray                = tray,
    position            = position,
    nr                  = nr,
    mass_4_5_mg         = `X4.5mg`,           # Excel "4-5mg"
    comment             = comment,
    
    mg_n                = `mg.N.in.sample`,
    mg_c                = `mg.C.in.sample`,
    
    pct_n               = `X.N`,              # "%N"
    pct_c               = `X.C`,              # "%C"
    c_to_n              = `C.N`,              # "C/N"
    
    d15n_korr           = d15Nkorr,
    d13c_korr           = d13Ckorr,
    
    nat_abun_15n_atm_pct = `Nat.abun..15N.`,  # "Nat abun (15N)"
    delta15n_korr        = `δ15Nkorr`,
    atom_pct_15n         = `atom...15N.`,     # "atom% (15N)"
    ape_pct_15n          = `APE.`,            # "APE%"
    n_per_dw_mg_per_g    = `N.pr.DW`,         # "N pr DW"
    n15_per_dw_ug_per_g  = `X15N.pr.DW`,      # "15N pr DW"
    n15_per_n_ug_per_g   = `X15N.pr.N`,       # "15N pr N"
    
    nat_abun_13c_atm_pct = `Nat.abun..13C.`,  # "Nat abun (13C)"
    delta13c_korr        = `δ13Ckorr`,
    atom_pct_13c         = `atom...13C.`,     # "atom% (13C)"
    ape_pct_13c          = `APE..1`,          # second "APE%" for 13C side
    
    c_per_dw_mg_per_g    = `C.pr.DW`,         # "C pr DW"
    c13_per_dw_ug_per_g  = `X13C.pr.DW`,      # "13C pr DW"
    c13_per_c_ug_per_g   = `X13C.pr.C`        # "13C pr C"
  ) %>%
  mutate(
    # vegetation type from treatment prefix
    veg = case_when(
      grepl("S", treatment) ~ "S",
      grepl("B", treatment) ~ "B",
      TRUE ~ NA_character_
    )
  )

# ---- 2) Compute veg × diameter averages of natural abundance (run == 1) ----
avg_table <- irms_da %>%
  filter(run == 1) %>%
  group_by(veg, diameter) %>%
  summarise(
    avg_d15 = mean(d15n_korr, na.rm = TRUE),
    avg_d13 = mean(d13c_korr, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 3) Replace nat abundances with veg × diameter–specific values ----
#    Formulae:
#      atom%15N = 100 * R15 * ( (δ/1000 + 1) / (1 + R15 * (δ/1000 + 1)) ), where R15 = 0.003676
#      atom%13C = 100 * R13 * ( (δ/1000 + 1) / (1 + R13 * (δ/1000 + 1)) ), where R13 = 0.011237
irms_data <- irms_da %>%
  left_join(avg_table, by = c("veg", "diameter")) %>%
  mutate(
    nat_abun_15n_atm_pct = 100 * 0.003676 * ((avg_d15 / 1000) + 1) /
      (1 + 0.003676 * ((avg_d15 / 1000) + 1)),
    nat_abun_13c_atm_pct = 100 * 0.011237 * ((avg_d13 / 1000) + 1) /
      (1 + 0.011237 * ((avg_d13 / 1000) + 1))
  ) %>%
  select(-avg_d15, -avg_d13, -X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6, -tray, -position, -PEAnA, -PEA15N, - PEAnA.1, -PEA13C, -Name, -Sample.Number, -Height..nA., -Height..nA..1, -gnsnSTD.C.weight, -gnsnSTD.N.weight, -X1nA.to.x.mgC.STD, -X1nA.to.x.mgN.STD)  # keep veg/diameter; drop only helpers

# ---- 4) Ensure numeric types where needed ----
irms_data <- irms_data %>%
  mutate(
    across(
      c(atom_pct_15n, atom_pct_13c,
        nat_abun_15n_atm_pct, nat_abun_13c_atm_pct,
        pct_n, pct_c,
        n_per_dw_mg_per_g, c_per_dw_mg_per_g,
        mg_n, mg_c),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# ---- 5) Corrected APE (after) ----
irms_da4 <- irms_data %>%
  mutate(
    ape_pct_15n_correct = atom_pct_15n - nat_abun_15n_atm_pct,
    ape_pct_13c_correct = atom_pct_13c - nat_abun_13c_atm_pct
  )

# ---- 6) Build 'irms_full' with BEFORE and AFTER variables side-by-side ----
irms_full <- irms_da4 %>%
  mutate(
    # APE: before (original from file) vs after (corrected)
    ape_15n_before = as.numeric(ape_pct_15n),
    ape_13c_before = as.numeric(ape_pct_13c),
    ape_15n_after  = as.numeric(ape_pct_15n_correct),
    ape_13c_after  = as.numeric(ape_pct_13c_correct),
    
    # Quantities: BEFORE (using original APE)
    n15_dw_before = n_per_dw_mg_per_g * (ape_15n_before / 100) * 1000,
    n15_n_before  = (((mg_n * (ape_15n_before / 100)) / mg_n) * 1000) * 1000,
    c13_dw_before = c_per_dw_mg_per_g * (ape_13c_before / 100) * 1000,
    c13_c_before  = (((mg_c * (ape_13c_before / 100)) / mg_c) * 1000) * 1000,
    
    # Quantities: AFTER (using corrected APE)
    n15_dw_after = n_per_dw_mg_per_g * (ape_15n_after / 100) * 1000,
    n15_n_after  = (((mg_n * (ape_15n_after / 100)) / mg_n) * 1000) * 1000,
    c13_dw_after = c_per_dw_mg_per_g * (ape_13c_after / 100) * 1000,
    c13_c_after  = (((mg_c * (ape_13c_after / 100)) / mg_c) * 1000) * 1000
  ) %>%
  mutate(
    veg      = factor(veg, levels = c("S", "B")),
    diameter = fct_relevel(as.factor(diameter), "fine", "coarse")
  )

#library(writexl)
#write_xlsx(irms_full, "./data_files/IRMS/IRMS_roots_corrected.xlsx")

#----------------------------------------------------------
#---- 7) Tidy: long format for plotting ----

long_root_df <- irms_full %>%
  select(
    veg, diameter,
    ape_15n_before, ape_15n_after,
    ape_13c_before, ape_13c_after,
    n15_dw_before,  n15_dw_after,
    n15_n_before,   n15_n_after,
    c13_dw_before,  c13_dw_after,
    c13_c_before,   c13_c_after
  ) %>%
  pivot_longer(
    cols = -c(veg, diameter),
    names_to = c("metric_raw", "version"),
    names_pattern = "(.*)_(before|after)$",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric_raw,
                    "ape_15n" = "APE 15N (%)",
                    "ape_13c" = "APE 13C (%)",
                    "n15_dw"  = "µg 15N / g DW",
                    "n15_n"   = "µg 15N / g N",
                    "c13_dw"  = "µg 13C / g DW",
                    "c13_c"   = "µg 13C / g C"
    ),
    version = factor(version, levels = c("before", "after"))
  )

# ---- 8) Summaries: mean ± SE ----
se_fun <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

cmp_root_summ <- long_root_df %>%
  group_by(metric, veg, diameter, version) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    se   = se_fun(value),
    .groups = "drop"
  )

# ---- 9) Plots: before vs after with SE ----
# Option A — facet by metric (columns), rows = diameter, x = veg
p_by_diameter <- ggplot(cmp_root_summ, aes(x = veg, y = mean, fill = version)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
  geom_errorbar(
    aes(ymin = pmax(mean - se, NA_real_), ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.4
  ) +
  facet_grid(metric~diameter , scales = "free_y") +
  labs(
    x = "Vegetation type", y = "Mean (± SE)", fill = "Version",
    title = "Before vs After: veg × diameter–specific natural abundance corrections"
  ) +
  scale_fill_manual(values = c("before" = "#7FB3D5", "after" = "#E59866")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank()
  )

# Option B — facet by metric (columns), rows = veg, x = diameter
p_by_veg <- ggplot(cmp_root_summ, aes(x = diameter, y = mean, fill = version)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
  geom_errorbar(
    aes(ymin = pmax(mean - se, NA_real_), ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.4
  ) +
  facet_grid(metric~veg , scales = "free_y") +
  labs(
    x = "Root diameter class", y = "Mean (± SE)", fill = "Version",
    title = "Before vs After: veg × diameter–specific natural abundance corrections"
  ) +
  scale_fill_manual(values = c("before" = "#7FB3D5", "after" = "#E59866")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank()
  )

# Print
p_by_diameter
p_by_veg

# before and after raw

p_version_raw <- ggplot() +
  # raw values
  geom_jitter(
    data = long_root_df,
    aes(x = version, y = value, color = version),
    width = 0.15, alpha = 0.5, size = 1.8
  ) +
  # mean ± SE
  geom_point(
    data = cmp_root_summ,
    aes(x = version, y = mean, fill = version),
    size = 3, shape = 21, color = "black"
  ) +
  geom_errorbar(
    data = cmp_root_summ,
    aes(x = version, ymin = mean - se, ymax = mean + se),
    width = 0.10, linewidth = 0.5
  ) +
  facet_grid(metric ~ veg + diameter, scales = "free_y") +
  scale_color_manual(values = c("before" = "#7FB3D5", "after" = "#E59866")) +
  scale_fill_manual(values = c("before" = "#7FB3D5", "after" = "#E59866")) +
  labs(
    x = "version",
    y = "Value",
    title = "before vs after (raw data + mean ± SE)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA)
  )

p_version_raw

# labelled vs unlabelled

irms_full2 <- irms_full %>%
mutate(
  label_status = case_when(
    run == 1 ~ "unlabelled",
    run %in% 2:6 ~ "labelled",
    TRUE ~ NA_character_
  ),
  label_status = factor(label_status, levels = c("unlabelled", "labelled"))
)

# long format

long_label <- irms_full2 %>%
  select(
    veg, diameter, label_status,
    n15_dw_after, n15_n_after,
    c13_dw_after, c13_c_after
  ) %>%
  pivot_longer(
    cols = -c(veg, diameter, label_status),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    "n15_dw_after" = "µg 15N / g DW",
                    "n15_n_after"  = "µg 15N / g N",
                    "c13_dw_after" = "µg 13C / g DW",
                    "c13_c_after"  = "µg 13C / g C"
    )
  )

# summary ( mean +- SE)

label_summ <- long_label %>%
  group_by(metric, veg, diameter, label_status) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# plot

# raw

p_label_raw <- ggplot() +
  # raw values
  geom_jitter(
    data = long_label,
    aes(x = label_status, y = value, color = label_status),
    width = 0.15, alpha = 0.5, size = 2
  ) +
  # mean ± SE
  geom_point(
    data = label_summ,
    aes(x = label_status, y = mean, fill = label_status),
    size = 3, shape = 21, color = "black"
  ) +
  geom_errorbar(
    data = label_summ,
    aes(x = label_status, ymin = mean - se, ymax = mean + se),
    width = 0.20, linewidth = 0.5
  ) +
  facet_grid(metric ~ veg + diameter, scales = "free_y") +
  scale_color_manual(values = c("unlabelled" = "#7FB3D5", "labelled" = "#E59866")) +
  scale_fill_manual(values = c("unlabelled" = "#7FB3D5", "labelled" = "#E59866")) +
  labs(
    x = "Label status",
    y = "Value",
    title = "Labelled vs Unlabelled (raw data + mean ± SE)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    strip.background = element_rect(fill = "grey92", color = NA)
  )

p_label_raw

# point

p_label_points <- ggplot(label_summ,
                         aes(x = label_status, y = mean, color = label_status)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.15, linewidth = 0.5,
    position = position_dodge(width = 0.5)
  ) +
  facet_grid(metric ~ veg + diameter, scales = "free_y") +
  scale_color_manual(values = c("unlabelled" = "#7FB3D5", "labelled" = "#E59866")) +
  labs(
    x = "Label status",
    y = "Mean (± SE)",
    title = "Labelled vs Unlabelled (points only, no bars)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA)
  )

p_label_points

# bar

p_label_bar <- ggplot(label_summ,
                  aes(x = label_status, y = mean, fill = label_status)) +
  geom_col(width = 0.7, color = "grey30") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.25, linewidth = 0.4
  ) +
  facet_grid(metric ~ veg + diameter, scales = "free_y") +
  scale_fill_manual(values = c("unlabelled" = "#7FB3D5", "labelled" = "#E59866")) +
  labs(
    x = "Label status",
    y = "Mean (± SE)",
    title = "Labelled vs Unlabelled: 15N and 13C metrics"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank()
  )

p_label_bar


# treatments


