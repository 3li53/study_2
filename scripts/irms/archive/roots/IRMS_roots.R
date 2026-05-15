#script is used to try things out. Use IRMS_roots_good for a cleaner script.


#pipes_data_soil <- read.csv("./data_files/pipes_data/roots/pipes_data_roots.csv", header = TRUE, sep = ",", skip = 2)

#IRMS_roots <- data.frame(
#  nr = c(1:288),
#  run = rep(1:6, each = 48),
#  treatment = rep(c("BCD", "BCW", "BUD", "BUW", "SCD", "SCW", "SUD", "SUW"), each = , times = 6),
#  layer = rep(c("top", "middle", "bottom"), each = 2, times = 48),
#  diameter = rep(c("fine", "coarse"), each = 1,  times = 144)
#)

# library(writexl)
# write_xlsx(IRMS_roots, "./data_files/IRMS/IRMS_roots.xlsx")

library(dplyr)

irms_raw <- read.csv("./data_files/IRMS/IRMS_roots_raw.csv", header = TRUE)

irms_d <- irms_raw %>%
  filter(comment != "PEACH" | is.na(comment)) %>% 
  slice(-1) # remove second header (-1), and first four rows (-2:-5)

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
    delta15n_korr        = `δ15Nkorr`,        # "δ15Nkorr" (text column in CSV)
    atom_pct_15n         = `atom...15N.`,     # "atom% (15N)"
    ape_pct_15n          = `APE.`,            # "APE%"
    n_per_dw_mg_per_g    = `N.pr.DW`,         # "N pr DW"
    n15_per_dw_ug_per_g  = `X15N.pr.DW`,      # "15N pr DW"
    n15_per_n_ug_per_g   = `X15N.pr.N`,       # "15N pr N"
    
    nat_abun_13c_atm_pct = `Nat.abun..13C.`,  # "Nat abun (13C)"
    delta13c_korr        = `δ13Ckorr`,        # "δ13Ckorr" (text column in CSV)
    atom_pct_13c         = `atom...13C.`,     # "atom% (13C)"
    ape_pct_13c          = `APE..1`,          # second "APE%" for 13C side
    
    c_per_dw_mg_per_g    = `C.pr.DW`,         # "C pr DW"
    c13_per_dw_ug_per_g  = `X13C.pr.DW`,      # "13C pr DW"
    c13_per_c_ug_per_g   = `X13C.pr.C`        # "13C pr C"
  )

#new column for vegetation types
irms_dat <- irms_da %>%
  mutate(veg = ifelse(grepl("^S", treatment), "S",
                      ifelse(grepl("^B", treatment), "B", NA)))

# Compute averages for natural abundance rows (run == 1)
avg_table <- irms_dat %>%
  filter(run == 1) %>%
  group_by(veg, diameter) %>%
  summarise(
    avg_d15 = mean(d15n_korr, na.rm = TRUE),
    avg_d13 = mean(d13c_korr, na.rm = TRUE),
    .groups = "drop"
  )

avg_table

#replace average column with vegetation and diameter specific averages
irms_data <- irms_dat %>%
  left_join(avg_table, by = c("veg", "diameter")) %>%
  mutate(
    # replace nat abundance for N based on veg × diameter
    nat_abun_15n_atm_pct = 100 * 0.003676 * ((avg_d15 / 1000) + 1)/(1 + 0.003676*(avg_d15 / 1000 + 1)), # = 100 * 0.003676*((avg_d15/1000)+1)/(1+0.003676*(avg_d15/1000+1))
    
    # replace nat abundance for C
    nat_abun_13c_atm_pct = 100 * 0.011237 * ((avg_d13 / 1000) + 1)/(1 + 0.011237*(avg_d13 / 1000 + 1))
  ) %>%
  select(-avg_d15, -avg_d13, -veg, -X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6)   # clean up helper columns



# APE 

# first, change data type from chr to num
irms_data <- irms_data %>%
  mutate(
    atom_pct_15n  = as.numeric(atom_pct_15n),
    atom_pct_13c  = as.numeric(atom_pct_13c),
    nat_abun_15n_atm_pct = as.numeric(nat_abun_15n_atm_pct),
    nat_abun_13c_atm_pct = as.numeric(nat_abun_13c_atm_pct)
  )


# calculate correct APE for 15N and 13C
irms_da4 <- irms_data %>%
  mutate(
    ape_pct_15n_correct = atom_pct_15n - nat_abun_15n_atm_pct,
    ape_pct_13c_correct = atom_pct_13c - nat_abun_13c_atm_pct
  )


# correct the datatypes

irms_da4 <- irms_da4 %>%
  mutate(
    pct_n = as.numeric(pct_n),
    pct_c = as.numeric(pct_c),
    
    n_per_dw_mg_per_g = as.numeric(n_per_dw_mg_per_g),
    c_per_dw_mg_per_g = as.numeric(c_per_dw_mg_per_g),
    
    ape_pct_15n = as.numeric(ape_pct_15n),
    ape_pct_13c = as.numeric(ape_pct_13c),
    
    ape_pct_15n_correct = as.numeric(ape_pct_15n_correct),
    ape_pct_13c_correct = as.numeric(ape_pct_13c_correct)
  )


# correct isotope pr DW and isotope pr atom 

irms_full <- irms_da4 %>%
  mutate(
    
    # -----------------------
    # BEFORE (original APE)
    # -----------------------
    ape_15n_before = ape_pct_15n,
    ape_13c_before = ape_pct_13c,
    
    n15_dw_before = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000,
    n15_n_before  = (((mg_n * (ape_pct_15n / 100)) / mg_n) * 1000) * 1000,
    
    c13_dw_before = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,
    c13_c_before  = (((mg_c * (ape_pct_13c / 100)) / mg_c) * 1000) * 1000,
    
    # -----------------------
    # AFTER (corrected APE)
    # -----------------------
    ape_15n_after = ape_pct_15n_correct,
    ape_13c_after = ape_pct_13c_correct,
    
    n15_dw_after = n_per_dw_mg_per_g * (ape_pct_15n_correct / 100) * 1000,
    n15_n_after  = (((mg_n * (ape_pct_15n_correct / 100)) / mg_n) * 1000) * 1000,
    
    c13_dw_after = c_per_dw_mg_per_g * (ape_pct_13c_correct / 100) * 1000,
    c13_c_after  = (((mg_c * (ape_pct_13c_correct / 100)) / mg_c) * 1000) * 1000
  )

# ignore below--------------------

irms_da5 <- irms_da4 %>%
  mutate(
    
    # 15N
    n15_per_dw_ug_per_g_specific = n_per_dw_mg_per_g * (ape_pct_15n_correct / 100) * 1000, # µg 15N per g DW =(AQ7*AP7%)*1000
    n15_per_n_ug_per_g_specific  = ((((mg_n * (ape_pct_15n_correct / 100)) / mg_n) * 1000) * 1000),           # µg 15N per g N =(((AF7*AP7%)/AF7)*1000)*1000
    
    # 13C
    c13_per_dw_ug_per_g_specific = c_per_dw_mg_per_g * (ape_pct_13c_correct / 100) * 1000,      # µg 13C per g DW
    c13_per_c_ug_per_g_specific  = ((((mg_c * (ape_pct_13c_correct / 100)) / mg_c) * 1000) * 1000)              # µg 13C per g C
  )

#original

irms_da5 <- irms_da4 %>%
  mutate(
    
    # 15N
    n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000, # µg 15N per g DW =(AQ7*AP7%)*1000
    n15_per_n_ug_per_g  = ((((mg_n * (ape_pct_15n / 100)) / mg_n) * 1000) * 1000),           # µg 15N per g N =(((AF7*AP7%)/AF7)*1000)*1000
    
    # 13C
    c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,      # µg 13C per g DW
    c13_per_c_ug_per_g  = ((((mg_c * (ape_pct_13c / 100)) / mg_c) * 1000) * 1000)              # µg 13C per g C
  )


# see the difference ------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Ensure veg is available (if you don't already have it here)
# (You created it earlier; re-create here for safety if needed)
add_veg <- function(df) {
  if (!"veg" %in% names(df)) {
    df <- df %>%
      mutate(veg = ifelse(grepl("^S", treatment), "S",
                          ifelse(grepl("^B", treatment), "B", NA)))
  }
  df
}

irms_da4 <- add_veg(irms_full)
irms_da5 <- add_veg(irms_da5)

# Pick which metrics you want to compare (you can edit this vector)
# Options present from your objects:
# - APE:      ape_pct_15n vs ape_pct_15n_correct
#            ape_pct_13c vs ape_pct_13c_correct
# - Quant.:   n15_per_dw_ug_per_g vs n15_per_dw_ug_per_g_specific
#            n15_per_n_ug_per_g  vs n15_per_n_ug_per_g_specific
#            c13_per_dw_ug_per_g vs c13_per_dw_ug_per_g_specific
#            c13_per_c_ug_per_g  vs c13_per_c_ug_per_g_specific
metrics_to_compare <- list(
  "APE 15N (%)"           = c(before = "ape_pct_15n",                after = "ape_pct_15n_correct"),
  "APE 13C (%)"           = c(before = "ape_pct_13c",                after = "ape_pct_13c_correct"),
  "µg 15N / g DW"         = c(before = "n15_per_dw_ug_per_g",        after = "n15_per_dw_ug_per_g_specific"),
  "µg 15N / g N"          = c(before = "n15_per_n_ug_per_g",         after = "n15_per_n_ug_per_g_specific"),
  "µg 13C / g DW"         = c(before = "c13_per_dw_ug_per_g",        after = "c13_per_dw_ug_per_g_specific"),
  "µg 13C / g C"          = c(before = "c13_per_c_ug_per_g",         after = "c13_per_c_ug_per_g_specific")
)

# Helper to safely pull columns if they exist in a data frame
pull_if_exists <- function(df, col) if (col %in% names(df)) df[[col]] else NA_real_

# Build a long table with veg, diameter, and metric values for before/after
# We’ll take "before" from irms_da5 (original) and "after" from irms_da4 (corrected APE & derived)
cmp_long <- purrr::imap_dfr(
  metrics_to_compare,
  function(cols, metric_nm) {
    tibble(
      veg      = irms_full$veg,
      diameter = irms_full$diameter,
      before   = pull_if_exists(irms_full, cols["before"]),
      after    = pull_if_exists(irms_full, cols["after"])
    ) %>%
      pivot_longer(cols = c(before, after), names_to = "version", values_to = "value") %>%
      mutate(metric = metric_nm)
  }
) %>%
  # optional: keep only rows with both veg & diameter present
  filter(!is.na(veg), !is.na(diameter))

# (Optional) Clean diameter ordering if needed
cmp_long <- cmp_long %>%
  mutate(
    veg      = factor(veg, levels = c("S", "B")),
    diameter = fct_relevel(as.factor(diameter), "fine", "coarse"),
    version  = factor(version, levels = c("before", "after"))
  )


se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

cmp_summ <- cmp_long %>%
  group_by(metric, veg, diameter, version) %>%
  summarise(
    n      = sum(!is.na(value)),
    mean   = mean(value, na.rm = TRUE),
    se     = se(value),
    .groups = "drop"
  )


p_a <- ggplot(cmp_summ, aes(x = veg, y = mean, fill = version)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.4
  ) +
  facet_grid(diameter ~ metric, scales = "free") +
  labs(
    x = "Vegetation type",
    y = "Mean (± SE)",
    fill = "Version",
    title = "Before vs After: Effect of veg×diameter-specific natural abundance corrections"
  ) +
  scale_fill_manual(values = c("before" = "#7FB3D5", "after" = "#E59866")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank()
  )

p_a


p_b <- ggplot(cmp_summ, aes(x = diameter, y = mean, fill = version)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.4
  ) +
  facet_grid(veg ~ metric, scales = "free") +
  labs(
    x = "Root diameter class",
    y = "Mean (± SE)",
    fill = "Version",
    title = "Before vs After: Effect of veg×diameter-specific natural abundance corrections"
  ) +
  scale_fill_manual(values = c("before" = "#7FB3D5", "after" = "#E59866")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank()
  )

p_b
