
# IRMS ROOTS CLEAN PIPELINE

library(tidyverse)
library(writexl)

# Import

irms_raw <- read.csv("./data_files/IRMS/IRMS_roots_raw.csv", header = TRUE)
glimpse(irms_raw)

# CLEANUP

cols_to_remove <- c(
  "tray", "position", "X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6",
  "Sample.Number", "Name", "Height..nA.", "N15", "Height..nA..1",
  "X13C", "Weight", "PEAnA", "PEA15N", "PEAnA.1", "PEA13C",
  "gnsnSTD.C.weight", "gnsnSTD.N.weight",
  "X1nA.to.x.mgC.STD", "X1nA.to.x.mgN.STD"
)

numeric_cols <- c(
  "nr", "X4.5mg",
  "mg.N.in.sample",   "mg.C.in.sample",
  "X.N", "X.C", "C.N",
  "d15Nkorr", "d13Ckorr",
  "Nat.abun..15N.", "δ15Nkorr", "atom...15N.", "APE.",
  "N.pr.DW", "X15N.pr.DW", "X15N.pr.N",
  "Nat.abun..13C.", "δ13Ckorr", "atom...13C.", "APE..1",
  "C.pr.DW", "X13C.pr.DW", "X13C.pr.C"
)

factor_cols <- c("run", "treatment", "layer", "diameter")

# ---- CLEAN DATA ----

irms_clean <- irms_raw %>%
  filter(!comment %in% c("PEACH", "too small") | is.na(comment)) %>%
  slice(-1) %>%             # safe and simple
  mutate(
    across(any_of(factor_cols), as.factor),
    across(any_of(numeric_cols), ~ suppressWarnings(as.numeric(.))),
    comment = as.character(comment)
  ) %>%
  select(-any_of(cols_to_remove))

glimpse(irms_clean)

# RENAME COLUMNS

irms <- irms_clean %>%
  rename(
    mass_4_5_mg            = X4.5mg,
    mg_n                   = mg.N.in.sample,
    mg_c                   = mg.C.in.sample,
    pct_n                  = X.N,
    pct_c                  = X.C,
    c_to_n                 = C.N,
    d15n_korr              = d15Nkorr,
    d13c_korr              = d13Ckorr,
    nat_abun_15n_atm_pct   = Nat.abun..15N.,
    delta15n_korr          = δ15Nkorr,
    atom_pct_15n           = atom...15N.,
    ape_pct_15n            = APE.,
    n_per_dw_mg_per_g      = N.pr.DW,
    n15_per_dw_ug_per_g    = X15N.pr.DW,
    n15_per_n_ug_per_g     = X15N.pr.N,
    nat_abun_13c_atm_pct   = Nat.abun..13C.,
    delta13c_korr          = δ13Ckorr,
    atom_pct_13c           = atom...13C.,
    ape_pct_13c            = APE..1,
    c_per_dw_mg_per_g      = C.pr.DW,
    c13_per_dw_ug_per_g    = X13C.pr.DW,
    c13_per_c_ug_per_g     = X13C.pr.C
  )
glimpse(irms)

# TREATMENT PARSING + BERICHT

irms <- irms %>%
  mutate(
    veg = str_extract(treatment, "[SB]"),
    cut = str_extract(treatment, "[UC]"),
    wet = str_extract(treatment, "[WD]"),
    beriget = run != 1
  )

# OUTLIER DETECTION FUNCTION

iqr_outlier <- function(df, group_var, numeric_vars) {
  df %>%
    group_by(across(all_of(group_var))) %>%
    mutate(across(all_of(numeric_vars),
                  ~ (. < quantile(., 0.25, na.rm = TRUE) - 1.5 * IQR(., na.rm = TRUE)) |
                    (. > quantile(., 0.75, na.rm = TRUE) + 1.5 * IQR(., na.rm = TRUE)),
                  .names = "outlier_{.col}"
    )) %>%
    ungroup()
}

# NATURAL ABUNDANCE BASELINE (RUN 1)

natabun <- irms %>% filter(!beriget)

natabun_means <- natabun %>%
  group_by(diameter, layer) %>%
  summarise(
    avg_d15 = mean(d15n_korr, na.rm = TRUE),
    avg_d13 = mean(d13c_korr, na.rm = TRUE),
    .groups = "drop"
  )

# ISOTOPE FORMULA HELPERS

R15 <- 0.003676   # air N2
R13 <- 0.01118    # VPDB

atom_pct <- function(delta, Rstd) {
  Rsample <- Rstd * (1 + delta / 1000)
  100 * (Rsample / (1 + Rsample))
}

# APPLY BASELINE CORRECTIONS

irms <- irms %>%
  left_join(natabun_means, by = c("diameter", "layer")) %>%
  mutate(
    nat_abun_15n_atm_pct = atom_pct(avg_d15, R15),
    nat_abun_13c_atm_pct = atom_pct(avg_d13, R13),
    
    atom_pct_15n = atom_pct(d15n_korr, R15),
    atom_pct_13c = atom_pct(d13c_korr, R13),
    
    ape_pct_15n = atom_pct_15n - nat_abun_15n_atm_pct,
    ape_pct_13c = atom_pct_13c - nat_abun_13c_atm_pct,
    
    n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000,
    n15_per_n_ug_per_g  = mg_n              * (ape_pct_15n / 100) * 1000,
    
    c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,
    c13_per_c_ug_per_g  = mg_c              * (ape_pct_13c / 100) * 1000
  )

# OUTLIERS ON CORRECTED DATA

outlier_vars <- c(
  "d15n_korr", "d13c_korr",
  "n15_per_dw_ug_per_g", "c13_per_dw_ug_per_g",
  "n15_per_n_ug_per_g", "c13_per_c_ug_per_g"
)

irms_outliers <- iqr_outlier(irms, group_var = c("diameter", "layer"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

# PLOT THE DATA

# boxplot

plot_box <- function(df, x, y, colour = "beriget", title = NULL) {
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    theme_minimal(base_size = 16) +
    labs(title = title, x = x, y = y)
}

plot_box(irms, "layer", "n15_per_dw_ug_per_g")
plot_box(irms, "layer", "c13_per_dw_ug_per_g")
plot_box(irms, "diameter", "d13c_korr")
plot_box(irms, "veg", "n15_per_dw_ug_per_g")

# scatterplot

plot_scatter <- function(df, x = "d13c_korr", y = "d15n_korr",
                         colour = "diameter", title = NULL) {
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 16) +
    labs(title = title, x = x, y = y)
}

plot_scatter(irms, title = "δ13C vs δ15N")

# with outliers 

plot_scatter_outliers <- function(df, x, y, colour = "diameter") {
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(data = df %>% filter(outlier),
               colour = "black", shape = 4, size = 5, stroke = 1.2) +
    geom_text(data = df %>% filter(outlier),
              aes(label = nr), vjust = -0.6, size = 4) +
    theme_minimal(base_size = 16) +
    labs(title = paste(y, "vs", x))
}

plot_scatter_outliers(irms_outliers, "d13c_korr", "d15n_korr")
plot_scatter_outliers(irms_outliers, "diameter", "c13_per_dw_ug_per_g")
