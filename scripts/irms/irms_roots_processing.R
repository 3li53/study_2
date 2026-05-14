################ IRMS ROOTS PIPELINE — PROCESSING SCRIPT #######################
load_or_install(c("tidyverse", "tidyr", "dplyr", "writexl", "readr", "dplyr", "stringr", "ggplot2"))

# ---------------------- IMPORT & clean ----------------------------------------
irms_raw <- read.csv("./data_files/IRMS/IRMS_roots_raw.csv", header = TRUE)

cols_to_remove <- c(
  "tray", "position", "X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6",
  "Sample.Number", "Name", "Height..nA.", "N15", "Height..nA..1",
  "X13C", "Weight", "PEAnA", "PEA15N", "PEAnA.1", "PEA13C",
  "gnsnSTD.C.weight", "gnsnSTD.N.weight",
  "X1nA.to.x.mgC.STD", "X1nA.to.x.mgN.STD",
  "δ15Nkorr", "δ13Ckorr" 
)

numeric_cols <- c(
  "nr", "X4.5mg",
  "mg.N.in.sample", "mg.C.in.sample",
  "X.N", "X.C", "C.N",
  "d15Nkorr", "d13Ckorr",
  "Nat.abun..15N.", "atom...15N.", "APE.",
  "N.pr.DW", "X15N.pr.DW", "X15N.pr.N",
  "Nat.abun..13C.", "atom...13C.", "APE..1",
  "C.pr.DW", "X13C.pr.DW", "X13C.pr.C"
)

factor_cols <- c("run", "treatment", "diameter", "layer")

irms <- irms_raw %>%
  filter(!comment %in% c("PEACH", "too small") | is.na(comment)) %>%
  slice(-1) %>%
  mutate(
    across(any_of(factor_cols), as.factor),
    across(any_of(numeric_cols), ~ suppressWarnings(as.numeric(.))),
    comment = as.character(comment)
  ) %>%
  select(-any_of(cols_to_remove)) %>%
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
    atom_pct_15n           = atom...15N.,
    ape_pct_15n            = APE.,
    n_per_dw_mg_per_g      = N.pr.DW,
    n15_per_dw_ug_per_g    = X15N.pr.DW,
    n15_per_n_ug_per_g     = X15N.pr.N,
    nat_abun_13c_atm_pct   = Nat.abun..13C.,
    atom_pct_13c           = atom...13C.,
    ape_pct_13c            = APE..1,
    c_per_dw_mg_per_g      = C.pr.DW,
    c13_per_dw_ug_per_g    = X13C.pr.DW,
    c13_per_c_ug_per_g     = X13C.pr.C
  ) %>%
  mutate(
    veg = str_extract(treatment, "[SB]"),
    cut = str_extract(treatment, "[UC]"),
    wet = str_extract(treatment, "[WD]"),
    beriget = !(run == 1)
  )
glimpse(irms)

#---------------------------- natural abundance outliers -----------------------
natabun <- irms %>% filter(!beriget)

outlier_vars <- outlier_vars <- c(
  "d15n_korr", "d13c_korr", ""
)

natabun_outliers <- iqr_outlier(natabun, group_var = c("layer", "diameter"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

colour_vars <- c(
  "layer"
)

plot_scatter_outliers(natabun_outliers, "d13c_korr", "d15n_korr", "outlier")

# ---------------------- REMOVE OUTLIERS ---------------------------------------
# manually in both dfs based on visual inspection

# remove_outlier_natabun <- c(46)
# natabun <- natabun %>% filter(!nr %in% remove_outlier_natabun)
# irms <- irms %>% filter(!nr %in% remove_outlier_natabun)

#----------------------- calculate baselines -----------------------------------
natabun_means <- calc_isotope_means(natabun, diameter)

# ---------------------- apply corrections on dataframe ------------------------
irms_corrected <- apply_baseline_correction(
  irms_df           = irms,
  natabun_means_df  = natabun_means,
  group_var         = diameter
)

# -------------------- check for outliers full dataset -------------------------
irms_outliers <- iqr_outlier(irms, group_var = c("layer", "diameter"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

colour_vars <- c("diameter")
plot_scatter_outliers(irms_outliers, "d13c_korr", "d15n_korr", "outlier")

# ---------------------- REMOVE OUTLIERS ---------------------------------------
remove_outlier_enriched <- c(239)
irms <- irms %>% filter(!nr %in% remove_outlier_enriched)

# ---------------------- BIOMASS SECTION ---------------------------------------
root_biomass_raw <- read.csv(
  "./data/raw/biomass/pipes_data_roots.csv",
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

h_pos  <- root_biomass_raw[1, ] |> as.character()  # mid / top / bottom
h_frac <- root_biomass_raw[2, ] |> as.character()  # fine / coarse
h_var  <- root_biomass_raw[3, ] |> as.character()  # pipe nr / bag / fresh / dry

fill_right <- function(x) {
  for (i in 2:length(x)) {
    if (is.na(x[i]) || x[i] == "") {
      x[i] <- x[i - 1]
    }
  }
  x
}

h_pos  <- fill_right(h_pos)
h_frac <- fill_right(h_frac)

new_names <- character(length(h_var))

for (i in seq_along(h_var)) {
  if (i <= 4) {
    new_names[i] <- h_var[i]
  } else {
    new_names[i] <- paste(h_pos[i], h_frac[i], h_var[i], sep = "_")
  }
}

root_biomass_raw <- root_biomass_raw[-c(1:3), ]
colnames(root_biomass_raw) <- new_names

root_biomass_raw[ , -c(1:4)] <-
  lapply(root_biomass_raw[ , -c(1:4)], as.numeric)

write_csv(root_biomass_raw, "./data/raw/biomass/roots.csv")

# ------------------------- calculate weights ----------------------------------
root_biomass_raw <- root_biomass_raw |>
  mutate(
    mid_fine    = mid_fine_dry    - mid_fine_bag,
    mid_coarse  = mid_coarse_dry  - mid_coarse_bag,
    top_fine    = top_fine_dry    - top_fine_bag,
    top_coarse  = top_coarse_dry  - top_coarse_bag,
    bottom_fine = bottom_fine_dry - bottom_fine_bag,
    bottom_coarse = bottom_coarse_dry - bottom_coarse_bag
  )
