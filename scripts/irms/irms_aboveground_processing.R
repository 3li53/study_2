###########################################################
# IRMS ROOTS PIPELINE — PROCESSING SCRIPT

library(tidyverse)
library(writexl)

# ---------------------- IMPORT ------------------------------------------------

irms_raw <- read.csv("./data_files/IRMS/IRMS_above_ground_raw.csv", header = TRUE)
glimpse(irms_raw)

# ---------------------- CLEANUP -----------------------------------------------

cols_to_remove <- c(
  "tray", "placement", "X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6",
  "Sample.Number", "Name", "Height..nA.", "N15", "Height..nA..1",
  "X13C", "Weight", "PEAnA", "PEA15N", "PEAnA.1", "PEA13C",
  "gnsnSTD.C.weight", "gnsnSTD.N.weight",
  "X1nA.to.x.mgC.STD", "X1nA.to.x.mgN.STD",
  "δ15Nkorr", "δ13Ckorr" 
)

numeric_cols <- c(
  "nr", "X3.4.mg",
  "mg.N.in.sample", "mg.C.in.sample",
  "X.N", "X.C", "C.N",
  "d15Nkorr", "d13Ckorr",
  "Nat.abun..15N.", "atom...15N.", "APE.",
  "N.pr.DW", "X15N.pr.DW", "X15N.pr.N",
  "Nat.abun..13C.", "atom...13C.", "APE..1",
  "C.pr.DW", "X13C.pr.DW", "X13C.pr.C"
)

factor_cols <- c("run", "treatment", "aboveground")

irms_clean <- irms_raw %>%
  filter(!comment %in% c("PEACH", "too small") | is.na(comment)) %>%
  slice(-1) %>%
  mutate(
    across(any_of(factor_cols), as.factor),
    across(any_of(numeric_cols), ~ suppressWarnings(as.numeric(.))),
    comment = as.character(comment)
  ) %>%
  select(-any_of(cols_to_remove))

glimpse(irms_clean)

# ---------------------- RENAME COLUMNS ----------------------------------------

irms <- irms_clean %>%
  rename(
    mass_3_4_mg            = X3.4.mg,
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
  )
# -----------------------
# -#-#-#-#-#-#-#-#-#-# combine vascular group (new column) -#-#-#-#-#-#-#-#-#-#-

# irms <- irms %>%
#  mutate(
#    aboveground2 = if_else(
#      aboveground %in% c("equisetum", "graminoid", "salix_leaf", "salix_stem"),
#      "vascular",
#      as.character(aboveground)
#    ),
#    aboveground2 = factor(aboveground2)
#  )
# ???????????????????????????????????????????

# ---------------------- TREATMENT PARSING -------------------------------------

irms <- irms %>%
  mutate(
    veg = str_extract(treatment, "[SB]"),
    cut = str_extract(treatment, "[UC]"),
    wet = str_extract(treatment, "[WD]"),
    beriget = !(run == 1 | cut == "C")
  )
glimpse(irms)

# ---------------------- OUTLIER FUNCTION --------------------------------------

iqr_outlier <- function(df, group_var, numeric_vars) {
  df %>%
    group_by(across(all_of(group_var))) %>%
    mutate(across(
      all_of(numeric_vars),
      ~ (. < quantile(., 0.25, na.rm = TRUE) - 1.5 * IQR(., na.rm = TRUE)) |
        (. > quantile(., 0.75, na.rm = TRUE) + 1.5 * IQR(., na.rm = TRUE)),
      .names = "outlier_{.col}"
    )) %>%
    ungroup()
}

# ---------------------- NATURAL ABUNDANCE outliers-----------------------------

natabun <- irms %>% filter(!beriget)

outlier_vars <- c(
  "d15n_korr", "d13c_korr",
  "n15_per_dw_ug_per_g", "c13_per_dw_ug_per_g",
  "n15_per_n_ug_per_g", "c13_per_c_ug_per_g"
)

natabun_outliers <- iqr_outlier(natabun, group_var = c("aboveground"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

# ---------------------- REMOVE natabun OUTLIERS -------------------------------
# manually in both dfs based on visual inspection

remove_outlier_natabun <- c(176, 158, 159, 163, 157, 42,
                    210, 1, 2, 217, 228)

natabun <- natabun %>% filter(!nr %in% remove_outlier_natabun)
irms <- irms %>% filter(!nr %in% remove_outlier_natabun)

# ---------------------- BASELINE ----------------------------------------------

natabun_means <- natabun %>%
  group_by(aboveground) %>%
  summarise(
    avg_d15 = mean(d15n_korr, na.rm = TRUE),
    se_d15  = sd(d15n_korr, na.rm = TRUE) / sqrt(n()),
    ymin    = (mean(d15n_korr, na.rm = TRUE)) - (sd(d15n_korr, na.rm = TRUE) / sqrt(n())),
    ymax    = (mean(d15n_korr, na.rm = TRUE)) + (sd(d15n_korr, na.rm = TRUE) / sqrt(n())),
    avg_d13 = mean(d13c_korr, na.rm = TRUE),
    se_d13  = sd(d13c_korr, na.rm = TRUE) / sqrt(n()),
    xmin    = (mean(d13c_korr, na.rm = TRUE)) - (sd(d13c_korr, na.rm = TRUE) / sqrt(n())),
    xmax    = (mean(d13c_korr, na.rm = TRUE)) + (sd(d13c_korr, na.rm = TRUE) / sqrt(n())),
    .groups = "drop"
  )

R15 <- 0.003676
R13 <- 0.011237

# define function to calculate atom pct
atom_pct <- function(delta, Rstd) {
  Rsample <- Rstd * (1 + (delta / 1000))
  100 * (Rsample / (1 + Rsample))
}

# ---------------------- APPLY BASELINE CORRECTION -----------------------------

irms <- irms %>%
  left_join(natabun_means, by = c("aboveground")) %>%
  mutate(
    nat_abun_15n_atm_pct = atom_pct(avg_d15,   R15),
    nat_abun_13c_atm_pct = atom_pct(avg_d13,   R13),
            atom_pct_15n = atom_pct(d15n_korr, R15),
            atom_pct_13c = atom_pct(d13c_korr, R13),
             ape_pct_15n = atom_pct_15n - nat_abun_15n_atm_pct,
             ape_pct_13c = atom_pct_13c - nat_abun_13c_atm_pct,
     n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000,
     n15_per_n_ug_per_g  = mg_n              * (ape_pct_15n / 100) * 1000,
     c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,
     c13_per_c_ug_per_g  = mg_c              * (ape_pct_13c / 100) * 1000,
  ) %>% 
  select(-c(avg_d15, ymin, ymax, avg_d13, se_d13, xmin, xmax))

# ---------------------- OUTLIERS ON CORRECTED DATA ----------------------------

irms_outliers <- iqr_outlier(irms, group_var = c("aboveground"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

# ---------------------- REMOVE ADDITIONAL OUTLIERS ----------------------------

remove_outlier_enriched <- c(164, 165, 182, 113, 160, 161, 162, 107, 156)

irms <- irms %>% filter(!nr %in% remove_outlier_enriched)

# ---------------------- BIOMASS SECTION ---------------------------------------
#this section includes the biomass data into further calculations

# Attach pipe.nr into IRMS
irms <- irms %>%
  left_join(
    aboveground_biomass_raw %>%
      transmute(
        run,
        veg = as.character(veg_type),
        cut = as.character(cut_treat),
        wet = as.character(wet_treat),
        pipe.nr
      ),
    by = c("run", "veg", "cut", "wet")
  )

aboveground <- left_join(
  irms,
  aboveground_biomass,
  by = "pipe.nr"
) %>%
  mutate(
    run = run.x,
    veg = veg.x,
    cut = cut.x,
    wet = wet.x,
    beriget = !(run == "1" | cut == "C")
  ) %>%
  select(-matches("\\.x$|\\.y$")) %>%
  mutate(
    treatment = paste0(veg, cut, wet)
  )
# add a vascular column for added weight of equisetum, graminoids, salix

aboveground <- aboveground %>%
  mutate(
    vascular_weight = rowSums(
      across(c(graminoid_weight,
               equ_weight,
               stem_weight,
               leaves_weight)),
      na.rm = TRUE
    )
  )

# ---------------------- pools and nitrogen in total biomass -------------------

aboveground <- aboveground %>%
  mutate(
    component_weight = case_when(
      aboveground == "bryophyte"   ~ bryophyte_weight,
      aboveground == "graminoid"   ~ graminoid_weight,
      aboveground == "equisetum"   ~ equ_weight,
      aboveground == "lichen"      ~ lichen_weight,
      aboveground == "salix_stem"  ~ stem_weight,
      aboveground == "salix_leaf"  ~ leaves_weight,
      #-#-#-#-#-#-#-#-#-#???????????
    ###  aboveground2 == "vascular"    ~ vascular_weight,
      #-#-#-#-#-#-#-#-#-#??????????
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    pulje_15N_ug = n15_per_dw_ug_per_g * component_weight,
    pulje_13C_ug = c13_per_dw_ug_per_g * component_weight,
    pulje_15N_mg = pulje_15N_ug / 1000,
    pulje_13C_mg = pulje_13C_ug / 1000,
    pulje_N_mg   = (pct_n / 100) * component_weight
  )

# ---------------------- BIOMASS OUTLIERS --------------------------------------

biomass_outliers <- iqr_outlier(
  aboveground, group_var = c("aboveground"), numeric_vars = "component_weight"
) %>% mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

remove_outlier_biomass <- c(97, 24, 105, 98, 65, 8, 177, 175)
aboveground <- aboveground %>% filter(!nr %in% remove_outlier_biomass)

# ----------------------- pools by area ----------------------------------------

aboveground <- aboveground %>% 
  left_join(
    read.csv("./data_files/pipes_data/aboveground/pipes_data_aboveground_field.csv") %>%
      filter(run != "trial2") %>%
      rename(pipe_diameter_cm = Ø) %>%
      select(nr, pipe_diameter_cm),
    by = c("pipe.nr" = "nr")
  ) %>% 
  mutate(
    pipe_area_m2 = pi * (pipe_diameter_cm / 2)^2 / 10000
  )

aboveground <- aboveground %>%
  mutate(
    pulje_15N_mg_m2 = pulje_15N_mg / pipe_area_m2,
    pulje_13C_mg_m2 = pulje_13C_mg / pipe_area_m2,
    pulje_N_mg_m2   = pulje_N_mg   / pipe_area_m2
  )

#------------------------- % RECOVERY ------------------------------------------

added_15N_mg_m2 <- 88   # <-- put your TRUE added 15N here
added_13C_mg_m2 <- 152.5   # <-- put your TRUE added 13C here

aboveground <- aboveground %>%
  mutate(
    recovery_15N_pct = (pulje_15N_mg_m2 / added_15N_mg_m2) * 100,
    recovery_13C_pct = (pulje_13C_mg_m2 / added_13C_mg_m2) * 100
  )

