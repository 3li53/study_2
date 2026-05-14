
# IRMS ROOTS CLEAN PIPELINE

library(tidyverse)
library(writexl)

# Import

irms_raw <- read.csv("./data_files/IRMS/IRMS_above_ground_raw.csv", header = TRUE)
glimpse(irms_raw)

# CLEANUP

cols_to_remove <- c(
  "tray", "placement", "X", "X.1", "X.2", "X.3", "X.4", "X.5", "X.6",
  "Sample.Number", "Name", "Height..nA.", "N15", "Height..nA..1",
  "X13C", "Weight", "PEAnA", "PEA15N", "PEAnA.1", "PEA13C",
  "gnsnSTD.C.weight", "gnsnSTD.N.weight",
  "X1nA.to.x.mgC.STD", "X1nA.to.x.mgN.STD"
)

numeric_cols <- c(
  "nr", "X3.4.mg",
  "mg.N.in.sample",   "mg.C.in.sample",
  "X.N", "X.C", "C.N",
  "d15Nkorr", "d13Ckorr",
  "Nat.abun..15N.", "δ15Nkorr", "atom...15N.", "APE.",
  "N.pr.DW", "X15N.pr.DW", "X15N.pr.N",
  "Nat.abun..13C.", "δ13Ckorr", "atom...13C.", "APE..1",
  "C.pr.DW", "X13C.pr.DW", "X13C.pr.C"
)

factor_cols <- c("run", "treatment", "aboveground")

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
    mass_4_5_mg            = X3.4.mg,
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
    beriget = !(run == 1 | cut == "C")
  )
glimpse(irms)


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

# NATURAL ABUNDANCE

natabun <- irms %>% filter(!beriget)
glimpse(natabun)

# OUTLIERS ON natabun DATA

outlier_vars <- c(
  "d15n_korr", "d13c_korr",
  "n15_per_dw_ug_per_g", "c13_per_dw_ug_per_g",
  "n15_per_n_ug_per_g", "c13_per_c_ug_per_g"
)

natabun_outliers <- iqr_outlier(natabun, group_var = c("aboveground"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))
glimpse(natabun_outliers)

# plot natural abundance
# with outliers 

plot_scatter_outliers <- function(df, x, y, outlier_col, colour = "aboveground") {
  
  # check that the outlier column exists
  if (!outlier_col %in% names(df)) {
    stop(paste("Outlier column", outlier_col, "not found in dataframe"))
  }
  
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_point(size = 3, alpha = 0.7) +
    # highlight outliers
    geom_point(
      data = df[df[[outlier_col]] == TRUE, ],
      colour = "black", shape = 4, size = 5, stroke = 1.2
    ) +
        geom_text(
      data = df[df[[outlier_col]] == TRUE, ],
      aes(label = nr),
      vjust = -0.6, size = 4, colour = "black"
    ) +
    theme_minimal(base_size = 16)
}

plot_scatter_outliers(natabun_outliers, x = "d13c_korr", y = "d15n_korr", outlier_col = "outlier")
plot_scatter_outliers(natabun_outliers, x = "aboveground", y = "d15n_korr", outlier_col = "outlier_d15n_korr")
plot_scatter_outliers(natabun_outliers, x = "aboveground", y = "d13c_korr", outlier_col = "outlier_d13c_korr")

# remove specific outliers

remove_outlier <- c(176, 158, 159, 163, 157, 42,
                    210, 1, 2, 217, 228)
natabun <- natabun_outliers %>% filter(!nr %in% remove_outlier)
irms <- irms %>% filter(!nr %in% remove_outlier)

# NATURAL ABUNDANCE BASELINE

natabun_means <- natabun %>%
  group_by(aboveground) %>%
  summarise(
    avg_d15 = mean(d15n_korr, na.rm = TRUE),
    avg_d13 = mean(d13c_korr, na.rm = TRUE),
    .groups = "drop"
  )
glimpse(natabun_means)

# ISOTOPE FORMULA HELPERS

R15 <- 0.003676   # air N2
R13 <- 0.011237    # VPDB

atom_pct <- function(delta, Rstd) {
  Rsample <- Rstd * (1 + (delta / 1000))
  100 * (Rsample / (1 + Rsample))
}

# APPLY BASELINE CORRECTIONS

irms <- irms %>%
  left_join(natabun_means, by = c("aboveground")) %>%
  mutate(
    nat_abun_15n_atm_pct = atom_pct(avg_d15, R15), # Replace nat abundance atm pct with above ground–specific values
    nat_abun_13c_atm_pct = atom_pct(avg_d13, R13), # Replace nat abundance atm pct with above ground–specific values
    
    atom_pct_15n = atom_pct(d15n_korr, R15), # corrected atom %
    atom_pct_13c = atom_pct(d13c_korr, R13), # corrected atom %
    
    ape_pct_15n = atom_pct_15n - nat_abun_15n_atm_pct, # correct APE
    ape_pct_13c = atom_pct_13c - nat_abun_13c_atm_pct, # correct APE
    
    n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000, #15N pr DW
    n15_per_n_ug_per_g  = mg_n              * (ape_pct_15n / 100) * 1000, #15N pr N
    
    c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000, #13C pr DW
    c13_per_c_ug_per_g  = mg_c              * (ape_pct_13c / 100) * 1000  #13C pr C
  )

# OUTLIERS ON CORRECTED DATA

irms_outliers <- iqr_outlier(irms, group_var = c("aboveground"), numeric_vars = outlier_vars) %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

# PLOT THE DATA

plot_scatter_outliers(irms_outliers, x = "c13_per_dw_ug_per_g", y = "n15_per_dw_ug_per_g", outlier_col = "outlier")
plot_scatter_outliers(irms_outliers, x = "aboveground", y = "c13_per_dw_ug_per_g", outlier_col = "outlier_c13_per_dw_ug_per_g")
plot_scatter_outliers(irms_outliers, x = "aboveground", y = "n15_per_dw_ug_per_g", outlier_col = "outlier_n15_per_dw_ug_per_g")

# boxplot

plot_box <- function(df, x, y, colour = "beriget", title = NULL) {
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    theme_minimal(base_size = 16) +
    labs(title = title, x = x, y = y)
}

plot_box(irms, "aboveground", "n15_per_dw_ug_per_g")
plot_box(irms, "aboveground", "c13_per_dw_ug_per_g")

# remove specific outliers

remove_outlier <- c(164, 165, 182, 113, 160, 161, 162, 107, 156)
irms <- irms_outliers %>% filter(!nr %in% remove_outlier)

plot_scatter_outliers(irms, x = "aboveground", y = "c13_per_dw_ug_per_g", outlier_col = "outlier_c13_per_dw_ug_per_g")
plot_scatter_outliers(irms, x = "aboveground", y = "n15_per_dw_ug_per_g", outlier_col = "outlier_n15_per_dw_ug_per_g")

# summary df

irms_summary <- irms %>%
  group_by(aboveground, beriget) %>%
  summarise(
    mean_n15_per_dw_ug_per_g = mean(n15_per_dw_ug_per_g, na.rm = TRUE),
    se_n15_per_dw_ug_per_g   = sd(n15_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    mean_c13_per_dw_ug_per_g = mean(c13_per_dw_ug_per_g, na.rm = TRUE),
    se_c13_per_dw_ug_per_g   = sd(c13_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# truncate negative values only for plotting

irms_plot <- irms %>%
  mutate(n15_per_dw_ug_per_g = pmax(n15_per_dw_ug_per_g, 0))

library(ggplot2)
library(patchwork)

# plot 1: natural abundance
p1 <- ggplot(subset(irms_summary, beriget == FALSE),
             aes(x = aboveground,
                 y = mean_n15_per_dw_ug_per_g,
                 fill = beriget)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(
      ymin = mean_n15_per_dw_ug_per_g - se_n15_per_dw_ug_per_g,
      ymax = mean_n15_per_dw_ug_per_g + se_n15_per_dw_ug_per_g
    ),
    width = 0.2
  ) +
  scale_fill_manual(values = "grey70", guide = "none") +
  labs(
    x = "Vegetation",
    y = "Natural abundance (15N per DW, µg/g)"
  ) +
  theme_classic(base_size = 20)

# plot 2: enriched
p2 <- ggplot(subset(irms_summary, beriget == TRUE),
             aes(x = aboveground,
                 y = mean_n15_per_dw_ug_per_g,
                 fill = beriget)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(
      ymin = mean_n15_per_dw_ug_per_g - se_n15_per_dw_ug_per_g,
      ymax = mean_n15_per_dw_ug_per_g + se_n15_per_dw_ug_per_g
    ),
    width = 0.2
  ) +
  scale_fill_manual(values = "grey30", guide = "none") +
  labs(
    x = "Vegetation",
    y = "Enriched (15N per DW, µg/g)"
  ) +
  theme_classic(base_size = 20)

# layout: stacked vertically with separate y axes
p1 | p2



ggplot(irms_summary, aes(x = aboveground, y = mean_n15_per_dw_ug_per_g, fill = beriget)) +
  geom_col(colour = "black", position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = mean_n15_per_dw_ug_per_g - se_n15_per_dw_ug_per_g, ymax = mean_n15_per_dw_ug_per_g + se_n15_per_dw_ug_per_g),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  theme_classic(base_size = 20) +
  scale_fill_manual(
    name = "Stable isotopes",
    values = c(
      "FALSE" = "grey100",
      "TRUE"  = "grey60"
    ),
    labels = c("Natural abundance", "Enriched"),
  ) +
  labs(x = "Vegetation", y = "15N pr DW (ug per g)") +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0.5, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = c("bryophyte" = "Bryophytes",
                              "equisetum" = "Equisetum",
                              "graminoid" = "Graminoids",
                              "lichen" = "Lichen",
                              "salix_leaf" = "Salix leaves",
                              "salix_stem" = "Salix stems"))



                        # POOLS STOCKS PULJER
# The total amount of tracer isotope (¹⁵N or ¹³C) contained in a biological compartment, usually expressed per unit dry mass.

# Stock15N = ( μg 15N / g DW ) × (g DW of the compartment)
# Stock13C = ( μg 13C / g DW ) × (g DW of the compartment)


library(tidyverse)

aboveground_biomass_raw <- read.csv(
  "./data_files/pipes_data/aboveground/pipes_data_aboveground.csv",
  header = TRUE,
  skip = 1
)

aboveground_biomass_raw <- aboveground_biomass_raw %>%
  filter(!is.na(run)) %>% # Remove garbage rows 
  select(-aboveground_weight_tot_g) %>%  # Remove the WRONG duplicate total-weight column (character)
  mutate( # Ensure joining keys match IRMS types
    run = as.character(run),
    veg_type = as.character(veg_type),
    cut_treat = as.character(cut_treat),
    wet_treat = as.character(wet_treat)
  )

# Create clean biomass table with ONLY needed columns
aboveground_biomass <- aboveground_biomass_raw %>%
  mutate(
    veg      = veg_type,
    cut      = cut_treat,
    wet      = wet_treat,
    beriget  = !(run == "1" | cut == "C"),
    treatment = paste0(veg, cut, wet)
  ) %>%
  rename(
    aboveground_weight_tot_g = aboveground_weight_tot_g.1,
    weight_sorted_excl_bag_g = weight_sorted_excl.bag_g
  ) %>%
  select(
    pipe.nr, run, veg, cut, wet, treatment, beriget,
    aboveground_weight_tot_g,
    equ_weight, graminoid_weight, bryophyte_weight,
    lichen_weight, stem_weight, leaves_weight
  )

# Attach pipe.nr using unique key: run + veg + cut + wet
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
  # overwrite with the IRMS treatment variables
  mutate(
    run = run.x,
    veg = veg.x,
    cut = cut.x,
    wet = wet.x,
    beriget = !(run == "1" | cut == "C")
  ) %>%
  select(-matches("\\.x$|\\.y$")) %>%   # remove suffix columns
  mutate(
    treatment = paste0(veg, cut, wet)
  )

aboveground <- aboveground %>%
  mutate(
    component_weight = case_when(
      aboveground == "bryophyte"   ~ bryophyte_weight,
      aboveground == "graminoid"   ~ graminoid_weight,
      aboveground == "equisetum"   ~ equ_weight,
      aboveground == "lichen"      ~ lichen_weight,
      aboveground == "salix_stem"  ~ stem_weight,
      aboveground == "salix_leaf"  ~ leaves_weight,
      TRUE ~ NA_real_
    )
  ) %>% 
  mutate(
    pulje_15N_ug = n15_per_dw_ug_per_g * component_weight,
    pulje_13C_ug = c13_per_dw_ug_per_g * component_weight,
    pulje_15N_mg = pulje_15N_ug / 1000,
    pulje_13C_mg = pulje_13C_ug / 1000
  )

# FIND THE BIOMASS OUTLIERS

biomass_outliers <- iqr_outlier(aboveground, group_var = c("aboveground"), numeric_vars = "component_weight") %>%
  mutate(outlier = if_any(starts_with("outlier_"), ~ .x))

plot_scatter_outliers(biomass_outliers, x = "aboveground", y = "component_weight", outlier_col = "outlier_component_weight")

remove_outlier <- c(97, 24, 105, 98, 65, 8, 177, 175)
aboveground <- biomass_outliers %>% filter(!nr %in% remove_outlier)

plot_scatter_outliers(aboveground, x = "aboveground", y = "component_weight", outlier_col = "outlier_component_weight")

# summary df

pool_summary <- aboveground %>%
  group_by(aboveground, beriget) %>%
  summarise(
    mean_15N = mean(pulje_15N_ug, na.rm = TRUE),
    se_15N   = sd(pulje_15N_ug, na.rm = TRUE) / sqrt(n()),
    mean_13C = mean(pulje_13C_ug, na.rm = TRUE),
    se_13C   = sd(pulje_13C_ug, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 15N Pools
ggplot(pool_summary, aes(x = aboveground, y = mean_15N, fill = beriget)) +
  geom_col(colour = "black", position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = mean_15N - se_15N, ymax = mean_15N + se_15N),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(y = "15N pool (µg)") +
  theme_minimal(base_size = 16)

# 13C Pools
ggplot(pool_summary, aes(x = aboveground, y = mean_13C, fill = beriget)) +
  geom_col(colour = "black", position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = mean_13C - se_13C, ymax = mean_13C + se_13C),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(y = "13C pool (µg)") +
  theme_minimal(base_size = 16)


# RAW POINTS OVERLAY (optional)

ggplot(aboveground, aes(x = aboveground, y = pulje_15N_ug, colour = beriget)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  stat_summary(fun = mean, geom = "col", fill = "grey80", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(y = "15N pool (µg)") +
  theme_minimal(base_size = 16)

ggplot(aboveground, aes(x = aboveground, y = pulje_13C_ug, colour = beriget)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  stat_summary(fun = mean, geom = "col", fill = "grey80", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(y = "13C pool (µg)") +
  theme_minimal(base_size = 16)

# Total pools per pipe

total_aboveground_pools <- aboveground %>%
  group_by(pipe.nr, run, veg, cut, wet, treatment, beriget) %>%
  summarise(
    total_15N_ug = sum(pulje_15N_ug, na.rm = TRUE),
    total_13C_ug = sum(pulje_13C_ug, na.rm = TRUE),
    total_15N_mg = sum(pulje_15N_mg, na.rm = TRUE),
    total_13C_mg = sum(pulje_13C_mg, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(total_aboveground_pools,
       aes(x = treatment, y = total_15N_ug, fill = beriget)) +
  stat_summary(fun = mean, geom = "col", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 16) +
  labs(y = "Total aboveground 15N pool (µg)")
ggplot(total_aboveground_pools,
       aes(x = treatment, y = total_13C_ug, fill = beriget)) +
  stat_summary(fun = mean, geom = "col", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 16) +
  labs(y = "Total aboveground 13C pool (µg)")
ggplot(total_aboveground_pools,
       aes(x = treatment, y = total_15N_ug, colour = beriget)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "col", fill = "grey85", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 16) +
  labs(y = "Total aboveground 15N (µg)")
ggplot(total_aboveground_pools,
       aes(x = treatment, y = total_13C_ug, colour = beriget)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "col", fill = "grey85", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 16) +
  labs(y = "Total aboveground 13C (µg)")

