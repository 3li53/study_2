irms_roots_import <- read.csv("./data_files/IRMS/IRMS_roots_raw.csv", 
                           header = TRUE)

library(dplyr)

summary(irms_roots_import)
head(irms_roots_import)
glimpse(irms_roots_import)

# ---------- cleanup ----------

irms_roots_raw <- irms_roots_import %>%
  filter(!comment %in% c("PEACH", "too small") | is.na(comment)) %>%
  slice(-1) %>%
  mutate(
    # factors
    across(c(run, treatment, layer, diameter), as.factor),
    # numeric
    across(c(nr, X4.5mg, 
             mg.N.in.sample, mg.C.in.sample,
             X.N, X.C, C.N, d15Nkorr, d13Ckorr, 
             Nat.abun..15N., δ15Nkorr, atom...15N., APE., N.pr.DW, X15N.pr.DW, X15N.pr.N, 
             Nat.abun..13C., δ13Ckorr, atom...13C., APE..1, C.pr.DW, X13C.pr.DW, X13C.pr.C),
           ~ as.numeric(.)
    ),
    # character
    across(c(comment), as.character)
  ) %>% 
  select(
    -tray, -position, -X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6,
    -Sample.Number, -Name, -Height..nA., -N15, -Height..nA..1, -X13C, -Weight,
    -PEAnA, -PEA15N, - PEAnA.1, -PEA13C, -gnsnSTD.C.weight, -gnsnSTD.N.weight, -X1nA.to.x.mgC.STD, -X1nA.to.x.mgN.STD
  )

glimpse(irms_roots_raw)

# rename column names and add treatment columns and label status
irms_roots_names <- irms_roots_raw %>%
  rename(
    nr                  = nr,
    run                 = run,
  # tinytag             = tt,
    treatment           = treatment,
    diameter            = diameter,
    layer               = layer,
    
    mg_n                = mg.N.in.sample,
    mg_c                = mg.C.in.sample,
    
    mass_4_5_mg         = `X4.5mg`, 
    
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
    # treatment specific columns
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
  mutate(
    # label status
    beriget = !(run == 1)
  )

# look at the data

library(ggplot2)

ggplot(irms_roots_names, aes(x = layer, y = d15n_korr, colour = beriget)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(title = "Distribution of δ15N", y = "δ15N", x = "soil layer")
ggplot(irms_roots_names, aes(x = layer, y = d13c_korr, colour = beriget)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(title = "Distribution of δ13C", y = "δ13C", x = "soil layer")

ggplot(irms_roots_names, aes(x = diameter, y = d15n_korr, colour = beriget)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(title = "Distribution of δ15N", y = "δ15N", x = "soil layer")
ggplot(irms_roots_names, aes(x = diameter, y = d13c_korr, colour = beriget)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(title = "Distribution of δ13C", y = "δ13C", x = "soil layer")

# ---------- find nat_abun outliers before correcting variables ----------

natabun <- irms_roots_names %>%
  filter(!beriget)

ggplot(natabun, aes(x = layer, y = d15n_korr)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) 
ggplot(natabun, aes(x = layer, y = d13c_korr)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) 

ggplot(irms_roots_names, aes(x = diameter, y = d15n_korr)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) 
ggplot(irms_roots_names, aes(x = diameter, y = d13c_korr)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) 

outliers_iqr <- natabun %>%
  group_by(diameter) %>%
  mutate(
    d15n_low  = quantile(d15n_korr, 0.25) - 1.5 * IQR(d15n_korr),
    d15n_high = quantile(d15n_korr, 0.75) + 1.5 * IQR(d15n_korr),
    d13c_low  = quantile(d13c_korr, 0.25) - 1.5 * IQR(d13c_korr),
    d13c_high = quantile(d13c_korr, 0.75) + 1.5 * IQR(d13c_korr),
    outlier_d15n = d15n_korr < d15n_low | d15n_korr > d15n_high,
    outlier_d13c = d13c_korr < d13c_low | d13c_korr > d13c_high
  )

outliers_iqr %>% 
  filter(outlier_d15n | outlier_d13c)

# Add a single outlier column "outlier"
outliers_iqr <- outliers_iqr %>%
  mutate(outlier = outlier_d15n | outlier_d13c)

# ---------- plot the natural abundance ----------

ggplot(outliers_iqr,
       aes(x = d13c_korr, y = d15n_korr, colour = diameter)) +
  # Non-outliers: above ground specific coloured circles
  geom_point(data = subset(outliers_iqr, !outlier),
             size = 3, alpha = 0.9) +
  # Outliers: root specific coloured x
  geom_point(data = subset(outliers_iqr, outlier),
             aes(colour = diameter),
             shape = 4, size = 3.5, stroke = 1.2) +
  # Label only outliers (toggle to label all)
  geom_text(
    data = subset(outliers_iqr, outlier),
    aes(label = nr),
    nudge_y = 0.3, size = 5, show.legend = FALSE, colour = "black"
  ) +
  labs(
    title = "Scatterplot of unlabelled outliers δ13C vs δ15N                      X = Tukey's fences outliers",
    x = "δ13C",
    y = "δ15N",
    colour = "diameter"
  ) +
  theme_minimal()+
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

ggplot(outliers_iqr,
       aes(x = diameter, y = d13c_korr, colour = beriget)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(outliers_iqr, outlier_d13c),
    aes(colour = beriget),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(outliers_iqr, outlier_d13c),
    aes(label = nr),
    nudge_y = .4, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
ggplot(outliers_iqr,
       aes(x = diameter, y = d15n_korr, colour = beriget)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(outliers_iqr, outlier_d15n),
    aes(colour = beriget),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(outliers_iqr, outlier_d15n),
    aes(label = nr),
    nudge_y = .3, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
################################################################################
# ---------- correct root-specific natural abundances ----------

# create nat_abundance averages from run 1 
avg_table_natabun <- natabun %>%
  group_by(diameter, layer) %>%
  summarise(
    avg_d15 = mean(d15n_korr, na.rm = TRUE),
    n_d15   = sum(!is.na(d15n_korr)),
    avg_d13 = mean(d13c_korr, na.rm = TRUE),
    n_d13   = sum(!is.na(d13c_korr)),
    n_both  = sum(!is.na(d15n_korr) & !is.na(d13c_korr)),
    n_rows  = n(),
    .groups = "drop"
  )

# ---------- correct δ → atom% ----------
# Correct isotope ratios:
R15 <- 0.003676       # 15N/14N ratio of air N2
R13 <- 0.01118        # 13C/12C ratio for VPDB

# helper function for atom%
atom_pct <- function(delta, Rstd) {
  Rsample <- Rstd * (1 + delta/1000)
  atom_frac <- Rsample / (1 + Rsample)
  atom_frac * 100
}

# apply correct baseline atom% values
irms_roots_natabund_corrected <- irms_roots_names %>%
  left_join(avg_table_natabun, by = c("diameter", "layer")) %>%
  mutate(
    nat_abun_15n_atm_pct = atom_pct(avg_d15, R15),
    nat_abun_13c_atm_pct = atom_pct(avg_d13, R13)
  )

# ---------- correct atom % for samples ----------

irms_roots_atompct_corrected <- irms_roots_natabund_corrected %>%
  mutate(
    atom_pct_15n = atom_pct(d15n_korr, R15),
    atom_pct_13c = atom_pct(d13c_korr, R13)
  )

# ---------- correct APE ----------

irms_roots_ape_corrected <- irms_roots_atompct_corrected %>%
  mutate(
    ape_pct_15n = atom_pct_15n - nat_abun_15n_atm_pct,
    ape_pct_13c = atom_pct_13c - nat_abun_13c_atm_pct
  )

# ---------- correct full dataframe: convert APE → µg/g ----------

irms_roots_corrected <- irms_roots_ape_corrected %>%
  mutate(
    n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000,
    n15_per_n_ug_per_g  = mg_n * (ape_pct_15n / 100) * 1000,
    
    c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,
    c13_per_c_ug_per_g  = mg_c * (ape_pct_13c / 100) * 1000
  )
##########################################################################################



library(writexl)
#write_xlsx(irms_roots_corrected, "./data_files/IRMS/IRMS_roots_corrected.xlsx")

# ---------- chek for outliers ----------

irms_outliers_iqr <- irms_roots_corrected %>%
  group_by(diameter) %>%
  mutate(
    d15n_low  = quantile(d15n_korr, 0.25) - 1.5 * IQR(d15n_korr),
    d15n_high = quantile(d15n_korr, 0.75) + 1.5 * IQR(d15n_korr),
    
    d13c_low  = quantile(d13c_korr, 0.25) - 1.5 * IQR(d13c_korr),
    d13c_high = quantile(d13c_korr, 0.75) + 1.5 * IQR(d13c_korr),
    
    n15_dw_low  = quantile(n15_per_dw_ug_per_g, 0.25) - 1.5 * IQR(n15_per_dw_ug_per_g),
    n15_dw_high = quantile(n15_per_dw_ug_per_g, 0.75) + 1.5 * IQR(n15_per_dw_ug_per_g),
    
    c13_dw_low  = quantile(c13_per_dw_ug_per_g, 0.25) - 1.5 * IQR(c13_per_dw_ug_per_g),
    c13_dw_high = quantile(c13_per_dw_ug_per_g, 0.75) + 1.5 * IQR(c13_per_dw_ug_per_g),
    
    n15_n_low   = quantile(n15_per_n_ug_per_g, 0.25) - 1.5 * IQR(n15_per_n_ug_per_g),
    n15_n_high  = quantile(n15_per_n_ug_per_g, 0.75) + 1.5 * IQR(n15_per_n_ug_per_g),
    
    c13_c_low   = quantile(c13_per_c_ug_per_g, 0.25) - 1.5 * IQR(c13_per_c_ug_per_g),
    c13_c_high  = quantile(c13_per_c_ug_per_g, 0.75) + 1.5 * IQR(c13_per_c_ug_per_g),
    
    outlier_d15n = d15n_korr < d15n_low | d15n_korr > d15n_high,
    outlier_d13c = d13c_korr < d13c_low | d13c_korr > d13c_high, 
    
    outlier_n15_dw = n15_per_dw_ug_per_g < n15_dw_low | n15_per_dw_ug_per_g > n15_dw_high,
    outlier_c13_dw = c13_per_dw_ug_per_g < c13_dw_low | c13_per_dw_ug_per_g > c13_dw_high,  
    
    outlier_n15_n = n15_per_n_ug_per_g < n15_n_low | n15_per_n_ug_per_g > n15_n_high,
    outlier_c13_c = c13_per_c_ug_per_g < c13_c_low | c13_per_c_ug_per_g > c13_c_high
  )

irms_outliers_iqr %>% 
  filter(outlier_d15n | outlier_d13c | outlier_n15_dw | outlier_c13_dw | outlier_n15_n |outlier_c13_c)

# Add a single outlier column "outlier"
irms_outliers_iqr <- irms_outliers_iqr %>%
  mutate(outlier = outlier_d15n | outlier_d13c | outlier_n15_dw | outlier_c13_dw | outlier_n15_n | outlier_c13_c)
irms_outliers_iqr %>% 
  filter(outlier_d15n)

# ---------- plot outliers----------

ggplot(irms_roots_corrected,
       aes(x = diameter, y = c13_per_dw_ug_per_g, colour = beriget)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(irms_outliers_iqr, outlier_c13_dw),
    aes(colour = beriget),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_c13_dw),
    aes(label = nr),
    nudge_y = .7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
