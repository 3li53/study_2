irms_above_ground_import <- read.csv("./data_files/IRMS/IRMS_above_ground_raw.csv",
                      sep = ",", header = TRUE)

library(dplyr)

summary(irms_above_ground_import)
head(irms_above_ground_import)
glimpse(irms_above_ground_import)

# ---------- cleanup ----------

irms_above_ground_raw <- irms_above_ground_import %>%
  filter(!comment %in% c("PEACH", "too small") | is.na(comment)) %>%
  slice(-1) %>%
  mutate(
    # factors
    across(c(run, tt, treatment, aboveground), as.factor),
    # numeric
    across(c(nr, X3.4.mg, 
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
    -tray, -placement, -X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6,
    -Sample.Number, -Name, -Height..nA., -N15, -Height..nA..1, -X13C, -Weight,
    -PEAnA, -PEA15N, - PEAnA.1, -PEA13C, -gnsnSTD.C.weight, -gnsnSTD.N.weight, -X1nA.to.x.mgC.STD, -X1nA.to.x.mgN.STD
    )

glimpse(irms_above_ground_raw)

# rename column names and add treatment columns and label status
irms_above_ground_names <- irms_above_ground_raw %>%
  rename(
    nr                  = nr,
    run                 = run,
    tinytag             = tt,
    treatment           = treatment,
    above_ground        = aboveground,
    
    mg_n                = mg.N.in.sample,
    mg_c                = mg.C.in.sample,

    mass_3_4_mg         = `X3.4.mg`, 
    
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
    beriget = !(run == 1 | cut == "C")
  )

# look at the data

library(ggplot2)

ggplot(irms_above_ground_names, aes(x = above_ground, y = d15n_korr, colour = beriget)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(title = "Distribution of δ15N", y = "δ15N", x = "Above ground")
ggplot(irms_above_ground_names, aes(x = above_ground, y = d13c_korr, colour = beriget)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(title = "Distribution of δ13C", y = "δ13C", x = "Above ground")

# ---------- find nat_abun outliers before correcting variables ----------

natabun <- irms_above_ground_names %>%
  filter(!beriget)

outliers_iqr <- natabun %>%
  group_by(above_ground) %>%
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
       aes(x = d13c_korr, y = d15n_korr, colour = above_ground)) +
  # Non-outliers: above ground specific coloured circles
  geom_point(data = subset(outliers_iqr, !outlier),
             size = 3, alpha = 0.9) +
  # Outliers: above ground specific coloured x
  geom_point(data = subset(outliers_iqr, outlier),
             aes(colour = above_ground),
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
    colour = "Above ground"
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
       aes(x = d15n_korr, y = d13c_korr, colour = above_ground)) +
  # Non-outliers: above ground specific coloured circles
  geom_point(data = subset(outliers_iqr, !outlier),
             size = 3, alpha = 0.9) +
  # Outliers: above ground specific coloured x
  geom_point(data = subset(outliers_iqr, outlier),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  # Label only outliers (toggle to label all)
  geom_text(
    data = subset(outliers_iqr, outlier),
    aes(label = nr),
    nudge_y = 0.3, size = 5, show.legend = FALSE, colour = "black"
  ) +
  labs(
    title = "Scatterplot of unlabelled outliers δ13C vs δ15N                      X = Tukey's fences outliers",
    x = "δ15N",
    y = "δ13C",
    colour = "Above ground"
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
       aes(x = above_ground, y = d13c_korr, colour = beriget)) +
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
       aes(x = above_ground, y = d15n_korr, colour = beriget)) +
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


# remove outliers from natabun and irms_names

remove_outlier <- c(176, 158, 159, 163, 157)
natabun_clean <- natabun %>% filter(!nr %in% remove_outlier)
irms_above_ground_names <- irms_above_ground_names %>% filter(!nr %in% remove_outlier)

# ---------- correct above ground specific natural abundances ----------

# create nat_abundance averages from run 1 AND Cut treatment
avg_table_natabun <- natabun_clean %>%
  group_by(above_ground) %>%
  summarise(
    avg_d15  = mean(d15n_korr, na.rm = TRUE),
    n_d15    = sum(!is.na(d15n_korr)),
    avg_d13  = mean(d13c_korr, na.rm = TRUE),
    n_d13    = sum(!is.na(d13c_korr)),
    n_both   = sum(!is.na(d15n_korr) & !is.na(d13c_korr)),
    n_rows   = n(),
    .groups  = "drop"
  )

avg_table_natabun

# Replace nat abundance atm pct with above ground–specific values

#      atom%15N = 100 * R15 * ( (δ/1000 + 1) / (1 + R15 * (δ/1000 + 1)) ), where R15 = 0.003676
#      atom%13C = 100 * R13 * ( (δ/1000 + 1) / (1 + R13 * (δ/1000 + 1)) ), where R13 = 0.011237

irms_above_ground_natabund_corrected <- irms_above_ground_names %>%
  left_join(avg_table_natabun, by = c("above_ground")) %>%
  mutate(
    nat_abun_15n_atm_pct = 100 * 0.003676 * ((avg_d15 / 1000) + 1) /
      (1 + 0.003676 * ((avg_d15 / 1000) + 1)),
    nat_abun_13c_atm_pct = 100 * 0.011237 * ((avg_d13 / 1000) + 1) /
      (1 + 0.011237 * ((avg_d13 / 1000) + 1))
  )

# ---------- correct atom % ----------

irms_above_ground_atompct_corrected <- irms_above_ground_natabund_corrected %>% 
  mutate(
    atom_pct_15n = 100 * 0.003676 * (d15n_korr / 1000 + 1)/(1 + 0.003676 * (d15n_korr / 1000 + 1)), 
    atom_pct_13c = 100 * 0.011237 * (d13c_korr / 1000 + 1)/(1 + 0.011237 * (d13c_korr / 1000 + 1))
  )

# ---------- correct APE ----------

irms_above_ground_ape_corrected <- irms_above_ground_atompct_corrected %>%
  mutate(
    ape_pct_15n = atom_pct_15n - nat_abun_15n_atm_pct,
    ape_pct_13c = atom_pct_13c - nat_abun_13c_atm_pct
  )

# ---------- correct full dataframe ----------

irms_above_ground_corrected <- irms_above_ground_ape_corrected %>%
  mutate(
    # correct 
    n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000,                   #15N pr DW
    n15_per_n_ug_per_g  = (((mg_n           * (ape_pct_15n / 100)) / mg_n) * 1000) * 1000,  #15N pr N
    c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,                   #13C pr DW
    c13_per_c_ug_per_g  = (((mg_c           * (ape_pct_13c / 100)) / mg_c) * 1000) * 1000,  #13C pr C
  )

#library(writexl)
#write_xlsx(irms_above_ground_corrected, "./data_files/IRMS/IRMS_above_ground_corrected.xlsx")

# ---------- check for outliers ----------

irms_outliers_iqr <- irms_above_ground_corrected %>%
  group_by(above_ground) %>%
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

# ---------- plot outliers----------

ggplot(irms_outliers_iqr,
       aes(x = d13c_korr, y = d15n_korr, colour = above_ground)) +
  # Non-outliers: above ground specific coloured circles
  geom_point(data = subset(irms_outliers_iqr, !outlier),
             size = 3, alpha = 0.9) +
  # Outliers: above ground specific coloured x
  geom_point(data = subset(irms_outliers_iqr, outlier_d13c),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  geom_point(data = subset(irms_outliers_iqr, outlier_d15n),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  # Label only outliers (toggle to label all)
  geom_text(
    data = subset(irms_outliers_iqr, outlier_d13c),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_d15n),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  labs(
    title = "Scatterplot of all dataset outliers for d13C and d15N                 X = Tukey's fences outliers",
    x = "d13C",
    y = "d15N",
    colour = "Above ground"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

ggplot(irms_outliers_iqr,
       aes(x = c13_per_dw_ug_per_g, y = n15_per_dw_ug_per_g, colour = above_ground)) +
  # Non-outliers: above ground specific coloured circles
  geom_point(data = subset(irms_outliers_iqr, !outlier),
             size = 3, alpha = 0.9) +
  # Outliers: above ground specific coloured x
  geom_point(data = subset(irms_outliers_iqr, outlier_c13_dw),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  geom_point(data = subset(irms_outliers_iqr, outlier_n15_dw),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  # Label only outliers (toggle to label all)
  geom_text(
    data = subset(irms_outliers_iqr, outlier_c13_dw),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_n15_dw),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  labs(
    title = "Scatterplot of all dataset outliers for 13C and 15N pr DW ug pr gr                 X = Tukey's fences outliers",
    x = "13C pr DW ug pr gr",
    y = "15N pr DW ug pr gr",
    colour = "Above ground"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

ggplot(irms_outliers_iqr,
       aes(x = n15_per_dw_ug_per_g, y = c13_per_dw_ug_per_g, colour = above_ground)) +
  # Non-outliers: above ground specific coloured circles
  geom_point(data = subset(irms_outliers_iqr, !outlier),
             size = 3, alpha = 0.9) +
  # Outliers: above ground specific coloured x
  geom_point(data = subset(irms_outliers_iqr, outlier_c13_dw),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  geom_point(data = subset(irms_outliers_iqr, outlier_n15_dw),
             aes(colour = above_ground),
             shape = 4, size = 3.5, stroke = 1.2) +
  # Label only outliers (toggle to label all)
  geom_text(
    data = subset(irms_outliers_iqr, outlier_c13_dw),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_n15_dw),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  labs(
    title = "Scatterplot of all dataset outliers for 13C and 15N pr DW ug pr gr                 X = Tukey's fences outliers",
    x = "15N pr DW ug pr gr",
    y = "13C pr DW ug pr gr",
    colour = "Above ground"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

# ---------- look at the data ----------

ggplot(irms_above_ground_corrected,
       aes(x = above_ground, y = d13c_korr, colour = beriget)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(irms_outliers_iqr, outlier_d13c),
    aes(colour = beriget),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_d13c),
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

ggplot(irms_above_ground_corrected,
       aes(x = above_ground, y = n15_per_dw_ug_per_g, colour = beriget)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(irms_outliers_iqr, outlier_n15_dw),
    aes(colour = beriget),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_n15_dw),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

ggplot(irms_above_ground_corrected,
       aes(x = above_ground, y = c13_per_dw_ug_per_g, colour = beriget)) +
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
    nudge_y = 1.7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )


ggplot(irms_above_ground_corrected,
       aes(x = run, y = c13_per_dw_ug_per_g, colour = above_ground)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(irms_outliers_iqr, outlier_c13_dw),
    aes(colour = above_ground),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_c13_dw),
    aes(label = nr),
    nudge_y = 1.7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

ggplot(irms_above_ground_corrected,
       aes(x = run, y = n15_per_dw_ug_per_g, colour = above_ground)) +
  theme_bw() + 
  geom_point(shape = 21, size = 3, stroke = .5) +
  geom_point(
    data = subset(irms_outliers_iqr, outlier_n15_dw),
    aes(colour = above_ground),
    shape = 4, size = 4, stroke = 1.5
  ) +
  geom_text(
    data = subset(irms_outliers_iqr, outlier_n15_dw),
    aes(label = nr),
    nudge_y = 7, size = 5, show.legend = FALSE, colour = "black"
  ) +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
