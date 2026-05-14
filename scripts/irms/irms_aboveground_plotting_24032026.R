###############################################################
# IRMS ROOTS — PLOTTING SCRIPT
# Requires objects created in: irms_processing.R


library(tidyverse)
library(ggplot2)
library(ggbreak)
library(patchwork)

# ---------------------- HELPER PLOTS ------------------------------

plot_scatter_outliers <- function(df, x, y, outlier_col, colour = colour_vars) {
  if (!outlier_col %in% names(df)) {
    stop(paste("Outlier column", outlier_col, "not found"))
  }
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_point(size = 3, alpha = 0.7) +
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

colour_vars <- c(
  "aboveground"
)


plot_scatter_summary <- function(df,
                                 xmean, ymean,
                                 xmin, xmax,
                                 ymin, ymax,
                                 label_col,
                                 point_palette = c("grey30", "grey40", "grey50",
                                                   "grey60", "grey70", "grey80")) {
  ggplot(df,
         aes_string(x = xmean,
                    y = ymean,
                    colour = label_col)) +
    # points (greyscale)
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = point_palette) +
    # vertical error bars
    geom_errorbar(
      aes_string(ymin = ymin, ymax = ymax),
      width = 0.02
    ) +
    # horizontal error bars (orientation = "y")
    geom_errorbar(
      aes_string(xmin = xmin, xmax = xmax),
      height = 0.02,
      orientation = "y"
    ) +
    # labels just northeast of points
    geom_text(
      aes_string(label = label_col),
      hjust = -0.1,    # moves text to the right
      vjust = -0.6,    # moves text upward
      size = 4,
      colour = "black"
    ) +
    labs(x = expression("δ"^"13"*"C ‰"), y = expression("δ"^"15"*"N ‰"), title = ) +
    theme_classic(base_size = 16) +
    guides(colour = "none")  # hide legend since labels show group
}



plot_box <- function(df, x, y, colour = "beriget", title = NULL) {
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    theme_minimal(base_size = 16) +
    labs(title = title, x = x, y = y)
}

# ---------------------- NATURAL 13C and 15N natural ABUNDANCE PLOTS ------------------------------

plot_scatter_outliers(natabun_outliers, "d13c_korr", "d15n_korr", "outlier")
plot_scatter_outliers(natabun_outliers, "aboveground", "d15n_korr", "outlier_d15n_korr")
plot_scatter_outliers(natabun_outliers, "aboveground", "d13c_korr", "outlier_d13c_korr")
plot_scatter_outliers(biomass_outliers, "aboveground", "component_weight", "outlier")

plot_scatter_summary(
  natabun_means,
  xmean     = "avg_d13",
  ymean     = "avg_d15",
  xmin      = "xmin",
  xmax      = "xmax",
  ymin      = "ymin",
  ymax      = "ymax",
  label_col = "aboveground")

# ---------------------- CORRECTED DATA PLOTS ------------------------------

plot_scatter_outliers(irms_outliers, "c13_per_dw_ug_per_g", "n15_per_dw_ug_per_g", "outlier")
plot_scatter_outliers(irms_outliers, "aboveground", "c13_per_dw_ug_per_g", "outlier_c13_per_dw_ug_per_g")
plot_scatter_outliers(irms_outliers, "aboveground", "n15_per_dw_ug_per_g", "outlier_n15_per_dw_ug_per_g")

plot_box(irms, "aboveground", "n15_per_dw_ug_per_g")
plot_box(irms, "aboveground", "c13_per_dw_ug_per_g")
plot_box(irms, "aboveground", "n_per_dw_mg_per_g")

# After removing outliers
plot_scatter_outliers(irms, "aboveground", "c13_per_dw_ug_per_g", "outlier_c13_per_dw_ug_per_g")
plot_scatter_outliers(irms, "aboveground", "n15_per_dw_ug_per_g", "outlier_n15_per_dw_ug_per_g")

#

# ---------------------- SUMMARY DF concentrations by aboveground ------------------------------
#this setion makes summaries for plotting later in the plotting script

irms_summary <- irms %>%
  group_by(aboveground, beriget) %>%
  summarise(
    mean_n15_per_dw_ug_per_g = mean(n15_per_dw_ug_per_g, na.rm = TRUE),
    se_n15_per_dw_ug_per_g   = sd(n15_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    mean_c13_per_dw_ug_per_g = mean(c13_per_dw_ug_per_g, na.rm = TRUE),
    se_c13_per_dw_ug_per_g   = sd(c13_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    mean_n_per_dw_mg_per_g   = mean(n_per_dw_mg_per_g, na.rm = TRUE),
    se_n_per_dw_mg_per_g     = sd(n_per_dw_mg_per_g, na.rm = TRUE) / sqrt(n()), 
    .groups = "drop"
  )

irms_summary_N <- irms %>%
  group_by(aboveground) %>%
  summarise(
    mean_n_per_dw_mg_per_g   = mean(n_per_dw_mg_per_g, na.rm = TRUE),
    se_n_per_dw_mg_per_g     = sd(n_per_dw_mg_per_g, na.rm = TRUE) / sqrt(n()), 
    .groups = "drop"
  )

irms_plot <- irms %>%
  mutate(n15_per_dw_ug_per_g = pmax(n15_per_dw_ug_per_g, 0))   # used only for plotting

irms_summary_plot <- irms_plot %>%
  group_by(aboveground, beriget) %>%
  summarise(
    mean_n15_per_dw_ug_per_g = mean(n15_per_dw_ug_per_g, na.rm = TRUE),
    se_n15_per_dw_ug_per_g   = sd(n15_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    mean_c13_per_dw_ug_per_g = mean(c13_per_dw_ug_per_g, na.rm = TRUE),
    se_c13_per_dw_ug_per_g   = sd(c13_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


# ---------------------- POOL SUMMARY ------------------------------

# 
pool_summary <- aboveground %>%
  group_by(aboveground, beriget) %>%
  summarise(
    mean_15N_ug = mean(pulje_15N_ug, na.rm = TRUE),
    se_15N_ug   = sd(pulje_15N_ug, na.rm = TRUE) / sqrt(n()),
    mean_13C_ug = mean(pulje_13C_ug, na.rm = TRUE),
    se_13C_ug   = sd(pulje_13C_ug, na.rm = TRUE) / sqrt(n()),
    
    mean_15N_mg = mean(pulje_15N_mg, na.rm = TRUE),
    se_15N_mg   = sd(pulje_15N_mg, na.rm = TRUE) / sqrt(n()),
    mean_13C_mg = mean(pulje_13C_mg, na.rm = TRUE),
    se_13C_mg   = sd(pulje_13C_mg, na.rm = TRUE) / sqrt(n()),
    
    mean_15N_mg_m2 = mean(pulje_15N_mg_m2, na.rm = TRUE),
    se_15N_mg_m2   = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
    mean_13C_mg_m2 = mean(pulje_13C_mg_m2, na.rm = TRUE),
    se_13C_mg_m2   = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
    
    mean_N_mg_m2 = mean(pulje_N_mg_m2, na.rm = TRUE),
    se_N_mg_m2   = sd(pulje_N_mg_m2, na.rm = TRUE) / sqrt(n()),
    
    .groups = "drop"
  )

aboveground_plot <- aboveground %>%
  mutate(pulje_15N_ug = pmax(pulje_15N_ug, 0),
         pulje_15N_mg = pmax(pulje_15N_mg, 0)
  )   # used only for plotting

pool_summary_plot <- aboveground_plot %>%
  group_by(aboveground, beriget) %>%
  summarise(
    mean_15N_ug = mean(pulje_15N_ug, na.rm = TRUE),
    se_15N_ug   = sd(pulje_15N_ug, na.rm = TRUE) / sqrt(n()),
    mean_13C_ug = mean(pulje_13C_ug, na.rm = TRUE),
    se_13C_ug   = sd(pulje_13C_ug, na.rm = TRUE) / sqrt(n()),
    
    mean_15N_mg = mean(pulje_15N_mg, na.rm = TRUE),
    se_15N_mg   = sd(pulje_15N_mg, na.rm = TRUE) / sqrt(n()),
    mean_13C_mg = mean(pulje_13C_mg, na.rm = TRUE),
    se_13C_mg   = sd(pulje_13C_mg, na.rm = TRUE) / sqrt(n()),
    
    mean_15N_mg_m2 = mean(pulje_15N_mg_m2, na.rm = TRUE),
    se_15N_mg_m2   = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
    mean_13C_mg_m2 = mean(pulje_13C_mg_m2, na.rm = TRUE),
    se_13C_mg_m2   = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
    
    mean_N_mg_m2 = mean(pulje_N_mg_m2, na.rm = TRUE),
    se_N_mg_m2   = sd(pulje_N_mg_m2, na.rm = TRUE) / sqrt(n()),
    
    .groups = "drop"
  )

pool_summary_plot_N <- aboveground_plot %>%
  group_by(aboveground) %>%
  summarise(
    mean_N_mg_m2 = mean(pulje_N_mg_m2, na.rm = TRUE),
    se_N_mg_m2   = sd(pulje_N_mg_m2, na.rm = TRUE) / sqrt(n()),
    
    .groups = "drop"
  )

# ---------------------- TOTAL POOLS PER PIPE (Beriget) ------------------------------

# isotopes, beriget
total_aboveground_pools <- aboveground_plot %>%
  group_by(veg, cut, wet, treatment, beriget) %>%
  summarise(
    total_15N_ug = sum(pulje_15N_ug, na.rm = TRUE),
    se_15N_ug    = sd(pulje_15N_ug, na.rm = TRUE) / sqrt(n()),
    total_13C_ug = sum(pulje_13C_ug, na.rm = TRUE),
    se_13C_ug    = sd(pulje_13C_ug, na.rm = TRUE) / sqrt(n()),
    
    total_15N_mg = sum(pulje_15N_mg, na.rm = TRUE),
    se_15N_mg    = sd(pulje_15N_mg, na.rm = TRUE) / sqrt(n()),
    total_13C_mg = sum(pulje_13C_mg, na.rm = TRUE),
    se_13C_mg    = sd(pulje_13C_mg, na.rm = TRUE) / sqrt(n()),
    
    total_15N_mg_m2 = sum(pulje_15N_mg_m2, na.rm = TRUE),
    se_15N_mg_m2    = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
    total_13C_mg_m2 = sum(pulje_13C_mg_m2, na.rm = TRUE),
    se_13C_mg_m2    = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --------------------- total POOLS PER PIPE & ABOVEGROUND --------------------------

# isotopes, beriget
total_aboveground_pools_stacked <- aboveground %>% 
  group_by(veg, cut, wet, treatment, beriget, aboveground) %>% 
  summarise(
    total_15N_ug = sum(pulje_15N_ug, na.rm = TRUE),
    se_15N_ug    = sd(pulje_15N_ug, na.rm = TRUE) / sqrt(n()),
    total_13C_ug = sum(pulje_13C_ug, na.rm = TRUE),
    se_13C_ug    = sd(pulje_13C_ug, na.rm = TRUE) / sqrt(n()),
    
    total_15N_mg = sum(pulje_15N_mg, na.rm = TRUE),
    se_15N_mg    = sd(pulje_15N_mg, na.rm = TRUE) / sqrt(n()),
    total_13C_mg = sum(pulje_13C_mg, na.rm = TRUE),
    se_13C_mg    = sd(pulje_13C_mg, na.rm = TRUE) / sqrt(n()),
    
    total_15N_mg_m2 = sum(pulje_15N_mg_m2, na.rm = TRUE),
    se_15N_mg_m2    = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
    total_13C_mg_m2 = sum(pulje_13C_mg_m2, na.rm = TRUE),
    se_13C_mg_m2    = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

total_aboveground_pools_stacked_plot <- subset(aboveground_plot, beriget == TRUE) %>% 
  group_by(veg, cut, wet, treatment, beriget, aboveground) %>% 
  summarise(
    total_15N_ug = sum(pulje_15N_ug, na.rm = TRUE),
    se_15N_ug    = sd(pulje_15N_ug, na.rm = TRUE) / sqrt(n()),
    total_13C_ug = sum(pulje_13C_ug, na.rm = TRUE),
    se_13C_ug    = sd(pulje_13C_ug, na.rm = TRUE) / sqrt(n()),
    
    total_15N_mg = sum(pulje_15N_mg, na.rm = TRUE),
    se_15N_mg    = sd(pulje_15N_mg, na.rm = TRUE) / sqrt(n()),
    total_13C_mg = sum(pulje_13C_mg, na.rm = TRUE),
    se_13C_mg    = sd(pulje_13C_mg, na.rm = TRUE) / sqrt(n()),
    
    total_15N_mg_m2 = sum(pulje_15N_mg_m2, na.rm = TRUE),
    se_15N_mg_m2    = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
    total_13C_mg_m2 = sum(pulje_13C_mg_m2, na.rm = TRUE),
    se_13C_mg_m2    = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

total_aboveground_pools_stacked_N <- aboveground %>% 
  group_by(veg, cut, wet, treatment, aboveground) %>% 
  summarise(
    total_N_mg = sum(pulje_N_mg_m2, na.rm = TRUE),
    se_N_mg = sd(pulje_N_mg_m2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# END OF summaries
# Objects produced above are used in the plotting below
###########################################################

# ---------------------- NITROGEN CONCENTRATIONS INCL NATABUN ---------------------------

p1 <- ggplot(subset(irms_summary_N),
             aes(x = aboveground, y = mean_n_per_dw_mg_per_g)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_n_per_dw_mg_per_g,
        ymax = mean_n_per_dw_mg_per_g + se_n_per_dw_mg_per_g),
    width = 0.2
  ) +
 # scale_fill_manual(values = "grey50", guide = "none") +
  scale_x_discrete(
    # limits = "bryophyte", "equisetum", .... # order
    labels = c("bryophyte" = "Bryophytes",
               "equisetum" = "Equisetum",
               "graminoid" = "Graminoids",
               "lichen" = "Lichen",
               "salix_leaf" = "Salix leaves",
               "salix_stem" = "Salix stems"
    ),
    expand = expansion(mult = c(0.15, 0.1))) +
  labs(x = "", y = "N per DW, mg/g", title = "N concentration in Vegetation") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )
p1

# ---------------------- NITROGEN POOLS INCL NATABUN ------------------------------------

p1 <- ggplot(subset(pool_summary_plot_N),
             aes(x = aboveground, y = mean_N_mg_m2)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_N_mg_m2,
        ymax = mean_N_mg_m2 + se_N_mg_m2),
    width = 0.2
  ) +
 # scale_fill_manual(values = "grey50", guide = "none") +
  scale_x_discrete(
    # limits = "bryophyte", "equisetum", .... # order
    labels = c("bryophyte" = "Bryophytes",
               "equisetum" = "Equisetum",
               "graminoid" = "Graminoids",
               "lichen" = "Lichen",
               "salix_leaf" = "Salix leaves",
               "salix_stem" = "Salix stems"
    ),
    expand = expansion(mult = c(0.15, 0.1))) +
  labs(x = "", y = "N mg m-2", title = "N pools in above ground vegetation") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )
p1

# ---------------------- ENRICHED concentration BARPLOTS --------------------------------

# 13C & 15N pr DW (ug per g)

p1 <- ggplot(subset(irms_summary_plot, beriget == TRUE),
             aes(x = aboveground, y = mean_c13_per_dw_ug_per_g, fill = beriget)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_c13_per_dw_ug_per_g,
        ymax = mean_c13_per_dw_ug_per_g + se_c13_per_dw_ug_per_g),
    width = 0.2
  ) +
  scale_fill_manual(values = "grey50", guide = "none") +
  scale_x_discrete(
    # limits = "bryophyte", "equisetum", .... # order
    labels = c("bryophyte" = "Bryophytes",
               "equisetum" = "Equisetum",
               "graminoid" = "Graminoids",
               "lichen" = "Lichen",
               "salix_leaf" = "Salix leaves",
               "salix_stem" = "Salix stems"
    ),
    expand = expansion(mult = c(0.15, 0.1))) +
  labs(x = "", y = "13C per DW, µg/g", title = "Enriched 13C and 15N in Vegetation") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )

p2 <- ggplot(subset(irms_summary_plot, beriget == TRUE),
             aes(x = aboveground, y = mean_n15_per_dw_ug_per_g, fill = beriget)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_n15_per_dw_ug_per_g,
        ymax = mean_n15_per_dw_ug_per_g + se_n15_per_dw_ug_per_g),
    width = 0.2
  ) +
  scale_fill_manual(values = "grey30", guide = "none") +
  scale_x_discrete(
    # limits = "bryophyte", "equisetum", .... # order
    labels = c("bryophyte" = "Bryophytes",
               "equisetum" = "Equisetum",
               "graminoid" = "Graminoids",
               "lichen" = "Lichen",
               "salix_leaf" = "Salix leaves",
               "salix_stem" = "Salix stems"
               ),
    expand = expansion(mult = c(0.15, 0.1))) +
  labs(x = "", y = "15N per DW, µg/g") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )
  
p1 | p2


# Combined barplot
ggplot(irms_summary_plot, aes(x = aboveground, y = mean_n15_per_dw_ug_per_g, fill = beriget)) +
  geom_col(colour = "black", position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = mean_n15_per_dw_ug_per_g - se_n15_per_dw_ug_per_g,
        ymax = mean_n15_per_dw_ug_per_g + se_n15_per_dw_ug_per_g),
    position = position_dodge(0.9), width = 0.2
  ) +
  theme_classic(base_size = 20) +
  scale_fill_manual(
    name = "Stable isotopes",
    values = c("FALSE" = "grey100","TRUE"  = "grey60"),
    labels = c("Natural abundance", "Enriched")
  ) +
  labs(x = "Vegetation", y = "15N pr DW (ug per g)") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))


# ---------------------- ENRICHED POOLS PLOTS ------------------------------


p1 <- ggplot(subset(pool_summary_plot, beriget == TRUE),
             aes(x = aboveground, y = mean_13C_mg_m2, fill = beriget)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_13C_mg_m2,
        ymax = mean_13C_mg_m2 + se_13C_mg_m2),
    width = 0.2
  ) +
  scale_fill_manual(values = "grey50", guide = "none") +
  scale_x_discrete(
    # limits = "bryophyte", "equisetum", .... # order
    labels = c("bryophyte" = "Bryophytes",
               "equisetum" = "Equisetum",
               "graminoid" = "Graminoids",
               "lichen" = "Lichen",
               "salix_leaf" = "Salix leaves",
               "salix_stem" = "Salix stems"
    ),
    expand = expansion(mult = c(0.15, 0.1))) +
  labs(x = "", y = "13C mg m-2", title = "Enriched 13C and 15N pools in above ground vegetation") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )

p2 <- ggplot(subset(pool_summary_plot, beriget == TRUE),
             aes(x = aboveground, y = mean_15N_mg_m2, fill = beriget)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_15N_mg_m2,
        ymax = mean_15N_mg_m2 + se_15N_mg_m2),
    width = 0.2
  ) +
  scale_fill_manual(values = "grey30", guide = "none") +
  scale_x_discrete(
  # limits = "bryophyte", "equisetum", .... # order
    labels = c("bryophyte" = "Bryophytes",
               "equisetum" = "Equisetum",
               "graminoid" = "Graminoids",
               "lichen" = "Lichen",
               "salix_leaf" = "Salix leaves",
               "salix_stem" = "Salix stems"
    ),
    expand = expansion(mult = c(0.15, 0.1)))+
  labs(x = "", y = "15N mg m-2") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.margin = margin(t = 5, r = 45, b = 0, l = 5)
  )


p1 | p2

# ---------------------- TOTAL STACKED POOLS PER PIPE ISOTOPES ------------------------------

p1 <- ggplot(total_aboveground_pools_stacked_plot,
             aes(x = treatment, y = total_13C_mg_m2, fill = aboveground)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    data = mean_aboveground_pools,
    aes(x = treatment,
        ymin = total_13C_mg_m2 + se_13C_mg_m2,
        ymax = total_13C_mg_m2 + se_13C_mg_m2),
    width = 0.05,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = mean_aboveground_pools,
    aes(x = treatment,
        ymin = total_13C_mg_m2,
        ymax = total_13C_mg_m2 + se_13C_mg_m2),
    width = 0.00,
    inherit.aes = FALSE,
  ) +
  scale_fill_manual(
    values = c("grey30", "grey40", "grey50", "grey60", "grey70", "grey80"), guide = "none"
    ) +
  scale_x_discrete(
    expand = expansion(mult = c(.17, 0.1)))+
  labs(x = "", y = "13C mg m2", title = "Enriched Aboveground 15N and 13C pools") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.margin = margin(t = 5, r = 25, b = 0, l = 5),
        legend.position = c(0.85, 0.85)
  )

p2 <- ggplot(total_aboveground_pools_stacked_plot,
             aes(x = treatment, y = total_15N_mg_m2, fill = aboveground)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    data = mean_aboveground_pools,
    aes(x = treatment,
        ymin = total_15N_mg_m2 + se_15N_mg_m2,
        ymax = total_15N_mg_m2 + se_15N_mg_m2),
        width = 0.05,
        inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = mean_aboveground_pools,
    aes(x = treatment,
        ymin = total_15N_mg_m2,
        ymax = total_15N_mg_m2 + se_15N_mg_m2),
    width = 0.00,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("grey30", "grey40", "grey50", "grey60", "grey70", "grey80"),
    name = "Vegetation",
    labels = c(
      bryophyte   = "Bryophytes",
      equisetum   = "Equisetum",
      graminoid   = "Graminoids",
      salix_leaf  = "Salix leaves",
      salix_stem  = "Salix stems",
      lichen      = "Lichen"
      )) +
  scale_x_discrete(
    expand = expansion(mult = c(.17, 0.1)))+
  labs(x = "", y = "15N mg m-2", title = "") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
        legend.position = c(0.80, 0.90)
        )
p1|p2

# total stacked pool N

p1 <- ggplot(total_aboveground_pools_stacked_N,
             aes(x = treatment, y = total_N_mg_m2, fill = aboveground)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    data = mean_aboveground_pools,
    aes(x = treatment,
        ymin = total_N_mg_m2 + se_N_mg_m2,
        ymax = total_N_mg_m2 + se_N_mg_m2),
    width = 0.05,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = mean_aboveground_pools,
    aes(x = treatment,
        ymin = total_N_mg_m2,
        ymax = total_N_mg_m2 + se_N_mg_m2),
    width = 0.00,
    inherit.aes = FALSE,
  ) +
  scale_fill_manual(
    values = c("grey30", "grey40", "grey50", "grey60", "grey70", "grey80"), guide = "none"
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(.17, 0.1)))+
  labs(x = "", y = "N mg m2", title = "Aboveground N pools") +
  scale_y_continuous(expand = expansion(mult = c(0.03,0.1))) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.margin = margin(t = 5, r = 25, b = 0, l = 5),
        legend.position = c(0.85, 0.85)
  )
p1
# ---------------- 15N  and 13C recovery % added N and C --------------------

