################################################################################
# PLOTTING SCRIPT
library("patchwork")

################### Nitrogen concentration incl natabun ########################
# by aboveground vegetation

N_conc_aboveground <- aboveground %>%
  group_by(aboveground) %>%
  summarise(
    mean_N = mean(n_per_dw_mg_per_g, na.rm = TRUE),
    se_N   = sd(n_per_dw_mg_per_g, na.rm = TRUE) / sqrt(n()),
    n      = n(),
    .groups = "drop"
  )

ggplot(N_conc_aboveground,
       aes(x = aboveground, y = mean_N)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_N - se_N,
        ymax = mean_N + se_N),
    width = 0.2
  ) +
  labs(x = "", y = "N per DW (mg g⁻¹)",
       title = "Nitrogen concentration in aboveground vegetation") +
  theme_classic(base_size = 20)

# by treatment (stacked)

N_conc_treatment <- aboveground %>%
  group_by(treatment, aboveground) %>%
  summarise(
    mean_N = mean(n_per_dw_mg_per_g, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  )  %>%  filter(n >= 3)

ggplot(N_conc_treatment,
       aes(x = treatment, y = mean_N, fill = aboveground)) +
  geom_col(colour = "black") +
  labs(x = "", y = "N per DW (mg g⁻¹)",
       title = "Nitrogen concentration by treatment") +
  theme_classic(base_size = 20)

###################### Nitrgoen pools incl natabun #############################
# by aboveground

N_pool_aboveground <- aboveground %>%
  group_by(aboveground) %>%
  summarise(
    mean_N = mean(pulje_N_mg_m2, na.rm = TRUE),
    se_N   = sd(pulje_N_mg_m2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(N_pool_aboveground,
       aes(x = aboveground, y = mean_N)) +
  geom_col(colour = "black") +
  geom_errorbar(
    aes(ymin = mean_N - se_N,
        ymax = mean_N + se_N),
    width = 0.2
  ) +
  labs(x = "", y = "N pool (mg m⁻²)",
       title = "Aboveground nitrogen pools") +
  theme_classic(base_size = 20)

# by treatment (stacked)

N_pool_treatment <- aboveground %>%
  group_by(treatment, aboveground) %>%
  summarise(
    mean_N = mean(pulje_N_mg_m2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) # %>% filter(n >= 3)

ggplot(N_pool_treatment,
       aes(x = treatment, y = mean_N, fill = aboveground)) +
  geom_col(colour = "black") +
  labs(x = "", y = "N pool (mg m⁻²)",
       title = "Aboveground nitrogen pools by treatment") +
  theme_classic(base_size = 20)

################ Enriched concentrations (13C and 15N) #########################

# by aboveground

enriched_conc_aboveground <- aboveground %>%
  filter(beriget) %>%
  group_by(aboveground) %>%
  summarise(
    mean_15N = mean(n15_per_dw_ug_per_g, na.rm = TRUE),
    se_15N   = sd(n15_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    mean_13C = mean(c13_per_dw_ug_per_g, na.rm = TRUE),
    se_13C   = sd(c13_per_dw_ug_per_g, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p1 <- ggplot(enriched_conc_aboveground,
             aes(x = aboveground, y = mean_13C)) +
  geom_col(colour = "black", fill = "grey50") +
  geom_errorbar(
    aes(ymin = mean_13C ,
        ymax = mean_13C + se_13C),
    width = 0.2
  ) +
  labs(x = "", y = "¹³C per DW (µg g⁻¹)",
       title = "Enriched ¹³C and ¹⁵N concentrations") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )

p2 <- ggplot(enriched_conc_aboveground,
             aes(x = aboveground, y = mean_15N)) +
  geom_col(colour = "black", fill = "grey30") +
  geom_errorbar(
    aes(ymin = mean_15N,
        ymax = mean_15N + se_15N),
    width = 0.2
  ) +
  labs(x = "", y = "¹⁵N per DW (µg g⁻¹)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                   margin = margin(t = 5, r = 100, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
        plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
  )

p1 | p2

# by treatment

enriched_conc_treatment <- aboveground %>%
  filter(beriget) %>%
  group_by(treatment, aboveground) %>%
  summarise(
    mean_15N = mean(n15_per_dw_ug_per_g, na.rm = TRUE),
    mean_13C = mean(c13_per_dw_ug_per_g, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>% filter(n >= 3)

 p1 <- ggplot(enriched_conc_treatment,
       aes(x = treatment, y = mean_15N, fill = aboveground)) +
  geom_col(colour = "black") +
  labs(x = "", y = "¹⁵N per DW (µg g⁻¹)",
       title = "Enriched ¹⁵N concentrations by treatment") +
  theme_classic(base_size = 20)

 p2 <- ggplot(enriched_conc_treatment,
             aes(x = treatment, y = mean_13C, fill = aboveground)) +
  geom_col(colour = "black") +
  labs(x = "", y = "¹³C per DW (µg g⁻¹)",
       title = "Enriched ¹³C concentrations by treatment") +
  theme_classic(base_size = 20)
 
 p1 | p2
 
 ######################### enriched pools 13C and 15N ##########################
 # by aboveground 
 
 enriched_pools_aboveground <- aboveground %>%
   filter(beriget) %>%
   group_by(aboveground) %>%
   summarise(
     mean_13C = mean(pulje_13C_mg_m2, na.rm = TRUE),
     se_13C   = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
     mean_15N = mean(pulje_15N_mg_m2, na.rm = TRUE),
     se_15N   = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
     n = n(),
     .groups = "drop"
   )
 
 p1 <- ggplot(enriched_pools_aboveground,
              aes(x = aboveground, y = mean_13C)) +
   geom_col(colour = "black", fill = "grey50") +
   geom_errorbar(
     aes(ymin = mean_13C - se_13C,
         ymax = mean_13C + se_13C),
     width = 0.2
   ) +
   labs(x = "", y = "¹³C pool (mg m⁻²)",
        title = "Enriched aboveground pools") +
   theme_classic(base_size = 20) +
   theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                    margin = margin(t = 5, r = 100, b = 0, l = 0)),
         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
         axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
         plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
   )
 p2 <- ggplot(enriched_pools_aboveground,
              aes(x = aboveground, y = mean_15N)) +
   geom_col(colour = "black", fill = "grey30") +
   geom_errorbar(
     aes(ymin = mean_15N - se_15N,
         ymax = mean_15N + se_15N),
     width = 0.2
   ) +
   labs(x = "", y = "¹⁵N pool (mg m⁻²)") +
   theme_classic(base_size = 20) +
   theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                    margin = margin(t = 5, r = 100, b = 0, l = 0)),
         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
         axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
         plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
   )
 p1 | p2

 # by treatment
 
 enriched_pools_treatment <- subset(aboveground, beriget == TRUE) %>% 
   group_by(veg, cut, wet, treatment, beriget, aboveground) %>% 
   summarise(
     mean_13C = mean(pulje_13C_mg_m2, na.rm = TRUE),
     se_13C   = sd(pulje_13C_mg_m2, na.rm = TRUE) / sqrt(n()),
     mean_15N = mean(pulje_15N_mg_m2, na.rm = TRUE),
     se_15N   = sd(pulje_15N_mg_m2, na.rm = TRUE) / sqrt(n()),
     n = n(),
     .groups = "drop"
   )
 
 p1 <- ggplot(enriched_pools_treatment,
        aes(x = treatment, y = mean_15N, fill = aboveground)) +
   geom_bar(stat = "identity") +
   labs(x = "", y = "¹⁵N pool (mg m⁻²)",
        title = "Enriched aboveground pools by treatment") +
   theme_classic(base_size = 20)
 p2 <- ggplot(enriched_pools_treatment,
              aes(x = treatment, y = mean_13C, fill = aboveground)) +
   geom_bar(stat = "identity") +
   labs(x = "", y = "¹³C pool (mg m⁻²)",
        title = "") +
   theme_classic(base_size = 20)
 p1 | p2
 
 ########################### recovery % ########################################

 recovery_aboveground <- aboveground %>%
   filter(beriget) %>%
   group_by(aboveground) %>%
   summarise(
     mean_13C = mean(recovery_13C_pct, na.rm = TRUE),
     se_13C   = sd(recovery_13C_pct, na.rm = TRUE) / sqrt(n()),
     mean_15N = mean(recovery_15N_pct, na.rm = TRUE),
     se_15N   = sd(recovery_15N_pct, na.rm = TRUE) / sqrt(n()),
     .groups = "drop"
   )
 
 p1 <- ggplot(recovery_aboveground,
              aes(x = aboveground, y = mean_13C)) +
   geom_col(colour = "black", fill = "grey50") +
   geom_errorbar(
     aes(ymin = mean_13C - se_13C,
         ymax = mean_13C + se_13C),
     width = 0.2
   ) +
   labs(x = "", y = "¹³C recovery %",
        title = "Enriched aboveground recovery %") +
   theme_classic(base_size = 20) +
   theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                    margin = margin(t = 5, r = 100, b = 0, l = 0)),
         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
         axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
         plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
   )
 p2 <- ggplot(recovery_aboveground,
              aes(x = aboveground, y = mean_15N)) +
   geom_col(colour = "black", fill = "grey30") +
   geom_errorbar(
     aes(ymin = mean_15N - se_15N,
         ymax = mean_15N + se_15N),
     width = 0.2
   ) +
   labs(x = "", y = "¹⁵N recovery %") +
   theme_classic(base_size = 20) +
   theme(axis.text.x = element_text(angle = -30, hjust = 0, 
                                    margin = margin(t = 5, r = 100, b = 0, l = 0)),
         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
         axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 20)),
         plot.margin = margin(t = 5, r = 50, b = 0, l = 5)
   )
 p1 | p2
 
 # by treatment
 
 recovery_treatment <- subset(aboveground, beriget == TRUE) %>% 
   group_by(veg, cut, wet, treatment, beriget, aboveground) %>% 
   summarise(
     mean_13C = mean(recovery_13C_pct, na.rm = TRUE),
     se_13C   = sd(recovery_13C_pct, na.rm = TRUE) / sqrt(n()),
     mean_15N = mean(recovery_15N_pct, na.rm = TRUE),
     se_15N   = sd(recovery_15N_pct, na.rm = TRUE) / sqrt(n()),
     .groups = "drop"
   )
 
 p1 <- ggplot(recovery_treatment,
              aes(x = treatment, y = mean_15N, fill = aboveground)) +
   geom_bar(stat = "identity") +
   labs(x = "", y = "¹⁵N recovery %",
        title = "Enriched aboveground ¹⁵N recovery % by treatment") +
   theme_classic(base_size = 20)
 p2 <- ggplot(recovery_treatment,
              aes(x = treatment, y = mean_13C, fill = aboveground)) +
   geom_bar(stat = "identity") +
   labs(x = "", y = "¹³C recovery %",
        title = "Enriched aboveground ¹⁵N recovery % by treatment") +
   theme_classic(base_size = 20)
 p1 | p2
 