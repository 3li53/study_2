# function to check, install, and load packages
check_install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# load analysis packages
invisible(lapply(c("dplyr","stringr","car","emmeans","lme4","lmerTest","effectsize","broom","ggplot2","tidyr"),
                 check_install_load))

# Work on the full flux_data
# Exclude rows with missing flux or malformed treatment codes
flux_data3 <- flux_data %>%
  filter(!is.na(flux_mg_m2_hr), !is.na(treatment)) %>%
  mutate(
    vegetation_type = factor(str_sub(treatment, 1, 1), levels = c("B","S")),
    cut_treatment = factor(str_sub(treatment, 2, 2), levels = c("C","U")),
    wet_treatment = factor(str_sub(treatment, 3, 3), levels = c("D","W"))
  )

# Optional: give descriptive labels (replace with your actual treatment names)
levels(flux_data3$vegetation_type) <- c("Bryophyte", "Salix")
levels(flux_data3$cut_treatment) <- c("Cut", "Uncut")  
levels(flux_data3$wet_treatment) <- c("Dry", "Wet")

# Sanity check: counts per cell
flux_data3 %>% count(vegetation_type, cut_treatment, wet_treatment)

# Use Type III SS
options(contrasts = c("contr.sum", "contr.poly"))

# Base 3-way ANOVA on treatments only
lm3 <- lm(flux_mg_m2_hr ~ vegetation_type*cut_treatment*wet_treatment, data = flux_data3)
anova_type3 <- car::Anova(lm3, type = 3)
print(anova_type3)

# Effect sizes (partial eta-squared)
eta <- effectsize::eta_squared(anova_type3, partial = TRUE)
print(eta)

# Tidy table if you want a neat data frame
anova_tidy <- broom::tidy(anova_type3)
anova_tidy





#include hrs_after_incubation as fixed factor
lm3_time <- lm(flux_mg_m2_hr ~ vegetation_type*cut_treatment*wet_treatment + hrs_after_incubation, data = flux_data3)
anova_time <- car::Anova(lm3_time, type = 3)
anova_time
effectsize::eta_squared(anova_time, partial = TRUE)


# Mixed model: treatments + time as fixed; run as random intercept
lmer3 <- lmerTest::lmer(
  flux_mg_m2_hr ~ vegetation_type*cut_treatment*wet_treatment + hrs_after_incubation + (1 | run),
  data = flux_data3
)

summary(lmer3)

# Type III tests with Satterthwaite degrees of freedom
library(lmerTest)

# keep Type III SS and Satterthwaite df
options(contrasts = c("contr.sum", "contr.poly"))

lmer3 <- lmerTest::lmer(
  flux_mg_m2_hr ~ vegetation_type*cut_treatment*wet_treatment + hrs_after_incubation + (1 | run),
  data = flux_data3
)

# call the generic; lmerTest supplies the correct method
lmer3_anova <- anova(lmer3, type = 3, ddf = "Satterthwaite")
lmer3_anova


# Residual diagnostics for the fixed-effects model
par(mfrow = c(1, 2))
plot(lm3, which = 1)  # residuals vs fitted
plot(lm3, which = 2)  # QQ plot

# Global normality test (large n -> often overly sensitive; use plots + robust checks)
shapiro.test(residuals(lm3))

# Homogeneity of variance across the 2x2x2 cells
car::leveneTest(flux_mg_m2_hr ~ vegetation_type*cut_treatment*wet_treatment, data = flux_data3)
par(mfrow = c(1, 1))


# Mixed-model residual diagnostics
resid_df <- data.frame(
  fitted = fitted(lmer3),
  resid  = resid(lmer3)
)

ggplot(resid_df, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Residuals vs Fitted (mixed model)")


# For the fixed-effects model:
emm_cells <- emmeans::emmeans(lm3, ~ vegetation_type*cut_treatment*wet_treatment)
emm_cells

# Pairwise comparisons of all 8 cells (Tukey-adjusted)
pairs(emm_cells, adjust = "tukey")

# If a 2-way or 3-way interaction is significant, probe simple effects:
# Example: simple effects of F1 at each combination of F2 and F3
emm_veg_by_w_d <- emmeans::emmeans(lm3, ~ vegetation_type | cut_treatment*wet_treatment)
pairs(emm_veg_by_w_d, adjust = "tukey")

# Main-effect EMMs (if interaction not significant)
emm_veg <- emmeans::emmeans(lm3, ~ vegetation_type)
pairs(emm_veg, adjust = "tukey")

emm_c <- emmeans::emmeans(lm3, ~ cut_treatment)
pairs(emm_c, adjust = "tukey")

emm_w <- emmeans::emmeans(lm3, ~ wet_treatment)
pairs(emm_w, adjust = "tukey")


# main-effect EMMs for mixed model
emm_cells_mixed <- emmeans::emmeans(lmer3, ~ vegetation_type*cut_treatment*wet_treatment)
pairs(emm_cells_mixed, adjust = "tukey")



# Build a tidy data frame of EMMs
emm_df <- as.data.frame(emm_cells) %>%
  rename(emmean = emmean, SE = SE)  # keep consistent names

# Plot interaction: F1 vs F2, faceted by F3
ggplot(emm_df, aes(x = cut_treatment, y = emmean, color = vegetation_type, group = vegetation_type)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.2, position = position_dodge(width = 0.2)) +
  facet_wrap(~ wet_treatment) +
  theme_minimal() +
  labs(
    title = "Estimated Marginal Means: 3-way Treatment Structure",
    x = "cut_treatment",
    y = "Estimated mean flux (mg CO₂ m⁻² h⁻¹)",
    color = "vegetation_type)"
  )
