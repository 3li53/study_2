# load packages and data
install.packages("readxl")
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)

library(readxl)

dat <- c("run1", "run2","run3", "run4", "run5", "run6") %>%
  map_dfr(
    ~ read_excel(
  "EGM_data_for_fluxes.xlsx",
  sheet= .x)
)


dat <- dat %>%
  mutate(
    time_min = ifelse(
      grepl("^ppm_[0-9]+$", minute),
      as.numeric(sub("ppm_", "", minute)),
      NA_real_
    )
  ) %>%
  filter(!is.na(time_min))


# Pa to hPa
dat <- dat %>%
  mutate(p.air = p.air * 100)



# group measurements and nest
dat_grp <- dat %>%
  group_by(run, date, hrs_after_incubation, pipe_nr) %>%
  filter(n() >= 4) %>%   # mindestens 4 Zeitpunkte
  nest()

#regression and flux
R <- 8.314462618  # J mol-1 K-1

results <- dat_grp %>%
  mutate(
    fit = map(data, ~ lm(ppm ~ time_min, data = .x)),
    glance = map(fit, broom::glance),
    tidy   = map(fit, broom::tidy),
    
    slope_ppm_min = map_dbl(
      tidy, ~ .x$estimate[.x$term == "time_min"]
    ),
    
    r_squared = map_dbl(glance, "r.squared"),
    
    volume = map_dbl(data, ~ mean(.x$chamber_volume)),
    area   = map_dbl(data, ~ mean(.x$area)),
    temp_K = map_dbl(data, ~ mean(.x$t.air) + 273.15),
    press  = map_dbl(data, ~ mean(.x$p.air)),
    
    flux = slope_ppm_min * 1e-6 *
      press / (R * temp_K) *
      volume / area / 60
  ) %>%
  select(-data, -fit, -tidy, -glance)

# export to excel
install.packages("writexl")
library(writexl)
write_xlsx(
  results,
  path = "flux_results_allruns.xlsx"
)
getwd()

# statistical analysis
library(readxl)

result <- read_excel("flux_results_allruns.xlsx")

library(lme4)
library(lmerTest)

m_final <- lmer(
  flux ~ hrs_after_incubation * vegetation_type +
    (1 | run) +
    (1 | run:pipe_nr),
  data = result
)

# interpretation
anova(m_final)
summary(m_final)

# change rate per treatment
library(emmeans)

emtrends(
  m_final,
  specs = "vegetation_type",
  var = "hrs_after_incubation"
)

# ausreißer finden
result |> 
  arrange(desc(flux)) |> 
  head(10)

#remove outliers
max_value <- max(result$flux, na.rm = TRUE)

result_clean <- result |> 
  filter(flux != max_value)

# visualisation
library(ggplot2)

ggplot(result_clean,
       aes(hrs_after_incubation, flux, color = vegetation_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ run) +
  theme_bw()

#check model assumptions
plot(m_final)
qqnorm(resid(m_final))
qqline(resid(m_final))
