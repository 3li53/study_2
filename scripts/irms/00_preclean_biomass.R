#### ---- 00 preprocessing biomass script ----

### ---- roots ----

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

### ---- vegetation ----

aboveground_biomass_raw <- read.csv(
  "./data/raw/biomass/pipes_data_aboveground.csv",
  header = TRUE,
  skip = 1
)

aboveground_biomass_raw <- aboveground_biomass_raw %>%
  filter(!is.na(run)) %>%
  select(-aboveground_weight_tot_g) %>%
  mutate(
    run = as.character(run),
    veg_type = as.character(veg_type),
    cut_treat = as.character(cut_treat),
    wet_treat = as.character(wet_treat)
  )

aboveground_biomass <- aboveground_biomass_raw %>%
  mutate(
    veg = veg_type,
    cut = cut_treat,
    wet = wet_treat,
    beriget = !(run == "1" | cut == "C"),
    treatment = paste0(veg, cut, wet),
    aboveground_weight_tot_g = aboveground_weight_tot_g.1 - bag_g.1,
    weight_sorted_excl_bag_g = weight_sorted_excl.bag_g,
  ) %>%
  select(
    pipe.nr, run, veg, cut, wet, treatment, beriget, 
    aboveground_weight_tot_g, weight_sorted_excl_bag_g,
    equ_weight, graminoid_weight, bryophyte_weight,
    lichen_weight, stem_weight, leaves_weight
  )


library(readr)
write_csv(aboveground_biomass, "./data/raw/biomass/vegetation.csv")