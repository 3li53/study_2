#### ---- 00 preprocessing biomass script ----

### ---- roots ----

root_biomass_raw <- read.csv(
  "./data/raw/biomass/pipes_data_roots.csv",   # import raw csv without headers
  header = FALSE,                              # first rows contain metadata, not column names
  stringsAsFactors = FALSE,                    # keep all text as character (no factors)
  check.names = FALSE                          # preserve original column names (no auto-fixing)
)

h_pos  <- root_biomass_raw[1, ] |> as.character()  # row 1: vertical position (mid / top / bottom)
h_frac <- root_biomass_raw[2, ] |> as.character()  # row 2: root fraction (fine / coarse)
h_var  <- root_biomass_raw[3, ] |> as.character()  # row 3: variable names (pipe nr / bag / fresh / dry)

fill_right <- function(x) {                         # helper: propagate last non-empty value to the right
  for (i in 2:length(x)) {
    if (is.na(x[i]) || x[i] == "") {               # if missing or empty
      x[i] <- x[i - 1]                             # fill with previous value
    }
  }
  x                                                # return filled vector
}

h_pos  <- fill_right(h_pos)                        # fill missing position labels across columns
h_frac <- fill_right(h_frac)                       # fill missing fraction labels across columns

new_names <- character(length(h_var))              # pre-allocate vector for new column names

for (i in seq_along(h_var)) {
  if (i <= 4) {                                    # first 4 columns = metadata (keep original names)
    new_names[i] <- h_var[i]
  } else {
    new_names[i] <- paste(h_pos[i], h_frac[i], h_var[i], sep = "_")  # construct hierarchical names
  }
}

root_biomass_raw <- root_biomass_raw[-c(1:3), ]    # remove header rows used to build column names
colnames(root_biomass_raw) <- new_names            # assign cleaned and combined column names

root_biomass_raw[ , -c(1:5)] <-
  lapply(root_biomass_raw[ , -c(1:5)], as.numeric) # convert measurement columns to numeric

# calculate root biomass
root_biomass <- root_biomass_raw %>% 
  mutate( #weight = dry sample        - bag weight
    mid_fine      = mid_fine_dry      - mid_fine_bag,
    mid_coarse    = mid_coarse_dry    - mid_coarse_bag,
    top_fine      = top_fine_dry      - top_fine_bag,
    top_coarse    = top_coarse_dry    - top_coarse_bag,
    bottom_fine   = bottom_fine_dry   - bottom_fine_bag,
    bottom_coarse = bottom_coarse_dry - bottom_coarse_bag
  ) %>% 
  select(-matches("fresh|dry|bag")) # drop fresh, dry and bag weights

write_csv(root_biomass, "./data/raw/biomass/roots.csv")  # save cleaned dataset


### ---- vegetation ----

aboveground_biomass_raw <- read.csv(
  "./data/raw/biomass/pipes_data_aboveground.csv",  # import aboveground biomass csv
  header = TRUE,                                     # use first row as column names
  skip = 1                                           # skip extra header/description row
)

aboveground_biomass_raw <- aboveground_biomass_raw %>%
  filter(!is.na(run)) %>%                            # remove rows without run identifier
  mutate(
    run = as.character(run),                         # ensure run is treated as categorical
    veg_type = as.character(veg_type),               # ensure vegetation type is character
    cut_treat = as.character(cut_treat),             # ensure cutting treatment is character
    wet_treat = as.character(wet_treat)              # ensure moisture treatment is character
  )

aboveground_biomass <- aboveground_biomass_raw %>%
  mutate(
    veg = veg_type,                                  # create shorthand variable for vegetation type
    cut = cut_treat,                                 # create shorthand variable for cutting treatment
    wet = wet_treat,                                 # create shorthand variable for moisture treatment
    beriget = !(run == "1" | cut == "C"),            # logical flag for unlabelled (not run 1 or cut C)
    treatment = paste0(veg, cut, wet),               # combined treatment identifier
    aboveground_weight_tot_g = aboveground_weight_tot_g - bag_g,  # subtract bag weight from total
    weight_sorted_excl_bag_g = weight_sorted_excl.bag_g,          # rename sorted weight variable
  ) %>%
  select(
    pipe.nr, run, veg, cut, wet, treatment, beriget, # retain identifiers and treatment variables
    aboveground_weight_tot_g, weight_sorted_excl_bag_g, # retain processed biomass metrics
    equ_weight, graminoid_weight, bryophyte_weight,     # retain functional group weights
    lichen_weight, stem_weight, leaves_weight, vascular_weight_g  # retain additional component weights
  )




library(readr)
write_csv(aboveground_biomass, "./data/raw/biomass/vegetation.csv")

# remove from global environment to declutter
rm(new_names, h_frac, h_pos, h_var)
