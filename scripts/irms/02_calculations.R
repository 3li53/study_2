#### ---- 02 calculations for irms and biomass data ----

#### ---- 00 packages ----
load_or_install(c("", "", "", "", ""))

### ---- 01 import ----

extracts_clean   <- read.csv("./data/clean/irms/extracts.csv", header = TRUE) # import irms output data
roots_clean      <- read.csv("./data/clean/irms/roots.csv", header = TRUE)
soil_clean       <- read.csv("./data/clean/irms/soil.csv", header = TRUE)
vegetation_clean <- read.csv("./data/clean/irms/vegetation.csv")
biomass_roots_clean      <- read.csv("./data/clean/irms/biomass_roots.csv", header = TRUE)     # import biomass data
biomass_vegetation_clean <- read.csv("./data/clean/irms/biomass_vegetation.csv", header = TRUE)

dfs <- list(  # make a df list
  extracts   = extracts_clean,
  roots      = roots_clean,
  soil       = soil_clean,
  vegetation = vegetation_clean,
  biomass_roots      = biomass_roots_clean,
  biomass_vegetation = biomass_vegetation_clean
)

### ---- extracts ----

### ---- roots ----

# first get the natural abundance means grouped by diameter
natabun_means <- dfs$roots %>%
  filter(!beriget) %>%          # only use unlabelled samples
  calc_isotope_means(diameter)  # calculate isotopes 
natabun_means                   # print the object

# then apply new means to the df, along with downstream corrections
roots_korr <- apply_baseline_correction(
  df                = dfs$roots,
  natabun_means_df  = natabun_means,
  group_var         = diameter
)

# calculate root biomass
dfs$biomass_roots <- dfs$biomass_roots %>% 
  mutate( #weight = dry sample        - bag weight
    mid_fine      = mid_fine_dry      - mid_fine_bag,
    mid_coarse    = mid_coarse_dry    - mid_coarse_bag,
    top_fine      = top_fine_dry      - top_fine_bag,
    top_coarse    = top_coarse_dry    - top_coarse_bag,
    bottom_fine   = bottom_fine_dry   - bottom_fine_bag,
    bottom_coarse = bottom_coarse_dry - bottom_coarse_bag
  ) %>% 
  select(-matches("fresh|dry|bag")) # drop fresh, dry and bag weights

### ---- soil ----

### ---- vegetation ----

