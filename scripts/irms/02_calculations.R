#### ---- 02 calculations for irms and biomass data ----

#### ---- 00 packages ----
load_or_install(c("tidyverse"))

### ---- 01 import data ----
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

### ---- corrections ----
# excel corrections are needed to correct for natural abundance means of group variables

## root group variables = diameter. there are two groups: fine and coarse.

# first get the natural abundance means grouped by diameter
natabun_means_roots <- dfs$roots %>%
  filter(!beriget) %>%          # only use unlabelled samples
  calc_isotope_means(diameter)  # calculate isotopes 
natabun_means_roots             # print the object

# then apply new means to the df, along with downstream corrections
roots_korr <- apply_baseline_correction(
  df                = dfs$roots,           # apply to roots dataframe
  natabun_means_df  = natabun_means_roots, # which mean values to use
  group_var         = diameter             # group by root diameter
)
dfs$roots <- roots_korr # transfer the corrected dataframe into the dataframe list

## We now do the same for vegetation. Here the variable is "aboveground" this includes 6 groups: equisetum, bryophytes, lichen, graminoid, salix_stem, salix_leaf

# first get the natural abundance means grouped by aboveground
natabun_means_vegetation <- dfs$vegetation %>%
  filter(!beriget) %>%             # only use unlabelled samples
  calc_isotope_means(aboveground)  # calculate isotopes 
natabun_means_vegetation           # print the object

# then apply new means to the df, along with downstream corrections
vegetation_korr <- apply_baseline_correction(
  df                = dfs$vegetation,           # apply to vegetation dataframe
  natabun_means_df  = natabun_means_vegetation, # which mean values to use
  group_var         = aboveground               # group by aboveground
)
dfs$vegetation <- vegetation_korr # transfer the corrected dataframe into the dataframe list

# now the data for roots and vegetation are corrected according to the natural abundance values of distinct groups.

## biomass

# biomass data needs to be made ready for pools calculations. 



# calculate vegetation biomass



### ---- soil calculations ----

# relative recovery = (excess in pool) / (total tracer added)
str(dfs$soil)

### ---- vegetation calculations ---- 



