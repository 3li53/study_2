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

### ---- soil ----

### ---- vegetation ----

