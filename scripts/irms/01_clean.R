#### ---- 01 cleaning irms and biomass data ----

### ---- 00 packages ----
load_or_install(c("readr", "dplyr", "janitor", "purrr", "tibble", "plotly"))

### ---- 01 import ----
extracts_raw   <- read.csv("./data/raw/irms/extracts/irms_extracts_raw.csv", header = TRUE) # import irms output data
roots_raw      <- read.csv("./data/raw/irms/roots/irms_roots_raw.csv", header = TRUE)
soil_raw       <- read.csv("./data/raw/irms/soil/irms_soil_raw.csv", header = TRUE)
vegetation_raw <- read.csv("./data/raw/irms/vegetation/irms_vegetation_raw.csv", header = TRUE)
biomass_roots_raw      <- read.csv("./data/raw/biomass/roots.csv", header = TRUE)     # import biomass data
biomass_vegetation_raw <- read.csv("./data/raw/biomass/vegetation.csv", header = TRUE)

dfs <- list(  # make a df list
  extracts   = extracts_raw,
  roots      = roots_raw,
  soil       = soil_raw,
  vegetation = vegetation_raw,
  biomass_roots      = biomass_roots_raw,
  biomass_vegetation = biomass_vegetation_raw
)

### ---- 02 initial cleaning ----
## ---- A. clean names ----
# compare names across dfs
out <- compare_clean_names(dfs)            
dfs <- out$cleaned_data                    # isolate dfs new names from nested list

# remove unwanted columns from all dataframes
remove_cols <- c(                          # define unwanted columns
  "tray", "position", "placement",
  "gosha_s_id_nr", "height_n_a", "height_n_a_1", "name", 
  "pe_an_a", "pea15n", "pe_an_a_1", "pea13c",
  "gnsn_std_n_weight", "gnsn_std_c_weight", 
  "x1n_a_to_x_mg_n_std", "x1n_a_to_x_mg_c_std",
  "sample_number", "nr", "well",
  "d15nkorr_2", "d13ckorr_2", "weight"
) 
dfs <- remove_cols_all(dfs, remove_cols)   # remove from all dfs

# rename columns if present in dfs
dfs <- rename_many_if_present(             
  dfs,
  c(
    "ape"           = "ape_pct_15n",
    "ape_1"         = "ape_pct_13c",
    "x_c"           = "c_pct",
    "x_n"           = "n_pct",
    "c_n"           = "cn_ratio",
    "n15"           = "n15_raw",
    "x13c"          = "c13_raw",
    "x13c_doc_ug_g" = "c13_doc_ug_g",
    "x13c_pr_c"     = "c13_ug_pr_gc",
    "x13c_pr_dw"    = "c13_ug_pr_gdw",
    "x15n_dtn_ug_g" = "n15_dtn_ug_g",
    "x15n_pr_dw"    = "n15_ug_pr_gdw",
    "x15n_pr_n"     = "n15_ug_pr_gn",
    "x3_4_mg"       = "sample_weight",
    "x4_5mg"        = "sample_weight",
    "x6_8mg"        = "sample_weight",
    "nat_abun_15n"  = "natabun_15n_atm_pct",
    "nat_abun_13c"  = "natabun_13c_atm_pct",
    "atom_15n"      = "atom_pct_15n",
    "atom_13c"      = "atom_pct_13c",
    "c_pr_dw"       = "c_mg_pr_gdw",
    "n_pr_dw"       = "n_mg_pr_gdw",
    "x_id"          = "treatment",
    "x13c_in_c"     = "d13c_in_c",
    "x15n_in_n"     = "d15n_in_n" 
  )
)

# ensure all dfs have important cols
dfs <- add_from_lookup(                    # copy tiny tag col from extracts to roots, soil and biomass_veg
  dfs = dfs,
  source_df = dfs$extracts,
  new_col = "tt",
  by_cols = c("run", "treatment"),
  target_names = c("roots", "soil", "biomass_vegetation")
)

dfs <- add_from_lookup(                     # copy pipe_nr col to dfs
  dfs = dfs,
  source_df = dfs$biomass_roots,
  new_col = "pipe_nr",
  by_cols = c("run", "treatment"),
  target_names = c("extracts", "roots", "soil", "vegetation")
)

dfs <- add_from_lookup(                     # copy pipe_diameter col to dfs
  dfs = dfs,
  source_df = dfs$biomass_roots,
  new_col = "pipe_diameter",
  by_cols = "pipe_nr",
  target_names = c("extracts", "roots", "soil", "vegetation", "biomass_vegetation")
)

dfs <- add_from_lookup(                     # copy single treatment cols to dfs
  dfs = dfs,
  source_df = dfs$biomass_vegetation,
  new_col = c("veg", "wet", "cut"),
  by_cols = "pipe_nr",
  target_names = c("extracts", "roots", "soil", "vegetation", "biomass_roots")
)

dfs <- add_from_lookup(                     # copy beriget col to non-vegetation dfs only
  dfs = dfs,
  source_df = dfs$soil,
  new_col = "beriget",
  by_cols = "pipe_nr",
  target_names = c("extracts", "roots", "biomass_roots")
)

dfs <- add_from_lookup(                     # copy beriget col to vegetation df only
  dfs = dfs,
  source_df = dfs$biomass_vegetation,
  new_col = "beriget",
  by_cols = "pipe_nr",
  target_names = "vegetation"
)

## ---- B. remove non-data rows ----

# too many headers
dfs$roots <- dfs$roots[-1, ]               # remove second header in roots 
dfs$vegetation <- dfs$vegetation[-1, ]     # and vegetation
dfs$extracts <- dfs$extracts[-1, ]         # and extracts

# remove rows: PEACH standards, too small samples, and expanded ranges from csv conversion
rules <- list(                                               # rules on what to remove
  roots      = list(col = "comment", patterns = c("PEACH")), # remove all PEACH samples
  soil       = list(col = "comment", patterns = c("PEACH")), # remove all PEACH samples
  vegetation = list(col = "comment", patterns = c("PEACH", "too small")),#PEACH and too small samples
  extracts   = list(col = "run", patterns = c("^\\s*$"))
)
dfs <- remove_rows_by_rules(dfs, rules)         # remove according to rules above

## ---- C. fix data types ----
# compare names across dfs
out <- compare_clean_names(dfs)            
dfs <- out$cleaned_data                    # isolate dfs new names from nested list

type_spec <- list(                         # list cols in each their data type
  logical   = c("beriget"),
  character = c("aboveground", "comment", "diameter", "layer",
    "pipe_nr", "run", "treatment", "tt", "cut", "veg", "wet"
  ),
  numeric   = c(
    "aboveground_weight_tot_g", "ape_pct_13c", "ape_pct_15n", "atom_pct_13c", "atom_pct_15n",
    "bottom_coarse_bag", "bottom_coarse_dry", "bottom_coarse_fresh",
    "bottom_fine_bag", "bottom_fine_dry", "bottom_fine_fresh", 
    "bryophyte_weight", "c_pct", "c_mg_pr_gdw", "c13_doc_ug_g", "c13_ug_pr_gc",
    "c13_ug_pr_gdw", "c13_raw", "c13korr", "cn_ratio", "doc_ug_g", "dtn_ug_g",
    "equ_weight", "graminoid_weight", "leaves_weight", "lichen_weight",
    "loi", "mg_c_in_sample", "mg_n_in_sample", 
    "mic_13c_ng_g", "mic_15n_ng_g", "mic_c_n", "mic_c_ug_g", "mic_n_ug_g",
    "mid_coarse_bag", "mid_coarse_dry", "mid_coarse_fresh", "mid_fine_bag",
    "mid_fine_dry", "mid_fine_fresh",
    "n_pct","n_mg_pr_gdw", "n15_dtn_ug_g", "n15_ug_pr_gdw", "n15_ug_pr_gn", "n15_raw",
    "n15korr", "natabun_13c_atm_pct", "natabun_15n_atm_pct",
    "pipe_diameter", "sample_weight", "stem_weight",
    "top_coarse_bag", "top_coarse_dry", "top_coarse_fresh",
    "top_fine_bag",  "top_fine_dry",  "top_fine_fresh",
    "weight_sorted_excl_bag_g"
  )
)

dfs <- lapply(dfs, coerce_types, type_spec = type_spec)      # apply type_spec to all dfs
str(c(dfs$vegetation, dfs$extracts, dfs$soil, 
      dfs$roots, dfs$biomass_roots, dfs$biomass_vegetation)) # print to check

# ---- D. handle missing values ----



# E. remove duplicates

# ---- F. filter technical outliers ----

# check the biplots for natural abundance values
p_veg <- ggplot(subset(dfs$vegetation, beriget == "FALSE"), # subset to only plot natabun
                aes(d13ckorr, d15nkorr, colour = aboveground)) +  geom_point() # show veg types
ggplotly(p_veg) # plot using plotly for an interactive plot

dfs$vegetation <- dfs$vegetation[dfs$vegetation$d13ckorr >= -33, ] # anything below -35 for deltaC are considered technical errors according to Anders Michelsen, therefore, removed.

# roots
p_roots <- ggplot(subset(dfs$roots, beriget =="FALSE"), # subset to only plot natabun
                  aes(d13ckorr, d15nkorr, colour = diameter)) +  geom_point() # show vegetation types
ggplotly(p_roots) # plot using plotly for an interactive plot

# soil
p_soil <- ggplot(subset(dfs$soil, beriget == "FALSE"), # subset to only plot natabun
                 aes(d13ckorr, d15nkorr, colour = layer)) +  geom_point() # show soil depth
ggplotly(p_soil) # plot using plotly for an interactive plot

# extracts
p_extracts <- ggplot(subset(dfs$extracts, beriget == "FALSE"), # subset to only plot natabun
                     aes(d13ckorr, d15nkorr, colour = layer)) +  geom_point() # show soil depth
ggplotly(p_extracts) # plot using plotly for an interactive plot

# dfs$extracts <- dfs$extracts[dfs$extracts$d13ckorr >= -31, ] # remove one obvious outlier?

### ---- 03 quality control checks ----

### ---- 04 export clean ----

paths <- c(
  extracts           = "data/clean/irms/extracts.csv",
  roots              = "data/clean/irms/roots.csv",
  soil               = "data/clean/irms/soil.csv",
  vegetation         = "data/clean/irms/vegetation.csv",
  biomass_roots      = "data/clean/irms/biomass_roots.csv",
  biomass_vegetation = "data/clean/irms/biomass_vegetation.csv"
)
Map(readr::write_csv, dfs[names(paths)], paths)
