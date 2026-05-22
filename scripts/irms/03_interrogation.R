# 03 data interrogation
load_or_install(c("dplyr", "purrr", "ggplot2"))

# density plots

# extracts
density_facets(dfs$extracts, c(d13ckorr, d15nkorr), # density of dissolved N and organic C
               outfile = "output/figures/irms_interrogation/density_facets/extracts_density.pdf") 
density_facets(dfs$extracts, c(doc_ug_g, mic_c_ug_g), # density of dissolved N and organic C
               outfile = "output/figures/irms_interrogation/density_facets/extracts_density.pdf") 
density_facets(dfs$extracts, c(dtn_ug_g, mic_n_ug_g), # density of microbial n and c and cn ratio
  outfile = "output/figures/irms_interrogation/density_facets/mic_extracts_density.pdf")
density_facets(dfs$extracts, c(c13_doc_ug_g, mic_13c_ng_g), # density of dissolved c isotopes in soil and in microbes
  outfile = "output/figures/irms_interrogation/density_facets/extracts_density.pdf")
density_facets(dfs$extracts, c(n15_dtn_ug_g, mic_15n_ng_g), # density of dissolved isotopes
               outfile = "output/figures/irms_interrogation/density_facets/extracts_density.pdf")

# roots
density_facets(dfs$roots, c(d13ckorr, d15nkorr, cn_ratio),
               outfile = "output/figures/irms_interrogation/density_facets/roots_density.pdf")
density_facets(dfs$roots, c(c13_ug_pr_gc, n15_ug_pr_gn, cn_ratio),
               outfile = "output/figures/irms_interrogation/density_facets/roots_density.pdf")

# soil
density_facets(dfs$soil, c(d13ckorr, d15nkorr, cn_ratio),
               outfile = "output/figures/irms_interrogation/density_facets/soil_density.pdf")
density_facets(dfs$soil, c(ape, n15_ug_pr_gn, cn_ratio),
               outfile = "output/figures/irms_interrogation/density_facets/soil_density.pdf")

# vegetation
density_facets(dfs$vegetation, c(d13ckorr, d15nkorr, cn_ratio),
               outfile = "output/figures/irms_interrogation/density_facets/vegetation_density.pdf")
density_facets(dfs$vegetation, c(c13_ug_pr_gc, n15_ug_pr_gn, cn_ratio),
               outfile = "output/figures/irms_interrogation/density_facets/vegetation_density.pdf")


# biplot 

all_natabun_isotopes <- imap_dfr(dfs, ~
                           if(all(c("d13ckorr", "d15nkorr", "beriget") %in% names(.x))) {
                             .x %>% 
                               filter(beriget == FALSE) %>%       # keep only non-enriched
                               select(d13ckorr, d15nkorr) %>%     # select columns
                               mutate(dataset = .y)
                           } else {
                             NULL
                           }
)


ggplot(all_isotopes, aes(d13ckorr, d15nkorr, colour = dataset)) +
  geom_point()



# violin plots

# extracts

violin_jitter_plot(
  df     = dfs$extracts,
  x      = "veg",  y      = "doc_ug_g",
  colour = "wet",  facet  = "cut",
  x_lab  = NULL,   y_lab  = expression(delta^{13}*C),
  title  = expression("δ"^{13}*"C in DOC")
)
violin_jitter_plot(
  df     = dfs$extracts,
  x      = "veg",  y      = "mic_13c_ng_g",
  colour = "wet",  facet  = "cut",
  x_lab  = NULL,   y_lab  = expression(delta^{13}*C),
  title  = expression("Microbial δ"^{13}*"C")
)

violin_facet_plot(
  data        = dfs$extracts,
  pivot_cols  = c(mic_13c_ng_g, mic_15n_ng_g, mic_c_n),
  x           = "veg",
  colour      = "wet",
  facet_rows  = "Variable",   # <- string, because created inside function
  facet_cols  = "cut",
  title       = ""
)

# inorganic carbon //

# lipid effect //

# biplots


# isotopes by depth




# structure

# ---- ANALYSIS FRAMEWORK: GLYCINE TRACER EXPERIMENT ----
# Goal: quantify tracer uptake, processing, and allocation across soil–plant–microbe system

# =========================
# 1. DATA PREPARATION
# =========================
# - Convert δ13C and δ15N to atom fraction
# - Calculate atom % excess using controls
# - Use excess values for all downstream analyses (not raw δ)

# =========================
# 2. LAYER 1: PATTERN EXPLORATION (WHERE IS THE TRACER?)
# =========================
# Use violin / boxplots to assess distribution of enrichment

# Group data by:
# - pool: microbial, dissolved, roots, soil, vegetation, CO2
# - time: thaw stage (critical dimension)
# - depth / functional group where relevant

# Key questions:
# - which pools show strongest enrichment?
# - how does uptake change over time?
# - is there depth or functional structure?

# Note:
# - plot 13C and 15N separately
# - treat CO2 as a time series (line plot recommended)

# =========================
# 3. LAYER 2: BIPLOTS (C–N COUPLING; HOW IS TRACER PROCESSED?)
# =========================
# Plot excess13C vs excess15N

# Do NOT combine all pools → plot separately:
# - microbial + dissolved pools
# - roots / vegetation
# - soil
# - CO2 (optional)

# Interpretation:
# - strong linear relationship → coupled uptake (intact glycine)
# - slope deviation / scatter → partial processing
# - decoupling → C lost (respiration) or N retained

# =========================
# 4. LAYER 3: QUANTIFICATION
# =========================

# 4A. Coupling models
# - fit linear models (or LMMs) of excess13C ~ excess15N
# - compare slopes among pools
# → reveals differences in processing pathways

# 4B. Two-source mixing (glycine vs background)
# - estimate fraction of pool derived from glycine
# - apply per element (C and N separately)
# → gives relative contribution of tracer

# 4C. Mass balance (core analysis)
# - multiply excess by pool size → tracer amounts
# - sum across pools
# → quantify partitioning of tracer (absolute allocation)

# =========================
# 5. SYNTHESIS (WHAT DOES IT MEAN?)
# =========================
# Combine results to answer:

# 1. Where did the tracer go? → mass balance
# 2. How was it processed? → C–N coupling
# 3. Who uses it and how consistently? → distribution plots

# Conceptual interpretation:
# - microbes: early uptake + transformation
# - CO2: carbon loss (mineralisation)
# - plants/roots: uptake pathway (direct vs indirect)
# - soil: mixing / stabilisation pool

# =========================
# KEY PRINCIPLES
# =========================
# - Always analyse C and N together AND separately
# - Never interpret δ-space as “niche” (this is a tracer system)
# - Separate pools, processes, and time in all plots
# - Mass balance = quantitative backbone of analysis
