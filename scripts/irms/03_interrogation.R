# 03 data interrogation

# density plots

#d13C

density_facets(
  dfs$extracts,
  c(mic_13c_ng_g, mic_15n_ng_g, mic_c_n),
  outfile = "output/figures/irms_interrogation/density_facets/mic_extracts_density.pdf"
)

density_facets(
  dfs$extracts,
  c(n15_dtn_ug_g, c13_doc_ug_g, ),
  outfile = "output/figures/irms_interrogation/density_facets/extracts_density.pdf"
)

density_facets(
  dfs$extracts,
  c(n15_dtn_ug_g, c13_doc_ug_g, ),
  outfile = "output/figures/irms_interrogation/density_facets/extracts_density.pdf"
)

# violin plots

# inorganic carbon //

# lipid effect //

# biplots

# isotopes by depth

