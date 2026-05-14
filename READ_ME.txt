project/
│
├── data/
│   ├── raw/
│   │   ├── irms/
│   │   ├── thaw_depth/
│   │   ├── loggers/
│   │   ├── cfux/
│   │   └── metadata/
│   │
│   ├── processed/
│   │   ├── irms/
│   │   ├── thaw_depth/
│   │   ├── loggers/
│   │   ├── cfux/
│   │   └── merged/
│
├── scripts/
│   ├── irms/
│   │   ├── 01_clean.R
│   │   ├── 02_calculations.R
│   │
│   ├── thaw_depth/
│   │   ├── 01_clean.R
│   │
│   ├── loggers/
│   │   ├── 01_clean.R
│   │   ├── 02_aggregation.R
│   │
│   ├── cfux/
│   │   ├── 01_clean.R
│   │   ├── 02_flux_calculations.R
│   │
│   ├── integration/
│   │   ├── 01_merge_all.R
│   │   ├── 02_derive_variables.R
│   │
│   ├── analysis/
│   │   ├── 01_statistics.R
│   │   ├── 02_multivariate.R
│   │
│   ├── figures/
│   │   ├── 01_main_figures.R
│   │
├── R/
│   ├── functions_isotopes.R
│   ├── functions_flux.R
│   ├── functions_time_series.R
│
├── output/
│   ├── figures/
│   ├── tables/