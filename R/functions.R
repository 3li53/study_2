#### functions for irms/biomass scripts

#### ---- 01_clean functions ----
###  ---- packages ----
# function to load or install and load packages
load_or_install <- function(pkgs) {
  for (pkg in pkgs) {
    cat("----\nChecking:", pkg, "\n")
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {  # Try to load first
      cat("✓ Already installed and loaded:", pkg, "\n")
    } else {
      cat("Installing:", pkg, "\n")
      try(install.packages(pkg), silent = TRUE)                 # Try installing
      if (require(pkg, character.only = TRUE, quietly = TRUE)) {# Try loading again
        cat("✓ Successfully installed and loaded:", pkg, "\n")
      } else {
        cat("✗ FAILED to load:", pkg, "\n")
      }
    }
  }
  cat("----\nDone.\n")
}
load_or_install(c("dplyr", "janitor", "purrr", "tibble")) #for this script

###  ---- initial cleaning ----
##   ---- A. clean names ----

# compare names across dataframes
compare_clean_names <- function(dfs) {              # collect inputs as named list
  dfs_clean <- map(dfs, ~ .x %>% clean_names())     # clean names
  name_list <- map(dfs_clean, names)                # extract names
  all_names <- sort(unique(unlist(name_list)))      # get all unique column names
  comparison <- map_dfc(name_list, function(nms) {  # Build comparison table
    all_names %in% nms
  }) %>%
    setNames(names(name_list)) %>%
    mutate(column = all_names, .before = 1)
  print(comparison, n = Inf)                        # Print nicely
  return(list(
    cleaned_data = dfs_clean,
    name_table = comparison
  ))                                                # Also return cleaned data for reuse
}

# remove columns from multiple dataframes
remove_cols_all <- function(dfs, cols_to_remove) {
  map(dfs, ~ .x %>%
        clean_names() %>%
        select(-any_of(cols_to_remove)))
}

# lookup and copy to other dfs
add_from_lookup <- function(dfs, source_df, new_col, by_cols, target_names) {
  source_df <- source_df %>%
    mutate(across(all_of(by_cols), as.character))            # ensure consistent types in source
  lookup <- source_df %>%
    select(all_of(c(by_cols, new_col))) %>%
    distinct()                                               # create lookup table
  for (name in target_names) {                               # join into selected dataframes
    dfs[[name]] <- dfs[[name]] %>%
      mutate(across(all_of(by_cols), as.character)) %>%      # ensure same types
      left_join(lookup, by = by_cols)
  }
  return(dfs)
}

# rename multiple columns in dfs
rename_many_if_present <- function(dfs, rename_list) {                 # function that renames cols across multiple dataframes in list
  dfs <- purrr::map(dfs, function(df) {                                # apply a function to each dataframe in the list
    for (old_name in names(rename_list)) {                             # loop over each old column name in the rename list
      new_name <- rename_list[[old_name]]                              # get the corresponding new column name
      if (old_name %in% names(df)) {                                   # check if the current dataframe actually contains this column
        df <- dplyr::rename(df, !!new_name := !!rlang::sym(old_name))  # if yes, rename the column using tidy evaluation
      }
    }
    return(df)                                                         # return the modified dataframe
  })
  return(dfs)                                                          # return the updated list of dataframes
}

##   ---- B. remove non-data rows ----

# remove rows from specific dfs, cols and entries
remove_rows_by_rules <- function(dfs, rules) {         # function that removes rows from selected dataframes based on matching rules
  for (nm in names(rules)) {                           # loop over dfs name from rules
    rule <- rules[[nm]]                                # extract rule for current df
    dfs[[nm]] <- dfs[[nm]] %>%                         # update current dataframe in list
      dplyr::filter(!if_any(                           # keep rows where the condition below is false
        all_of(rule$col),                              # select column to check (as specified in the rule)
        ~ grepl(                                       # test if the column contains any of the specified patterns
          paste(rule$patterns, collapse = "|"),        # combine multiple patterns into a single regular expression
          .x,                                          # apply the pattern matching to the column values
          ignore.case = TRUE                           # ignore letter case when matching text
        )
      ))
  }
  return(dfs)                                          # return the updated list of dataframes
}

##   ---- C. fix data types ----

# apply datatype to specified columns for each df
coerce_types <- function(df, type_spec) {                      # apply declared datatypes
  for (type in names(type_spec)) {                             # loop over target types
    cols <- intersect(names(df), type_spec[[type]])            # keep existing columns only
    if (type == "logical")   df[cols] <- lapply(df[cols], as.logical)   # coerce logical
    if (type == "numeric")   df[cols] <- lapply(df[cols], as.numeric)   # coerce numeric
    if (type == "character") df[cols] <- lapply(df[cols], as.character) # coerce character
  }
  df                                                           # return modified df
}

## ---- D. handle missing values
## ---- E. remove duplicates
## ---- F. filter

#### ---- 02_calculations functions ----

### calculate isotope means for baseline corrections

# summary function for isotope means, sd, ymin and ymax, n
calc_isotope_means <- function(df, group_var) {            # summarize isotope data by group
  df %>%
    group_by({{ group_var }}) %>%                          # group by supplied variable
    summarise(
      n       = n(),
      avg_d15 = mean(d15nkorr, na.rm = TRUE),              # mean δ15N
      se_d15  = sd(d15nkorr, na.rm = TRUE) / sqrt(n),      # SE δ15N
      ymin    = avg_d15 - se_d15,                          # lower δ15N bound
      ymax    = avg_d15 + se_d15,                          # upper δ15N bound
      avg_d13 = mean(d13ckorr, na.rm = TRUE),              # mean δ13C
      se_d13  = sd(d13ckorr, na.rm = TRUE) / sqrt(n),      # SE δ13C
      xmin    = avg_d13 - se_d13,                          # lower δ13C bound
      xmax    = avg_d13 + se_d13,                          # upper δ13C bound
      
      .groups = "drop"                                     # return ungrouped tibble
    )
}

# atom pct function
atom_pct <- function(delta, Rstd) {         # atom pct helper function
  Rsample <- Rstd * (1 + (delta / 1000))
  100 * (Rsample / (1 + Rsample))
}

# correct downstream calculations based on natabun group means
apply_baseline_correction <- function(
    df,               # data with sample isotope values
    natabun_means_df, # group-wise natural abundance means
    group_var,        # grouping variable for join
    R15 = 0.003676,   # air 15N/14N standard ratio
    R13 = 0.011237    # VPDB 13C/12C standard ratio
) {
  df %>%
    left_join(
      natabun_means, 
      by = rlang::as_name(rlang::enquo(group_var))  # join baseline by group
      ) %>%
    mutate(
      natabun_15n_atm_pct = atom_pct(avg_d15, R15), # baseline atom % 15N
      natabun_13c_atm_pct = atom_pct(avg_d13, R13), # baseline atom % 13C
      atom_pct_15n = atom_pct(d15nkorr, R15),       # sample atom % 15N
      atom_pct_13c = atom_pct(d13ckorr, R13),       # sample atom % 13C
      ape_pct_15n  = atom_pct_15n - natabun_15n_atm_pct, # atom % excess 15N
      ape_pct_13c  = atom_pct_13c - natabun_13c_atm_pct, # atom % excess 13C
      n15_ug_pr_gdw = n_mg_pr_gdw * (ape_pct_15n / 100) * 1000,  # ug 15N excess per g DW
      n15_ug_pr_gn = (ape_pct_15n / 100) * 1e6,                 # µg 15N excess per g N
      c13_ug_pr_gdw = c_mg_pr_gdw * (ape_pct_13c / 100) * 1000,  # ug 13C excess per g DW
      c13_ug_pr_gc = (ape_pct_13c / 100) * 1e6                  # µg 13C excess per g C
    ) %>%
    select(-avg_d15, -se_d15, -ymin, -ymax, -avg_d13, -se_d13, -xmin, -xmax, -n) # drop baseline summaries
}








# ---------------------- OUTLIER FUNCTION --------------------------------------

iqr_outlier <- function(df, group_var, numeric_vars) {
  df %>%
    group_by(across(all_of(group_var))) %>%
    mutate(across(
      all_of(numeric_vars),
      ~ (. < quantile(., 0.25, na.rm = TRUE) - 1.5 * IQR(., na.rm = TRUE)) |
        (. > quantile(., 0.75, na.rm = TRUE) + 1.5 * IQR(., na.rm = TRUE)),
      .names = "outlier_{.col}"
    )) %>%
    ungroup()
}






















# ----------------------- SCATTER PLOT FOR OUTLIERS ----------------------------
plot_scatter_outliers <- function(df, x, y, outlier_col, colour = colour_vars) {
  if (!outlier_col %in% names(df)) {
    stop(paste("Outlier column", outlier_col, "not found"))
  }
  ggplot(df, aes_string(x = x, y = y, colour = colour)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(
      data = df[df[[outlier_col]] == TRUE, ],
      colour = "black", shape = 4, size = 5, stroke = 1.2
    ) +
    geom_text(
      data = df[df[[outlier_col]] == TRUE, ],
      aes(label = nr),
      vjust = -0.6, size = 4, colour = "black"
    ) +
    theme_minimal(base_size = 16)
}

colour_vars <- c(
  "aboveground"
)