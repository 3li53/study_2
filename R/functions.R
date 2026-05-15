#### functions for irms/biomass scripts

#### ---- 01_clean functions ----
### ---- 00 packages ----
# functions script for irms data
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

load_or_install(c("dplyr", "janitor", "purrr", "tibble"))

### ---- 02 initial cleaning ----
## ---- A. clean names ----

### compare names across dataframes
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

### remove columns from multiple dataframes
remove_cols_all <- function(dfs, cols_to_remove) {
  map(dfs, ~ .x %>%
        clean_names() %>%
        select(-any_of(cols_to_remove)))
}

### lookup and copy to other dfs

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

### rename multiple columns in dfs

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

## ---- B. remove non-data rows ----

### remove rows from specific dfs, cols and entries

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

## ---- C. fix data types ----

coerce_types <- function(df, type_spec) {                      # apply declared datatypes
  for (type in names(type_spec)) {                             # loop over target types
    cols <- intersect(names(df), type_spec[[type]])            # keep existing columns only
    if (type == "logical")   df[cols] <- lapply(df[cols], as.logical)   # coerce logical
    if (type == "numeric")   df[cols] <- lapply(df[cols], as.numeric)   # coerce numeric
    if (type == "character") df[cols] <- lapply(df[cols], as.character) # coerce character
  }
  df                                                           # return modified df
}



















# ----------------------- natural abundance baseline ---------------------------
calc_isotope_means <- function(df, group_var) {
  df %>%
    group_by({{ group_var }}) %>%
    summarise(
      avg_d15 = mean(d15n_korr, na.rm = TRUE),
      se_d15  = sd(d15n_korr, na.rm = TRUE) / sqrt(n()),
      ymin    = avg_d15 - se_d15,
      ymax    = avg_d15 + se_d15,
      
      avg_d13 = mean(d13c_korr, na.rm = TRUE),
      se_d13  = sd(d13c_korr, na.rm = TRUE) / sqrt(n()),
      xmin    = avg_d13 - se_d13,
      xmax    = avg_d13 + se_d13,
      
      .groups = "drop"
    )
}
#-------------------- atm pct helper -------------------------------------------

atom_pct <- function(delta, Rstd) {
  Rsample <- Rstd * (1 + (delta / 1000))
  100 * (Rsample / (1 + Rsample))
}

apply_baseline_correction <- function(
    irms_df,
    natabun_means_df,
    group_var,
    R15 = 0.003676,
    R13 = 0.011237
) {
  irms_df %>%
    left_join(natabun_means, by = rlang::as_name(rlang::enquo(group_var))) %>%
    mutate(
      nat_abun_15n_atm_pct = atom_pct(avg_d15,   R15),
      nat_abun_13c_atm_pct = atom_pct(avg_d13,   R13),
      
      atom_pct_15n = atom_pct(d15n_korr, R15),
      atom_pct_13c = atom_pct(d13c_korr, R13),
      
      ape_pct_15n = atom_pct_15n - nat_abun_15n_atm_pct,
      ape_pct_13c = atom_pct_13c - nat_abun_13c_atm_pct,
      
      n15_per_dw_ug_per_g = n_per_dw_mg_per_g * (ape_pct_15n / 100) * 1000,
      n15_per_n_ug_per_g  = mg_n              * (ape_pct_15n / 100) * 1000,
      
      c13_per_dw_ug_per_g = c_per_dw_mg_per_g * (ape_pct_13c / 100) * 1000,
      c13_per_c_ug_per_g  = mg_c              * (ape_pct_13c / 100) * 1000
    ) %>%
    select(-avg_d15, -ymin, -ymax, -avg_d13, -se_d13, -xmin, -xmax)
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