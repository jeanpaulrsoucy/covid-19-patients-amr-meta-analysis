##### Functions for: Antimicrobial resistance in patients with COVID-19: a systematic review and meta-analysis #####
##### Langford et al. #####
##### Script by Jean-Paul R. Soucy #####
##### https://github.com/jeanpaulrsoucy/covid-19-patients-amr-meta-analysis #####

## get global dataset by population
get_dat_by_pop <- function(population) {
  if (population == "patients_eval") {
    dat_sample <- get("dat", envir = .GlobalEnv)
  } else if (population == "Pts_Denom") {
    dat_sample <- get("dat_amr", envir = .GlobalEnv) %>%
      filter(!is.na(Pts_Denom))
  } else if (population == "Orgs_Denom") {
    dat_sample <- get("dat_amr", envir = .GlobalEnv) %>%
      filter(!is.na(Orgs_Denom))
  } else {
    stop("Error with population.")
  }
  dat_sample
}

## run rma.glmm model for specified outcome w/ or w/out moderators
rma_mod <- function(dat, outcome, term, population) {
  rma.glmm(
    xi = dat[[outcome]],
    ni = dat[[population]],
    measure = "PLO",
    method = "ML",
    mods = as.formula(paste0("~", paste(term, collapse = " + "))),
    data = dat)
}

## run rma_mod() after filtering dataset to complete cases (no NA)
run_rma_mod <- function(dat, outcomes, var, adjust, population, verbose = TRUE, adjusted = TRUE) {
  
  ## filter dataset
  if (adjusted) {
    ## filter dataset (use same dataset for adjusted and unadjusted analysis)
    lapply(outcomes, FUN = function(outcome) {
      ## filter to complete cases
      dat_sub <- dat[!is.na(dat[[outcome]]), ] # denominator: those with data for outcome
      n_before <- nrow(dat_sub)
      dat_sub <- dat_sub[, c(outcome, var, adjust, population)] %>%
        {.[complete.cases(.), ]}
      n_after <- nrow(dat_sub)
      
      ## report filtering results
      if (verbose) {
        cat(n_after, "/", n_before, " studies have complete data", sep = "", fill = TRUE)
      }
      
      ## run rma_mod for unadjusted and adjusted analysis
      vars_unadjusted <- c(var)
      vars_adjusted <- unique(c(var, adjust)) # unique in case of overlap between var and adjust
      var_list <- list(vars_unadjusted, vars_adjusted)
      lapply(var_list, function(model_vars) {
        tryCatch(
          rma_mod(dat = dat_sub, outcome = outcome, term = model_vars, population = population),
          error = function(e) {return(NA)}
        )
      })
    })
  } else {
    ## filter dataset (unadjusted only)
    lapply(outcomes, FUN = function(outcome) {
      ## filter to complete cases
      dat_sub <- dat[!is.na(dat[[outcome]]), ] # denominator: those with data for outcome
      n_before <- nrow(dat_sub)
      dat_sub <- dat_sub[, c(outcome, var, population)] %>%
        {.[complete.cases(.), ]}
      n_after <- nrow(dat_sub)
      
      ## report filtering results
      if (verbose) {
        cat(n_after, "/", n_before, " patients have complete data", sep = "", fill = TRUE)
      }
      
      ## run rma_mod for unadjusted analysis
      rma_mod(dat = dat_sub, outcome = outcome, term = var, population = population)
    })
  }
}

## extract prevalence odds ratio from rma.glmm model
extract_por <- function(var_vals, var) {
  # get var type
  var_type <- get_var_type(var)
  
  # get POR
  if (var_type == "continuous") {
    sprintf("%.2f", exp(var_vals * 10))
  } else if (var_type == "factor") {
    sprintf("%.2f", exp(var_vals))
  } else {
    "VALUE ERROR"
  }
}

## print prevalence with 95% CI
make_prev_ci <- function(mod) {
  out <- paste0(mod$estimate, " (", mod$conf.low, " to ", mod$conf.high, ")")
  return(out)
}

## get variable type
get_var_type <- function(var) {
  case_when(
    var %in% c("setting", "body_site", "study_end_month_grp",
               "Infection_Type", "who_region", "risk_of_bias",
               "income") ~ "factor",
    var %in% c("age", "percent_female", "percent_mechanical_vent",
               "percent_ards", "percent_icu", "percent_severe",
               "percent_smoker", "percent_copd", "percent_cvd",
               "percent_diabetes", "percent_malignancy", "percent_immunocompromised",
               "percent_corticosteroid", "percent_il6", "percent_abx",
               "end_month_int") ~ "continuous",
    TRUE ~ NA_character_
  )
}

## get outcome display name
get_outcome_name <- function(var) {
  case_when(
    var == "coinfection" ~ "Co-infection",
    var == "secondary_infection" ~ "Secondary infection",
    var == "bacterial_infection_unspecified" ~ "Unspecified bacterial infection",
    var == "Resistant_pts_total" ~ "Patient-level resistance",
    var == "Resistant_organisms" ~ "Organism-level resistance",
    TRUE ~ "MISSING OUTCOME NAME"
  )
}

## get variable display name
get_var_name <- function(var) {
  case_when(
    var == "setting" ~ "Setting",
    var == "body_site" ~ "Body site",
    var == "study_end_month_grp" ~ "Study end month",
    var == "Infection_Type" ~ "Infection type",
    var == "who_region" ~ "WHO region",
    var == "income" ~ "Income",
    var == "risk_of_bias" ~ "Risk of bias",
    var == "age" ~ "Age (10-year increase)",
    var == "percent_female" ~ "% Female (10% increase)",
    var == "percent_mechanical_vent" ~ "% Mechanical ventilation (10% increase)",
    var == "percent_ards" ~ "% ARDS (10% increase)",
    var == "percent_icu" ~ "% ICU (10% increase)",
    var == "percent_severe" ~ "% Severe (10% increase)",
    var == "percent_smoker" ~ "% Smoker (10% increase)",
    var == "percent_copd" ~ "% COPD (10% increase)",
    var == "percent_cvd" ~ "% CVD (10% increase)",
    var == "percent_diabetes" ~ "% Diabetes (10% increase)",
    var == "percent_malignancy" ~ "% Malignancy (10% increase)",
    var == "percent_immunocompromised" ~ "% Immunicompromised (10% increase)",
    var == "percent_corticosteroid" ~ "% Corticosteroid (10% increase)",
    var == "percent_il6" ~ "% IL-6 inhibitor (10% increase)",
    var == "percent_abx" ~ "% Antibiotics (10% increase)",
    var == "end_month_int" ~ "Study end month (10-month increase)",
    TRUE ~ "MISSING VARIABLE NAME"
  )
}

## get reference level of factor
get_factor_ref_level <- function(var) {
  # grab factor levels from global data frame (dat) by population
  dat <- get("dat", envir = .GlobalEnv)
  switch(
    var,
    "setting" = {levels(dat[["setting"]])[1]},
    "body_site" = {levels(dat[["body_site"]])[1]},
    "study_end_month_grp" = {levels(dat[["study_end_month_grp"]])[1]},
    "Infection_Type" = {levels(dat[["Infection_Type"]])[1]},
    "who_region" = {levels(dat[["who_region"]])[1]},
    "risk_of_bias" = {levels(dat[["risk_of_bias"]])[1]},
    "income" = {levels(dat[["income"]])[1]},
    "MISSING REFERENCE LEVEL"
  )
}

## extract data from a list of rma models
extract_rma_mod <- function(mod_list, var) {
  # get var type
  var_type <- get_var_type(var)
  
  # extract data
  lapply(mod_list, function(mod) {
    # NA if model cannot be estimated due to lack of data
    if (identical(mod, NA)) {
      # placeholder
      data.frame(
        term = "",
        estimate = "",
        conf.low = "",
        conf.high = ""
      )
    } else {
      tidy(mod, conf.int = TRUE, measure = "PLO", exponentiate = FALSE) %>%
        select(term, estimate, conf.low, conf.high) %>%
        filter(grepl(paste0("^", var), term)) %>%
        mutate(
          term = sub(paste0("^", var), "", term),
          estimate = extract_por(estimate, var),
          conf.low = extract_por(conf.low, var),
          conf.high = extract_por(conf.high, var)
        )
    }
  })
}

## make table row
make_table_row <- function(results_list, var, population) {
  
  # get variable type
  var_type <- get_var_type(var)
  match.arg(var_type,
            choices = c("factor", "continuous"),
            several.ok = FALSE)
  
  # make table row
  if (var_type == "factor") {
    # handle missing terms (e.g., if a factor level is missing for a particular outcome)
    # first, get all possible terms and order according to factor levels
    all_terms <- lapply(results_list, function(result) {
      result$term
    })
    all_terms <- unique(unlist(all_terms))
    dat_parent <- get_dat_by_pop(population)
    dat_parent[[var]] <- droplevels(dat_parent[[var]]) # drop empty levels
    all_terms <- all_terms[order(match(all_terms, levels(dat_parent[[var]])))]
    # then, add blank rows for missing terms
    results_list <- lapply(results_list, function(result) {
      missing_terms <- all_terms[!all_terms %in% result$term]
      if (length(missing_terms) > 0) {
        for (m in missing_terms) {
          result <- result %>%
            add_row(data.frame(term = m, estimate = "*", conf.low = "*", conf.high = "*"))
        }
      }
      # order rows
      result <- result[order(match(result$term, all_terms)), ]
      result
    })
    # continue building row
    term_col <- c(get_var_name(var), get_factor_ref_level(var), all_terms)
    value_cols <- lapply(results_list, function(result) {
      c("", "Reference", case_when(
        make_prev_ci(result) == "* (* to *)" ~ "(No data)",
        TRUE ~ make_prev_ci(result)))
    })
  } else {
    term_col <- get_var_name(var)
    value_cols <- lapply(results_list, function(result) {
      case_when(
        make_prev_ci(result) == " ( to )" ~ "(Insufficient data to estimate)",
        TRUE ~ make_prev_ci(result)
      )
    })
  }
  data.frame(
    terms = term_col,
    matrix(unlist(value_cols), ncol = length(value_cols), byrow = FALSE))
}

## make sample size columns
make_sample_columns <- function(dat, outcomes, var, adjust, population, adjusted = TRUE) {
  # grab sample sizes from global data frame (dat)
  dat_sample <- get_dat_by_pop(population)
  if (get_var_type(var) == "factor") {
    dat_sample[[var]] <- droplevels(dat_sample[[var]]) # drop empty levels
  }
  # add sample size columns
  if (adjusted) {
    for (i in 1:length(outcomes)) {
      dat_i <- dat_sample[, c(outcomes[i], var, adjust, population)] %>%
        {.[complete.cases(.), ]}
      sample_size <- nrow(dat_i)
      if (get_var_type(var) == "factor") {
        factor_levels_table <- table(dat_i[, var])
        # for study_end_month_grp: drop levels that do not appear in any result
        if (var == "study_end_month_grp")  {
          if (population == "patients_eval") {
            factor_levels_table <- factor_levels_table[names(factor_levels_table) != "Not Specified"]
          }
        }
        if (var == "who_region") {
          if (population == "Pts_Denom") {
            factor_levels_table <- factor_levels_table[names(factor_levels_table) %in% c(
              "Americas", "Europe", "Eastern Mediterranean", "South-East Asia", "Western Pacific")]
          } else if (population == "Orgs_Denom") {
            # factor_levels_table <- factor_levels_table[names(factor_levels_table) != "Not Specified"]
          }
        }
        if (var == "income") {
          if (population == "Pts_Denom") {
            factor_levels_table <- factor_levels_table[names(factor_levels_table) %in% c(
              "HIC", "LMIC")]
          } else if (population == "Orgs_Denom") {
            # placeholder
          }
        }
        if (var == "Infection_Type") {
          if (population == "Pts_Denom") {
            factor_levels_table <- factor_levels_table[!names(factor_levels_table) %in% c(
              "Co-infection")]
            } else if (population == "Orgs_Denom") {
              # placeholder
          }
        }
        factor_levels_n <- as.vector(factor_levels_table)
        sample_size_col <- c(sample_size, factor_levels_n)
      } else {
        sample_size_col <- sample_size
      }
      # print(dat)
      # print(sample_size_col)
      dat <- dat %>%
        add_column(sample_size_col, .after = 3 * i, .name_repair = "unique")
    }
  } else {
    for (i in 1:length(outcomes)) {
      dat_i <- dat_sample[, c(outcomes[i], var, population)] %>%
        {.[complete.cases(.), ]}
      sample_size <- nrow(dat_i)
      if (get_var_type(var) == "factor") {
        factor_levels_n <- as.vector(table(dat_i[, var]))
        sample_size_col <- c(sample_size, factor_levels_n)
      } else {
        sample_size_col <- sample_size
      }
      dat <- dat %>%
        add_column(sample_size_col, .after = 2 * i, .name_repair = "unique")
    }
  }
  
  # return data
  return(dat)
}

## make table header
make_table_headers <- function(tab, outcomes, population, adjusted = TRUE, exclude_not_specified = FALSE) {
  # grab sample sizes from global data frame (dat)
  dat_sample <- get_dat_by_pop(population)
  if (exclude_not_specified) {
    dat_sample <- dat_sample %>% filter(Type != "Not-specified")
  }
  # create header rows
  if (adjusted) {
    # create row 1
    row_1 <- c("Characteristic", rep("", each = 3 * length(outcomes)))
    # calculate sample size for each outcome
    for (i in 1:length(outcomes)) {
      outcome_n <- nrow(dat_sample[!is.na(dat_sample[[outcomes[i]]]), ])
      outcome_header <- paste0(get_outcome_name(outcomes[i]), " (n = ", outcome_n, ")")
      row_1[2 + 3 * (i - 1)] <- outcome_header
    }
    row_2 <- c("", c(rep(c("Unadjusted", "Adjusted", "Studies included"), times = length(outcomes))))
    header <- matrix(c(row_1, row_2), nrow = 2, byrow = TRUE,
                     dimnames = list(1:2, names(tab)))
  } else {
    # create row 1
    row_1 <- c("Characteristic", rep("", each = 2 * length(outcomes)))
    # calculate sample size for each outcome
    for (i in 1:length(outcomes)) {
      outcome_n <- nrow(dat_sample[!is.na(dat_sample[[outcomes[i]]]), ])
      outcome_header <- paste0(get_outcome_name(outcomes[i]), " (n = ", outcome_n, ")")
      row_1[2 + 2 * (i - 1)] <- outcome_header
    }
    row_2 <- c("", c(rep(c("Unadjusted", "Studies included"), times = length(outcomes))))
    header <- matrix(c(row_1, row_2), nrow = 2, byrow = TRUE,
                     dimnames = list(1:2, names(tab)))
  }
  
  ## return data with headers
  return(rbind(header, tab))
}

## replace rows in the summary table
replace_rows <- function(dat, outcomes, var, adjust, population, summary_table, characteristic, n_rows) {
  # generate new rows
  rows <- run_rma_mod(dat, outcomes, var, adjust, population) %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(var) %>%
    make_table_row(var, population) %>%
    make_sample_columns(outcomes, var, adjust, population)
  
  # find row number to begin replacement
  n_rows_begin <- grep(characteristic, summary_table$terms)
  
  # replace rows
  summary_table[n_rows_begin:(n_rows_begin + n_rows - 1), ] <- rows
  
  # return summary table
  return(summary_table)
}

## calculate forest plots for specified outcome w/ or w/out subgroups
forest_calc <- function(dat, outcome, type, population,
                        subgroup_n = TRUE # add the number of studies to the subgroup label
                        ) {
  
  if (type == "All") {
    # fit meta model
    metaprop(
      event = dat[[outcome]],
      n = dat[[population]],
      studlab = dat[["study"]],
      method = "GLMM",
      sm = "PLOGIT"
    )
  } else {
    # add number of studies to subgroup label
    if (subgroup_n) {
      table_subgroup <- as.integer(table(dat[[type]]))
      levels(dat[[type]]) <- paste0(levels(dat[[type]]), " (n = ", table_subgroup, ")")
    }
    # fit meta model
    metaprop(
      event = dat[[outcome]],
      n = dat[[population]],
      studlab = dat[["study"]],
      method = "GLMM",
      sm = "PLOGIT",
      subgroup = dat[[type]],
      subgroup.name = type
    )
  }
}

## plot forest plot
forest_plot <- function(dat,
                        xmin = 0, # x-axis minimum value
                        xmax = 100, # x-axis maximum value
                        pop_lab = "Total patients", # label for population (default: "Total patients")
                        order_subgroups = FALSE, # order subgroups from low to high prevalence
                        show_ind_studies = FALSE, # show individual studies?
                        out_png = NULL, # if specified, path to output png
                        out_pdf = NULL, # if specified, path to output pdf
                        out_width = 7, # width of output, in inches
                        out_height = 7 # height of output, in inches
){
  
  ## order subgroup results in decreasing order of prevalence
  if (order_subgroups == TRUE) {
    o <- order(dat$TE.random.w, decreasing = FALSE)
    for (var in c("bylevs", grep("\\.w$", names(dat)[!names(dat) %in% c("df.hakn.w", "df.Q.w")], value = TRUE))) {
      dat[[var]] <- dat[[var]][o]
    }
  }
  
  ## open graphics device to save plot
  if (!is.null(out_png)) {
    png(file = out_png, width = out_width, height = out_height, units = "in", res = 300)
  } else if (!is.null(out_pdf)) {
    pdf(file = out_pdf, width = out_width, height = out_height)
  }
  
  ## calculate number of subgroups (so that ALL subgroups are plotted, even those with 1 study)
  subgroup_true <- rep(TRUE, length(dat$TE.random.w))
  
  ## forest plot
  meta::forest(dat,
               xlim = c(xmin, xmax),
               pscale = 100,
               common = FALSE,
               rightcols = FALSE,
               leftcols = c("studlab", "n", "effect", "ci"),
               leftlabs = c("Subgroup",
                            pop_lab,
                            "Prevalence (%)",
                            "95% C.I."),
               xlab = "Prevalence (%)", smlab = "",
               weight.study = "random", squaresize = 0.5, col.square = "navy",
               col.square.lines = "navy",
               col.diamond = "maroon", col.diamond.lines = "maroon",
               pooled.totals = TRUE,
               comb.common = FALSE,
               fs.hetstat = 10,
               print.tau2 = TRUE,
               print.Q = TRUE,
               print.pval.Q = TRUE,
               print.I2 = TRUE,
               print.I2.ci = TRUE,
               digits.I2 = 1,
               digits = 1,
               digits.pval = 4,
               digits.pval.Q = 4,
               subgroup.name = "",
               col.by = "black",
               study.results = show_ind_studies,
               subgroup = subgroup_true
  )
  
  ## close graphics device
  if (!is.null(out_png) | !is.null(out_pdf)) {
    dev.off()
  }
  
}

# summary table of AMR prevalence in individual organisms
summary_individual_orgs <- function(dat) {
  # list of organisms with denominator and numerator(s)
  orgs <- list(
    list(org = "S. aureus", denom = "SA",
         num = c("SA_R"),
         num_name = c("MRSA")
         ),
    list(org = "Enterococcus spp.", denom = "Enterococcus",
         num = c("Enterococcus_R"),
         num_name = c("VRE")
         ),
    list(org = "S. pneumoniae", denom = "S_pneumo",
         num = c("S_pneumo_R", "S_pneumo_pen_R", "S_pneumo_FQ_R"),
         num_name = c("Any resistance", "Penicillin-Resistant S. pneumoniae", "Fluoroquinolone-Resistant S. pneumoniae")),
    list(org = "Pseudomonas spp.", denom = "Psa",
         num = c("Psa_R", "Psa_MDR", "Psa_CR", "Psa_ESBL", "Psa_ColR"),
         num_name = c("Any resistance", "Multi-drug Resistance", "Carbapenem-Resistance", "ESBL-producing/3GC resistant", "Colistin Resistance")
         ),
    list(org = "Klebsiella spp.", denom = "Kleb",
         num = c("Kleb_R", "Kleb_MDR", "Kleb_CR", "Kleb_ESBL", "Kleb_ColR"),
         num_name = c("Any resistance", "Multi-drug Resistance", "Carbapenem-Resistance", "ESBL-producing/3GC resistant", "Colistin Resistance")
         ),
    list(org = "E. coli", denom = "Ecoli",
         num = c("EColi_R", "Ecoli_MDR", "Ecoli_CR", "Ecoli_ESBL", "Ecoli_ColR"),
         num_name = c("Any resistance", "Multi-drug Resistance", "Carbapenem-Resistance", "ESBL-producing/3GC resistant", "Colistin Resistance")
         ),
    list(org = "Acinetobacter spp.", denom = "Acineto",
         num = c("Acineto_R", "Acineto_MDR", "Acineto_CR", "Acineto_ColR"),
         num_name = c("Any resistance", "Multi-drug Resistance", "Carbapenem-Resistance", "Colistin Resistance")
         ),
    list(org = "Enterobacter spp.", denom = "Ebacter",
         num = c("Ebacter_R_isolates", "Ebacter_MDR", "Ebacter_CR", "Ebacter_ESBL"),
         num_name = c("Any resistance", "Multi-drug Resistance", "Carbapenem-Resistance", "ESBL-producing/3GC resistant")
         ),
    list(org = "Stenotrophomonas spp.", denom = "Steno",
         num = c("Steno_MDR"),
         num_name = c("Multi-drug Resistance")
         ),
    list(org = "Serratia spp.", denom = "Serratia",
         num = c("Serratia_R", "Serratia_ESBL"),
         num_name = c("Any resistance", "ESBL-producing/3GC resistant")
         ),
    list(org = "Proteus spp.", denom = "Proteus",
         num = c("Proteus_R", "Proteus_MDR"),
         num_name = c("Any resistance", "Multi-drug Resistance")
         )
  )
  # create table for outputs
  out <- data.frame(
    org = rep(sapply(orgs, function(x) x[["org"]]), times = sapply(orgs, function(x) length(x[["num"]]))),
    denom_val = rep(sapply(orgs, function(x) x[["denom"]]), times = sapply(orgs, function(x) length(x[["num"]]))),
    num_val = unlist(sapply(orgs, function(x) x[["num"]])),
    num_name = unlist(sapply(orgs, function(x) x[["num_name"]])),
    n_studies = NA,
    num_denom = NA,
    prev = NA,
    hetero = NA
  )
  # count studies
  for (i in 1:length(orgs)) {
    out[out$org == orgs[[i]][["org"]] & (
      out$num_name %in% c("Any resistance", "MRSA", "VRE") | out$org == "Stenotrophomonas spp."),
        "n_studies"] <- dat %>%
      filter(!is.na(!!sym(orgs[[i]][["num"]][1]))) %>%
      nrow()
  }
  # extract outputs
  for (i in 1:nrow(out)) {
    # subset data
    d <- dat %>%
      transmute(
        study,
        num = !!sym(out[i, "num_val", drop = TRUE]),
        denom = !!sym(out[i, "denom_val", drop = TRUE])
      ) %>%
      filter(!is.na(num))
    # extract num and denom
    out[i, "num_denom"] <- paste(sum(d$num), sum(d$denom), sep = " / ")
    # fit model
    mod <- forest_calc(d, outcome = "num", type = "All", population = "denom")
    # extract prevalence estimate
    out[i, "prev"] <- paste0(sprintf("%.2f", inv.logit(mod$TE.random) * 100),
                             " (",
                             sprintf("%.2f", inv.logit(mod$lower.random) * 100),
                             "–",
                             sprintf("%.2f", inv.logit(mod$upper.random) * 100),
                             ")")
    # extract heterogeneity estimate
    out[i, "hetero"] <- paste0(sprintf("%.1f", mod$I2 * 100),
                               "% (",
                               sprintf("%.1f", mod$lower.I2 * 100),
                               "%–",
                               sprintf("%.1f", mod$upper.I2 * 100),
                               "%)")
  }
  # replace missing heterogeneity values/CIs
  out$hetero <- case_when(
    out$hetero == "0.0% (NA%–NA%)" ~ "0.0%",
    out$hetero == "NA% (NA%–NA%)" ~ "Not applicable",
    TRUE ~ out$hetero
  )
  
  # return table
  out
}
