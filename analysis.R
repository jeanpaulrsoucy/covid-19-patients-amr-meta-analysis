##### Antimicrobial resistance in patients with COVID-19: a systematic review and meta-analysis #####
##### Langford et al. #####
##### Script by Jean-Paul R. Soucy #####
##### https://github.com/jeanpaulrsoucy/covid-19-patients-amr-meta-analysis #####

# load libraries
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(broom)
library(boot)
library(meta)
library(metafor)

# load functions
source("funs.R")

# custom p-value formatting in meta::forest
source("formatPT.R")
assignInNamespace("formatPT", formatPT, "meta")

# create output directories (if they don't already exist)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# load data
dat <- read.csv("data.csv", stringsAsFactors = FALSE) %>%
  # delete blank rows
  filter(study != "")

# process variables
dat <- dat %>%
  mutate(
    patients_eval = as.integer(readr::parse_number(as.character(patients_eval))),
    Sample = as.integer(readr::parse_number(as.character(Sample))),
    End.Date = as.Date(End.Date, "%m/%d/%Y"),
    end_month_int = as.integer(case_when(
      year(End.Date) == 2020 ~ month(End.Date),
      year(End.Date) == 2021 ~ 12 + month(End.Date)
    ))
  )

# create factor variables
dat <- dat %>%
  mutate(
    setting = factor(
      case_when(
        Setting %in% c("Hospital", "hospital") ~ "Hospital",
        Setting == "Hospital ICU" ~ "Hospital ICU",
        Setting %in% c("Hospital/Outpatient", "Hospital/Outpatients") ~ "Hospital/Outpatient"
      ),
      levels = c("Hospital", "Hospital ICU", "Hospital/Outpatient")
    ),
    Age_Group = factor(
      case_when(
        as.numeric(Age) < 18 ~ "0-17",
        as.numeric(Age) >=18 & as.numeric(Age) < 66 ~ "18-65",
        as.numeric(Age) >= 66 ~ "66+",
        is.na(as.numeric(Age)) ~ "Not Specified"
      ),
      levels = c("0-17", "18-65", "66+", "Not Specified")
    ),
    body_site = factor(
      case_when(
        body_site == "resp" ~ "Resp",
        body_site == "blood" ~ "Blood",
        body_site == "multiple" ~ "Multiple",
        body_site == "unspecified" ~ "Not Specified"
      ),
      levels = c("Resp", "Blood", "Multiple", "Not Specified")
    ),
    quality_score = factor(
      case_when(
        Quality_score %in% 0:4 ~ "Low (0-4)",
        Quality_score %in% 5:7 ~ "Moderate (5-7)",
        Quality_score %in% 8:10 ~ "High (8-10)"
      ),
      levels = c("Low (0-4)", "Moderate (5-7)", "High (8-10)")
    ),
    study_end_month_grp = factor(
      case_when(
        End.Date < as.Date("2020-07-01") ~ "Jan-Jun 2020",
        End.Date < as.Date("2021-01-01") ~ "Jul-Dec 2020",
        End.Date < as.Date("2021-07-01") ~ "Jan-Jun 2021",
        is.na(End.Date) ~ "Not Specified"
    ),
    levels = c("Jan-Jun 2020", "Jul-Dec 2020", "Jan-Jun 2021", "Not Specified")
  ),
  who_region = factor(
    case_when(
      WHO_Region == "Africa" ~ "Africa",
      WHO_Region == "Eastern Mediterranean" ~ "Eastern Mediterranean",
      WHO_Region == "Europe" ~ "Europe",
      WHO_Region == "North America" ~ "Americas",
      WHO_Region == "South America" ~ "Americas",
      WHO_Region == "South-East Asia" ~ "South-East Asia",
      WHO_Region == "Western Pacific" ~ "Western Pacific",
      WHO_Region == "Multiple" ~ "Multiple"
    ),
    levels = c("Europe", "Africa", "Americas", "Eastern Mediterranean",
               "South-East Asia", "Western Pacific", "Multiple")
  ),
  income = factor(
    case_when(
      Income == "HIC" ~ "HIC",
      Income == "LMIC" ~ "LMIC",
      Income == "Both" ~ "Both"
    ),
    levels = c("HIC", "LMIC", "Both"),
  ),
  Infection_Type = factor(
    case_when(
      Infection_Type %in% c("Both", "Unspecified") ~ "Both/Unspecified",
      Infection_Type == "Co-infection" ~ "Co-infection",
      Infection_Type == "Secondary" ~ "Secondary"
    ),
    levels = c("Both/Unspecified", "Co-infection", "Secondary")
  ),
  Infection_Type_setting = factor(
    paste(Infection_Type, setting, sep = ", ")
  )
)

# create continuous variables
dat <- dat %>%
  mutate(
    age = as.numeric(Age),
    percent_female = as.integer(Female.patients..n.) / Sample * 100,
    percent_mechanical_vent = as.integer(dat$Mechanical.Ventilation) / Sample * 100,
    percent_ards = as.integer(ARDS) / Sample * 100,
    percent_icu = as.integer(ICU) / Sample * 100,
    percent_smoker = as.integer(Smokers) / Sample * 100,
    percent_copd = as.integer(COPD) / Sample * 100,
    percent_cvd = as.integer(CVD) / Sample * 100,
    percent_diabetes = as.integer(Diabetes) / Sample * 100,
    percent_malignancy = as.integer(Malignancy) / Sample * 100,
    percent_immunocompromised = as.integer(Immunocompromised) / Sample * 100,
    percent_corticosteroid = as.integer(corticosteroids) / Sample * 100,
    percent_il6 = as.integer(il_6_inhibitors) / Sample * 100,
    percent_abx = as.integer(abx) / Sample * 100
  ) %>%
  rowwise() %>%
  mutate(percent_severe = ifelse(
    all(is.na(percent_mechanical_vent), is.na(percent_ards), is.na(percent_icu)),
    NA,
    max(percent_mechanical_vent, percent_ards, percent_icu, na.rm = TRUE))
  ) %>%
  ungroup()

### 1: Prevalence of Bacterial Infection (patient-level) ###

# run plots for each of the 3 outcomes
for (outcome in c("coinfection", "secondary_infection", "bacterial_infection_unspecified")) {
  
  ## Plot 1: Stratified by setting (Setting)
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])), outcome, "setting", population = "patients_eval"), xmax = 70,
              out_png = paste0("figures/1_plot_1_", outcome, ".png"),
              out_width = 10, out_height = 5)
  # PDF version of plot
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])), outcome, "setting", population = "patients_eval"), xmax = 70,
              out_pdf = paste0("figures/1_plot_1_", outcome, ".pdf"),
              out_width = 10, out_height = 5)
  
  ## Plot 2: Stratified by age range (Age_Group)
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])), outcome, "Age_Group", population = "patients_eval"), xmax = 35,
              out_png = paste0("figures/1_plot_2_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 3: Stratified by body site (body_site)
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])), outcome, "body_site", population = "patients_eval"), xmax = 60,
              out_png = paste0("figures/1_plot_3_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 4: Stratified by quality score (Quality_score)
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])), outcome, "quality_score", population = "patients_eval"), xmax = 50,
              out_png = paste0("figures/1_plot_4_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 5: Sensitivity Analysis of Plot 1 - exclude studies with "Y" under "I_Prev_Sensi_exclude"
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])) %>%
                            filter(I_Prev_Sensi_exclude != "Y"), outcome, "setting", population = "patients_eval"), xmax = 70,
              out_png = paste0("figures/1_plot_5_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 6: Sensitivity Analysis of Plot 1 - exclude studies with "No" under "clinical_def"
  forest_plot(forest_calc(dat %>% filter(!is.na(dat[[outcome]])) %>%
                            filter(!clinical_def %in% c("No", "no")), outcome, "setting", population = "patients_eval"), xmax = 75,
              out_png = paste0("figures/1_plot_6_", outcome, ".png"),
              out_width = 10, out_height = 5)
}

### 2: Predictors of Bacterial Infection ###

## define outcomes
outcomes <- c("coinfection", "secondary_infection", "bacterial_infection_unspecified")

## define adjustment variables
adjust <- c("age", "percent_severe")

## run models
vars <- c("setting", "body_site", "study_end_month_grp", "who_region",
          "age", "percent_female", "percent_mechanical_vent", "percent_smoker",
          "percent_copd", "percent_cvd", "percent_diabetes", "percent_malignancy",
          "percent_immunocompromised", "percent_corticosteroid", "percent_il6")

## create summary table
summary_table_2 <- lapply(vars, function(v) {
  cat(v, fill = TRUE)
  run_rma_mod(dat, outcomes, v, adjust, population = "patients_eval") %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(v) %>%
    make_table_row(v, population = "patients_eval") %>%
    make_sample_columns(outcomes, v, adjust, population = "patients_eval")
}) %>%
  bind_rows() # combine all the rows

## replace some data (alternative adjustment sets)

## % Mechanical ventilation - don't control for percent_severe
summary_table_2 <- replace_rows(dat, outcomes,
                                var = "percent_mechanical_vent",
                                adjust = "age",
                                population = "patients_eval",
                                summary_table = summary_table_2,
                                characteristic = "% Mechanical ventilation",
                                n_rows = 1)

## Setting - don't control for percent_severe
summary_table_2 <- replace_rows(dat, outcomes,
                              var = "setting",
                              adjust = "age",
                              population = "patients_eval",
                              summary_table = summary_table_2,
                              characteristic = "Setting",
                              n_rows = 4)

## add table headers
summary_table_2 <- summary_table_2 %>%
  make_table_headers(outcomes, population = "patients_eval")

## write summary table
write.table(summary_table_2,
            "tables/summary_table_2.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)

### 3: Prevalence of Antibiotic Resistance ###

# subset rows with AMR_Prevalance_Evaluable == "Yes"
dat_amr <- dat %>% filter(AMR_Prev_Evaluable. %in% c("Yes", "yes"))

# run plots for patient-level and organism-level
lvls <- matrix(c(
  "Resistant_pts_total", "Pts_Denom", "Total patients",
  "Resistant_organisms", "Orgs_Denom", "Total organisms"
), byrow = TRUE, ncol = 3)
for (i in 1:nrow(lvls)) {
  
  ## set variables
  outcome = lvls[i, 1]
  population = lvls[i, 2]
  pop_lab = lvls[i, 3]
  
  ## Plot 1: Stratified by setting (Setting)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "setting", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_1_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 2: Stratified by age range (Age_Group)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "Age_Group", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_2_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 3: Stratified by body site (body_site)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "body_site", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_3_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 4: Stratified by WHO region (WHO_Region)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "who_region", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_4_", outcome, ".png"),
              out_width = 10, out_height = 8)
  
  ## Plot 5: Stratified by infection type (Infection_Type)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "Infection_Type", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_5_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 6: Stratified by quality score (quality_score)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "quality_score", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_6_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 7: Stratified by income (income)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "income", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_7_", outcome, ".png"),
              out_width = 10, out_height = 5)
  
  ## Plot 8: Stratified by infection type (Infection_Type) and setting (Setting)
  forest_plot(forest_calc(dat_amr %>% filter(!is.na(dat_amr[[outcome]])), outcome, "Infection_Type_setting", population), xmax = 100, pop_lab = pop_lab,
              out_png = paste0("figures/3_plot_8_", outcome, ".png"),
              out_width = 10, out_height = 7)
}

## funnel plots
mod_Pts_denom <- forest_calc(dat_amr %>% filter(!is.na(dat_amr[["Resistant_pts_total"]])), "Resistant_pts_total", "All", "Pts_Denom")
png("figures/3_funnel_plot_Resistant_pts_total.png", units = "in", width = 6, height = 6, res = 300)
funnel(mod_Pts_denom, random = TRUE, fixed = FALSE)
dev.off()
mod_Orgs_denom <- forest_calc(dat_amr %>% filter(!is.na(dat_amr[["Resistant_organisms"]])), "Resistant_organisms", "All", "Orgs_Denom")
png("figures/3_funnel_plot_Resistant_organisms.png", units = "in", width = 6, height = 6, res = 300)
funnel(mod_Orgs_denom, random = TRUE, fixed = FALSE)
dev.off()

### 4: Predictors of Antibiotic Resistance ###

## define adjustment variables
adjust <- c("age", "percent_severe")

## run models
vars <- c("setting", "body_site", "Infection_Type",
          "study_end_month_grp", "who_region", "income",
          "age", "percent_female", "percent_mechanical_vent", "percent_smoker",
          "percent_copd", "percent_cvd", "percent_diabetes", "percent_malignancy",
          "percent_immunocompromised", "percent_corticosteroid", "percent_il6",
          "percent_abx")

## create summary table

## patient-level
outcomes <- "Resistant_pts_total"
summary_table_4_patients <- lapply(vars, function(v) {
  cat(v, fill = TRUE)
  run_rma_mod(dat_amr, outcomes, v, adjust, population = "Pts_Denom") %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(v) %>%
    make_table_row(v, population = "Pts_Denom") %>%
    make_sample_columns(outcomes, v, adjust, population = "Pts_Denom")
}) %>%
  bind_rows() # combine all the rows

## replace some data (alternative adjustment sets)

## % Mechanical ventilation - don't control for percent_severe
summary_table_4_patients <- replace_rows(dat_amr, outcomes,
                                var = "percent_mechanical_vent",
                                adjust = "age",
                                population = "Pts_Denom",
                                summary_table = summary_table_4_patients,
                                characteristic = "% Mechanical ventilation",
                                n_rows = 1)

## Setting - don't control for percent_severe
summary_table_4_patients <- replace_rows(dat_amr, outcomes,
                                var = "setting",
                                adjust = "age",
                                population = "Pts_Denom",
                                summary_table = summary_table_4_patients,
                                characteristic = "Setting",
                                n_rows = 4)

## add table headers
summary_table_4_patients <- summary_table_4_patients %>%
  make_table_headers(outcomes, population = "Pts_Denom")

## organism-level
outcomes <- "Resistant_organisms"
summary_table_4_organisms <- lapply(vars, function(v) {
  cat(v, fill = TRUE)
  run_rma_mod(dat_amr, outcomes, v, adjust, population = "Orgs_Denom") %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(v) %>%
    make_table_row(v, population = "Orgs_Denom") %>%
    make_sample_columns(outcomes, v, adjust, population = "Orgs_Denom")
}) %>%
  bind_rows() # combine all the rows

## replace some data (alternative adjustment sets)

## % Mechanical ventilation - don't control for percent_severe
summary_table_4_organisms <- replace_rows(dat_amr, outcomes,
                                         var = "percent_mechanical_vent",
                                         adjust = "age",
                                         population = "Orgs_Denom",
                                         summary_table = summary_table_4_organisms,
                                         characteristic = "% Mechanical ventilation",
                                         n_rows = 1)

## Setting - don't control for percent_severe
summary_table_4_organisms <- replace_rows(dat_amr, outcomes,
                                         var = "setting",
                                         adjust = "age",
                                         population = "Orgs_Denom",
                                         summary_table = summary_table_4_organisms,
                                         characteristic = "Setting",
                                         n_rows = 4)

## add table headers
summary_table_4_organisms <- summary_table_4_organisms %>%
  make_table_headers(outcomes, population = "Orgs_Denom")

## add blank rows to patients table
summary_table_4_patients <- summary_table_4_patients %>%
  {add_row(
    .data = .,
    .after = which(.$terms == "Both/Unspecified"),
    data.frame(terms = "Co-infection", X1 = "(No data)", X2 = "(No data)", sample_size_col = "0"))} %>%
  {add_row(
    .data = .,
    .after = which(.$terms == "Europe"),
    data.frame(terms = "Africa", X1 = "(No data)", X2 = "(No data)", sample_size_col = "0"))} %>%
  {add_row(
    .data = .,
    .after = which(.$terms == "Western Pacific"),
    data.frame(terms = "Multiple", X1 = "(No data)", X2 = "(No data)", sample_size_col = "0"))} %>%
  {add_row(
    .data = .,
    .after = which(.$terms == "LMIC"),
    data.frame(terms = "Both", X1 = "(No data)", X2 = "(No data)", sample_size_col = "0"))}

## combine tables
# avoid double join on "Multiple" term
summary_table_4_patients[which(summary_table_4_patients$terms == "Multiple")[1], "terms"] <- "Multiple2"
summary_table_4_organisms[which(summary_table_4_organisms$terms == "Multiple")[1], "terms"] <- "Multiple2"
summary_table_4 <- full_join(summary_table_4_patients, summary_table_4_organisms, by = "terms")
summary_table_4$terms <- ifelse(summary_table_4$terms == "Multiple2", "Multiple", summary_table_4$terms)

## write summary table
write.table(summary_table_4,
            "tables/summary_table_4.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)

## create scatterplots for univariate continuous associations
vars_cont <- c("age", "percent_female", "percent_mechanical_vent", "percent_smoker",
               "percent_copd", "percent_cvd", "percent_diabetes", "percent_malignancy",
               "percent_immunocompromised", "percent_corticosteroid", "percent_il6",
               "percent_abx")

## patient-level
outcomes <- "Resistant_pts_total"
png(filename = "figures/3_scatterplots_unadjusted_Resistant_pts_total.png",
    width = 881, height = 1070)
par(mfrow = c(4, 3))
par(mar = c(5, 5, 1, 1))
scatterplots_patients <- lapply(vars_cont, function(v) {
  cat(v, fill = TRUE)
  if (v == "percent_mechanical_vent") {
    adjust <- "age" # alternative adjustment set
  } else {
    adjust <- c("age", "percent_severe")
  }
  mods <- run_rma_mod(dat_amr, outcomes, v, adjust = adjust, population = "Pts_Denom") %>%
    unlist(recursive = FALSE)
  mod <- mods[[1]] # unadjusted model
  xlab <- sub(" \\(10% increase)", "", get_var_name(v))
  xlab <- sub(" \\(10-year increase)", "", xlab)
  tryCatch({
    regplot(mod, v, xlab = xlab, psize = log(mod$ni + 1)/2, cex.axis = 1.5, cex.lab = 2)
    },
           error = function(e){warning("Skipped: ", v)}) # catch error if insufficient info to estimate model
})
dev.off()
par(mfrow = c(1, 1))
par(mar = c(1, 1, 1, 1))

## organism-level
outcomes <- "Resistant_organisms"
png(filename = "figures/3_scatterplots_unadjusted_Resistant_organisms.png",
    width = 881, height = 1070)
par(mfrow = c(4, 3))
par(mar = c(5, 5, 1, 1))
scatterplots_organisms <- lapply(vars_cont, function(v) {
  cat(v, fill = TRUE)
  if (v == "percent_mechanical_vent") {
    adjust <- "age" # alternative adjustment set
  } else {
    adjust <- c("age", "percent_severe")
  }
  mods <- run_rma_mod(dat_amr, outcomes, v, adjust, population = "Orgs_Denom") %>%
    unlist(recursive = FALSE)
  mod <- mods[[1]] # unadjusted model
  xlab <- sub(" \\(10% increase)", "", get_var_name(v))
  xlab <- sub(" \\(10-year increase)", "", xlab)
  tryCatch({
    regplot(mod, v, xlab = xlab, psize = log(mod$ni + 1)/2, cex.axis = 1.5, cex.lab = 2)
  },
  error = function(e){warning("Skipped: ", v)}) # catch error if insufficient info to estimate model
})
dev.off()
par(mfrow = c(1, 1))
par(mar = c(1, 1, 1, 1))

### 5: Bacterial infection Association with Mortality ###

# select potentially relevant studies
# Did this study evaluate bacterial co-infection or secondary infection outcomes compared to non-coinfection or non-secondary infection clinical outcomes? (i.e, mortality or LOS) -> Yes
# Adjusted? (at minimum age) -> Yes
dat_5 <- dat %>%
  filter(study %in% c(
    "Albelenda-Alonso, 2021", # HR
    "Bardi T, 2021", # OR
    # "Bhatt P, 2020", # "didn't provide OR or 95%CI, just p value" # may be able to extract raw numbers?
    "Cona A, 2021", # OR
    "Copaja-Corzo C, 2021", # HR
    "Goncalves Mendes Neto A, 2021", # OR
    "Moolla MS, 2021", # OR
    "Nasir N, 2021", # OR
    "Nebreda-Mayoral T, 2021", # OR
    "Nseir S, 2021", # HR
    "Petty L, 2021",  # OR
    "Ritter LA, 2021", # HR
    "Rouze A, 2021 (1)", # HR
    "Shafran N, 2021" # OR
  ))

# meta-analysis of odds ratios
dat_5_or <- dat_5 %>%
  filter(study %in% c(
    "Bardi T, 2021",
    "Cona A, 2021",
    "Goncalves Mendes Neto A, 2021",
    "Moolla MS, 2021",
    "Nasir N, 2021",
    "Nebreda-Mayoral T, 2021",
    "Petty L, 2021",
    "Shafran N, 2021"
  )) %>%
  select(study,
         setting,
         Mortality.OR.HR.Co.infection, LCL..Mortality.OR.HR.Co.infection., UCL..Mortality.OR.HR.Co.infection.,
         Mortality.OR.HR.Secondary.Infection, LCL..Mortality.OR.HR.Secondary.Infection., UCL..Mortality.OR.HR.Secondary.Infection.,
         Mortality.OR.HR.unspecified.both.types.Infection, LCL..Mortality.OR.HR.unspecified.both.types.Infection., UCL..Mortality.OR.HR.unspecified.both.types.Infection.)
names(dat_5_or) <- sub("LCL..Mortality.OR.HR.", "lcl_", names(dat_5_or))
names(dat_5_or) <- sub("UCL..Mortality.OR.HR.", "ucl_", names(dat_5_or))
names(dat_5_or) <- sub("Mortality.OR.HR.", "point_", names(dat_5_or))
dat_5_or <- dat_5_or %>%
  pivot_longer(
    cols = 3:11,
    names_to = "Infection_Type",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    value_type = case_when(
      grepl("^point_", Infection_Type) ~ "point",
      grepl("^lcl_", Infection_Type) ~ "lcl",
      grepl("^ucl_", Infection_Type) ~ "ucl"
    ),
    Infection_Type = factor(
      case_when(
        grepl("Secondary.Infection", Infection_Type) ~ "Secondary",
        grepl("Co.infection", Infection_Type) ~ "Co-infection",
        grepl("unspecified.both.types.Infection", Infection_Type) ~ "Both/Unspecified"
      ),
      levels = c("Secondary", "Co-infection", "Both/Unspecified")
    )
  ) %>%
  pivot_wider(
    id_cols = c(study, setting, Infection_Type),
    names_from = value_type,
    values_from = value
  )

# # no sub-groups
# png(file = "figures/5_plot_or.png", width = 10, height = 5, units = "in", res = 300)
# metagen(log(dat_5_or$point), lower = log(dat_5_or$lcl), upper = log(dat_5_or$ucl),
#         studlab = dat_5_or$study, sm = "HR", method.tau = "PM", fixed = FALSE) %>%
#   forest
# dev.off()

# sub-groups - by infection type
png(file = "figures/5_plot_or_by_infection_type.png", width = 10, height = 7, units = "in", res = 300)
metagen(log(dat_5_or$point), lower = log(dat_5_or$lcl), upper = log(dat_5_or$ucl),
        studlab = dat_5_or$study, sm = "OR", method.tau = "PM", fixed = FALSE,
        subgroup = dat_5_or$Infection_Type, subgroup.name = "Infection type") %>%
  forest(
    weight.study = "random",
    squaresize = 0.5,
    pooled.totals = TRUE,
    comb.fixed = FALSE,
    print.tau2 = TRUE,
    digits = 2,
    subgroup.name = "",
    col.by = "black"
  )
dev.off()

# # sub-groups - by setting
# png(file = "figures/5_plot_or_by_setting.png", width = 10, height = 7, units = "in", res = 300)
# metagen(log(dat_5_or$point), lower = log(dat_5_or$lcl), upper = log(dat_5_or$ucl),
#         studlab = dat_5_or$study, sm = "OR", method.tau = "PM", fixed = FALSE,
#         subgroup = dat_5_or$setting, subgroup.name = "Setting") %>%
#   forest(
#     weight.study = "random",
#     squaresize = 0.5,
#     pooled.totals = TRUE,
#     comb.fixed = FALSE,
#     print.tau2 = TRUE,
#     digits = 2,
#     subgroup.name = "",
#     col.by = "black"
#   )
# dev.off()

# sub-groups - by setting and infection type
dat_5_or <- arrange(dat_5_or, setting, Infection_Type) # sort so that subgroup order is sensible
png(file = "figures/5_plot_or_by_setting_infection_type.png", width = 10, height = 8, units = "in", res = 300)
metagen(log(dat_5_or$point), lower = log(dat_5_or$lcl), upper = log(dat_5_or$ucl),
        studlab = dat_5_or$study, sm = "OR", method.tau = "PM", fixed = FALSE,
        subgroup = paste(dat_5_or$setting, dat_5_or$Infection_Type, sep = " & "),
                         subgroup.name = "Setting and infection type") %>%
  forest(
    weight.study = "random",
    squaresize = 0.5,
    pooled.totals = TRUE,
    comb.fixed = FALSE,
    print.tau2 = TRUE,
    digits = 2,
    subgroup.name = "",
    col.by = "black"
  )
dev.off()

# meta-analysis of hazard ratios
dat_5_hr <- dat_5 %>%
  filter(study %in% c(
    "Albelenda-Alonso, 2021",
    "Copaja-Corzo C, 2021",
    "Nseir S, 2021",
    "Ritter LA, 2021",
    "Rouze A, 2021 (1)"
  )) %>%
  select(study,
         setting,
         Mortality.OR.HR.Co.infection, LCL..Mortality.OR.HR.Co.infection., UCL..Mortality.OR.HR.Co.infection.,
         Mortality.OR.HR.Secondary.Infection, LCL..Mortality.OR.HR.Secondary.Infection., UCL..Mortality.OR.HR.Secondary.Infection.,
         Mortality.OR.HR.unspecified.both.types.Infection, LCL..Mortality.OR.HR.unspecified.both.types.Infection., UCL..Mortality.OR.HR.unspecified.both.types.Infection.)
names(dat_5_hr) <- sub("LCL..Mortality.OR.HR.", "lcl_", names(dat_5_hr))
names(dat_5_hr) <- sub("UCL..Mortality.OR.HR.", "ucl_", names(dat_5_hr))
names(dat_5_hr) <- sub("Mortality.OR.HR.", "point_", names(dat_5_hr))
dat_5_hr <- dat_5_hr %>%
  pivot_longer(
    cols = 3:11,
    names_to = "Infection_Type",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    value_type = case_when(
      grepl("^point_", Infection_Type) ~ "point",
      grepl("^lcl_", Infection_Type) ~ "lcl",
      grepl("^ucl_", Infection_Type) ~ "ucl"
    ),
    Infection_Type = factor(
      case_when(
        grepl("Secondary.Infection", Infection_Type) ~ "Secondary",
        grepl("Co.infection", Infection_Type) ~ "Co-infection",
        grepl("unspecified.both.types.Infection", Infection_Type) ~ "Both/Unspecified"
      ),
      levels = c("Secondary", "Co-infection", "Both/Unspecified")
    )
  ) %>%
  pivot_wider(
    id_cols = c(study, setting, Infection_Type),
    names_from = value_type,
    values_from = value
  )

# # no sub-groups
# png(file = "figures/5_plot_hr.png", width = 10, height = 5, units = "in", res = 300)
# metagen(log(dat_5_hr$point), lower = log(dat_5_hr$lcl), upper = log(dat_5_hr$ucl),
#         studlab = dat_5_hr$study, sm = "HR", method.tau = "PM", fixed = FALSE) %>%
#   forest
# dev.off()

# sub-groups - by infection type
png(file = "figures/5_plot_hr_by_infection_type.png", width = 10, height = 5, units = "in", res = 300)
metagen(log(dat_5_hr$point), lower = log(dat_5_hr$lcl), upper = log(dat_5_hr$ucl),
        studlab = dat_5_hr$study, sm = "HR", method.tau = "PM", fixed = FALSE,
        subgroup = dat_5_hr$Infection_Type, subgroup.name = "Infection type") %>%
  forest(
    weight.study = "random",
    squaresize = 0.5,
    pooled.totals = TRUE,
    comb.fixed = FALSE,
    print.tau2 = TRUE,
    digits = 2,
    subgroup.name = "",
    col.by = "black"
  )
dev.off()

# # sub-groups - by setting
# png(file = "figures/5_plot_hr_by_setting.png", width = 10, height = 5, units = "in", res = 300)
# metagen(log(dat_5_hr$point), lower = log(dat_5_hr$lcl), upper = log(dat_5_hr$ucl),
#         studlab = dat_5_hr$study, sm = "HR", method.tau = "PM", fixed = FALSE,
#         subgroup = dat_5_hr$setting, subgroup.name = "Setting") %>%
#   forest(
#     weight.study = "random",
#     squaresize = 0.5,
#     pooled.totals = TRUE,
#     comb.fixed = FALSE,
#     print.tau2 = TRUE,
#     digits = 2,
#     subgroup.name = "",
#     col.by = "black"
#   )
# dev.off()

# sub-groups - by setting and infection type
dat_5_hr <- arrange(dat_5_hr, setting, Infection_Type) # sort so that subgroup order is sensible
png(file = "figures/5_plot_hr_by_setting_infection_type.png", width = 10, height = 5, units = "in", res = 300)
metagen(log(dat_5_hr$point), lower = log(dat_5_hr$lcl), upper = log(dat_5_hr$ucl),
        studlab = dat_5_hr$study, sm = "HR", method.tau = "PM", fixed = FALSE,
        subgroup = paste(dat_5_hr$setting, dat_5_hr$Infection_Type, sep = " & "),
        subgroup.name = "Setting and infection type") %>%
  forest(
    weight.study = "random",
    squaresize = 0.5,
    pooled.totals = TRUE,
    comb.fixed = FALSE,
    print.tau2 = TRUE,
    digits = 2,
    subgroup.name = "",
    col.by = "black"
  )
dev.off()

### 6: Prevalence of Antibiotic Resistance in Specific Organisms ###
summary_table_6 <- summary_individual_orgs(dat) %>%
  write.csv("tables/summary_table_6.csv", row.names = FALSE, na = "")
