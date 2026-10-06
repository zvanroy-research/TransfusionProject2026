# ==============================================================================
# script 0: cohort generation
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(lubridate)
})

summarize <- dplyr::summarize
filter    <- dplyr::filter

root_dir   <- "/Users/zvr/Downloads/Clinical2.0"
file_path  <- "/Users/zvr/Downloads/Transfusion_Data.xlsx"

output_dir <- file.path(root_dir, "processed_data")
table_dir  <- file.path(root_dir, "tables_output")
fig_dir    <- file.path(root_dir, "figures_output")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(table_dir))  dir.create(table_dir, recursive = TRUE)
if (!dir.exists(fig_dir))    dir.create(fig_dir, recursive = TRUE)

baseline_intake_cutoff_days <- 7

parse_date_robust <- function(x) {
  d <- suppressWarnings(as.Date(parse_date_time(
    as.character(x), 
    orders = c("ymd HMS", "ymd HM", "ymd", "mdy HMS", "mdy HM", "mdy", "dmy HMS", "dmy HM", "dmy")
  )))
  na_idx <- which(is.na(d) & !is.na(x) & nzchar(trimws(as.character(x))))
  if (length(na_idx) > 0) {
    num_vals <- suppressWarnings(as.numeric(as.character(x[na_idx])))
    valid_num <- !is.na(num_vals) & num_vals > 20000 & num_vals < 80000
    if (any(valid_num)) {
      d[na_idx[valid_num]] <- as.Date(num_vals[valid_num], origin = "1899-12-30")
    }
  }
  return(d)
}

df_demo_raw   <- read_excel(file_path, sheet = "Demographics")
df_trans_raw  <- read_excel(file_path, sheet = "Transfusion")
df_events_raw <- read_excel(file_path, sheet = "Patient Events")
df_screen_raw <- read_excel(file_path, sheet = "Screening")
df_mhc1_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-I")
df_mhc2_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-II")

descriptions_to_exclude <- c(
  "Kidney Living Donor", "HSC Donor", "Disease Association",
  "Altruistic Donor", "OTHER/RECIPIENT", "Liver Living Donor", "OTHER/DONOR", 
  "HSC Recipient", "HSC/Kidney Recipient"
)

eth_col_name <- names(df_demo_raw)[str_detect(names(df_demo_raw), regex("ethnic|race", ignore_case = TRUE))][1]

df_patients_clean <- df_demo_raw %>%
  mutate(clean_desc = str_trim(description)) %>%
  filter(!clean_desc %in% str_trim(descriptions_to_exclude)) %>%
  filter(sex %in% c("M", "F")) %>% 
  mutate(
    pat_id         = as.character(patient_number),
    listing_date   = as.Date(as.numeric(create_date), origin = "1899-12-30"),
    birth_date     = as.Date(as.numeric(dob), origin = "1899-12-30"),
    age_at_listing = as.numeric(difftime(listing_date, birth_date, units = "weeks")) / 52.25,
    sex            = factor(sex, levels = c("M", "F")),
    raw_desc       = clean_desc,
    # Directly capture ethnic_code as character:
    ethnicity      = as.character(ethnic_code),
    
    organ_group    = case_when(
      str_detect(clean_desc, regex("Isolated Kidney", ignore_case = TRUE)) ~ "Isolated Kidney",
      str_detect(clean_desc, regex("Isolated Heart", ignore_case = TRUE))  ~ "Isolated Heart",
      str_detect(clean_desc, regex("Isolated Lung", ignore_case = TRUE))   ~ "Isolated Lung",
      # Robust fallbacks for single organs if "Isolated" prefix isn't strictly used in some raw rows
      clean_desc %in% c("Kidney", "Kidney Recipient", "Renal")   ~ "Isolated Kidney",
      clean_desc %in% c("Heart", "Heart Recipient", "Cardiac")   ~ "Isolated Heart",
      clean_desc %in% c("Lung", "Lung Recipient")                ~ "Isolated Lung",
      # Everything else (Liver, Pancreas, Bowel, Multi-Visceral) drops here
      TRUE ~ "Multi-Organ & Other"
    ),
    organ_group    = factor(organ_group, levels = c("Isolated Kidney", "Isolated Heart", "Isolated Lung", "Multi-Organ & Other"))
    # ---------------------------------------------------------
    
  ) %>%
  select(pat_id, listing_date, age_at_listing, sex, organ_group, raw_desc, ethnicity)

message("[00_build_cohort] Processing Screening and Reflex SAB cutoffs...")
df_screen_clean <- df_screen_raw %>%
  rename(patient_number = 1) %>%
  mutate(
    pat_id          = as.character(patient_number),
    draw_date_clean = parse_date_robust(draw_date),
    screen_pos      = if_else(str_to_upper(str_trim(Class_I_Result)) == "POS" | 
                                str_to_upper(str_trim(Class_II_Result)) == "POS", 1, 0)
  ) %>%
  filter(!is.na(draw_date_clean)) %>%
  select(pat_id, draw_date_clean, screen_pos)

df_sab_positive_draws <- bind_rows(
  df_mhc1_raw %>% 
    rename(patient_number = 1) %>%
    mutate(pat_id = as.character(patient_number), draw_date_clean = parse_date_robust(draw_date)) %>%
    select(-patient_number, -draw_date) %>%
    pivot_longer(cols = -c(pat_id, draw_date_clean), names_to = "Target", values_to = "MFI") %>%
    mutate(MFI = suppressWarnings(as.numeric(MFI))) %>%
    filter(!is.na(MFI) & !is.na(draw_date_clean)) %>%
    filter(
      (str_detect(Target, regex("^C\\*", ignore_case = TRUE)) & MFI >= 5000) |
        (!str_detect(Target, regex("^C\\*", ignore_case = TRUE)) & MFI >= 3000)
    ) %>%
    select(pat_id, draw_date_clean),
  
  df_mhc2_raw %>% 
    rename(patient_number = 1) %>%
    mutate(pat_id = as.character(patient_number), draw_date_clean = parse_date_robust(draw_date)) %>%
    select(-patient_number, -draw_date) %>%
    pivot_longer(cols = -c(pat_id, draw_date_clean), names_to = "Target", values_to = "MFI") %>%
    mutate(MFI = suppressWarnings(as.numeric(MFI))) %>%
    filter(!is.na(MFI) & !is.na(draw_date_clean)) %>%
    filter(
      (str_detect(Target, regex("^DP", ignore_case = TRUE)) & MFI >= 5000) |
        (!str_detect(Target, regex("^DP", ignore_case = TRUE)) & MFI >= 1500)
    ) %>%
    select(pat_id, draw_date_clean)
) %>%
  distinct(pat_id, draw_date_clean)

df_sab_all_draws <- bind_rows(
  df_mhc1_raw %>% rename(patient_number = 1) %>% mutate(pat_id = as.character(patient_number), draw_date_clean = parse_date_robust(draw_date)) %>% select(pat_id, draw_date_clean),
  df_mhc2_raw %>% rename(patient_number = 1) %>% mutate(pat_id = as.character(patient_number), draw_date_clean = parse_date_robust(draw_date)) %>% select(pat_id, draw_date_clean)
) %>%
  filter(!is.na(draw_date_clean)) %>%
  distinct(pat_id, draw_date_clean)

df_all_confirmed_positives <- df_screen_clean %>%
  filter(screen_pos == 1) %>%
  inner_join(df_sab_positive_draws, by = "pat_id", suffix = c("_screen", "_sab"), relationship = "many-to-many") %>%
  mutate(day_diff = as.numeric(difftime(draw_date_clean_sab, draw_date_clean_screen, units = "days"))) %>%
  filter(day_diff >= 0 & day_diff <= 30) %>%
  mutate(confirmed_pos_date = draw_date_clean_screen) %>%
  select(pat_id, confirmed_pos_date) %>%
  distinct()

df_all_evals <- bind_rows(
  df_screen_clean %>% select(pat_id, draw_date_clean),
  df_sab_all_draws
) %>%
  distinct(pat_id, draw_date_clean)

df_last_lab_eval <- df_all_evals %>%
  group_by(pat_id) %>%
  summarize(last_eval_date = max(draw_date_clean, na.rm = TRUE), .groups = "drop")

message("[00_build_cohort] Qualifying baseline cohort...")
df_baseline_evals <- df_all_evals %>%
  inner_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(draw_date_clean <= (listing_date + baseline_intake_cutoff_days)) %>%
  group_by(pat_id) %>%
  summarize(baseline_anchor_date = max(draw_date_clean, na.rm = TRUE), .groups = "drop")

df_baseline_sensitized <- df_all_confirmed_positives %>%
  inner_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(confirmed_pos_date <= (listing_date + baseline_intake_cutoff_days)) %>%
  distinct(pat_id)

df_baseline_pos_screens <- df_screen_clean %>%
  inner_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(draw_date_clean <= (listing_date + baseline_intake_cutoff_days) & screen_pos == 1) %>%
  select(pat_id, screen_date = draw_date_clean)

df_pos_screens_with_paired_sab <- df_baseline_pos_screens %>%
  inner_join(df_sab_all_draws, by = "pat_id") %>%
  mutate(day_diff = as.numeric(difftime(draw_date_clean, screen_date, units = "days"))) %>%
  filter(day_diff >= 0 & day_diff <= 30) %>%
  distinct(pat_id)

df_unconfirmed_pos_screens <- df_baseline_pos_screens %>%
  anti_join(df_pos_screens_with_paired_sab, by = "pat_id") %>%
  distinct(pat_id)

df_total_baseline_exclusions <- bind_rows(df_baseline_sensitized, df_unconfirmed_pos_screens) %>% distinct(pat_id)
df_baseline_qualified        <- df_baseline_evals %>% anti_join(df_total_baseline_exclusions, by = "pat_id")

df_baseline_history <- df_events_raw %>%
  mutate(pat_id = as.character(patient_number), event_date_clean = as.Date(event_date)) %>%
  left_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(event_date_clean < listing_date) %>%
  group_by(pat_id) %>%
  summarize(
    has_prior_transplant = if_else(any(event_code == "OTX"), 1, 0), 
    has_prior_pregnancy  = if_else(any(event_code == "PRG"), 1, 0),
    .groups = "drop"
  )

df_waitlist_endpoints <- df_events_raw %>%
  mutate(pat_id = as.character(patient_number), event_date_clean = as.Date(event_date)) %>%
  left_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(!is.na(event_date_clean) & event_date_clean >= listing_date & event_code %in% c("OTX", "DEA")) %>%
  group_by(pat_id) %>%
  summarize(removal_date = min(event_date_clean, na.rm = TRUE), .groups = "drop")

study_freeze_date <- max(c(df_patients_clean$listing_date, 
                           df_waitlist_endpoints$removal_date, 
                           df_last_lab_eval$last_eval_date), na.rm = TRUE)

df_incident_events <- df_baseline_qualified %>%
  inner_join(df_all_confirmed_positives, by = "pat_id") %>%
  inner_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(confirmed_pos_date > (listing_date + baseline_intake_cutoff_days)) %>%
  group_by(pat_id) %>%
  summarize(first_pos_date = min(confirmed_pos_date, na.rm = TRUE), .groups = "drop")

df_master_prelim <- df_patients_clean %>%
  inner_join(df_baseline_qualified, by = "pat_id") %>%
  left_join(df_baseline_history, by = "pat_id") %>%
  left_join(df_waitlist_endpoints, by = "pat_id") %>%
  left_join(df_incident_events, by = "pat_id") %>%
  left_join(df_last_lab_eval, by = "pat_id") %>%
  mutate(
    has_prior_transplant = coalesce(has_prior_transplant, 0),
    has_prior_pregnancy  = coalesce(has_prior_pregnancy, 0),
    sex_pregnancy_status = case_when(
      sex == "M"                            ~ "Male",
      sex == "F" & has_prior_pregnancy == 1 ~ "Female_Prior_Pregnancy",
      sex == "F" & has_prior_pregnancy == 0 ~ "Female_Nulligravid"
    ),
    sex_pregnancy_status = factor(
      sex_pregnancy_status, 
      levels = c("Female_Nulligravid", "Female_Prior_Pregnancy", "Male")
    ),
    entry_date = pmax(listing_date, baseline_anchor_date),
    status     = if_else(!is.na(first_pos_date), 1, 0),
    final_date = case_when(
      status == 1                 ~ first_pos_date,
      last_eval_date > entry_date ~ pmin(last_eval_date, coalesce(removal_date, study_freeze_date)),
      TRUE                        ~ entry_date
    ),
    end_day = as.numeric(difftime(final_date, entry_date, units = "days"))
  ) %>%
  filter(!is.na(end_day) & end_day > 0)

prod_col_name <- names(df_trans_raw)[
  str_detect(names(df_trans_raw), regex("product|prod|component|type|desc", ignore_case = TRUE))
][1]
if (is.na(prod_col_name)) prod_col_name <- "PRODUCT"

excluded_trans_recipients <- df_trans_raw %>%
  mutate(
    pat_id     = as.character(`Pat ID`),
    trans_date = as.Date(BLOOD_START_INSTANT),
    raw_prod   = as.character(.data[[prod_col_name]]),
    is_plt     = str_detect(raw_prod, regex("PLT|Platelet|Apheresis|Pooled", ignore_case = TRUE)),
    is_cryo    = str_detect(raw_prod, regex("CRYO|Cryoprecipitate", ignore_case = TRUE))
  ) %>%
  filter(is_plt | is_cryo) %>%
  inner_join(select(df_master_prelim, pat_id, entry_date, end_day), by = "pat_id") %>%
  mutate(raw_trans_day = as.numeric(difftime(trans_date, entry_date, units = "days"))) %>%
  filter(!is.na(raw_trans_day) & raw_trans_day >= 0 & raw_trans_day < end_day) %>%
  pull(pat_id) %>%
  unique()

df_master_base <- df_master_prelim %>%
  filter(!pat_id %in% excluded_trans_recipients)

df_trans_clean <- df_trans_raw %>%
  mutate(pat_id = as.character(`Pat ID`), trans_date = as.Date(BLOOD_START_INSTANT), units_transfused = 1) %>%
  filter(!pat_id %in% excluded_trans_recipients) %>%
  inner_join(select(df_master_base, pat_id, entry_date, end_day), by = "pat_id") %>%
  mutate(raw_trans_day = as.numeric(difftime(trans_date, entry_date, units = "days"))) %>%
  filter(!is.na(raw_trans_day) & raw_trans_day >= 0 & raw_trans_day < end_day) %>%
  group_by(pat_id, trans_day = raw_trans_day) %>%
  summarize(units_transfused = sum(units_transfused), .groups = "drop_last") %>%
  arrange(trans_day) %>%
  mutate(cum_units = cumsum(units_transfused)) %>%
  ungroup() %>%
  select(pat_id, trans_day, cum_units)

df_trans_products <- df_trans_raw %>%
  mutate(
    pat_id     = as.character(`Pat ID`),
    trans_date = as.Date(BLOOD_START_INSTANT),
    raw_prod   = as.character(.data[[prod_col_name]]),
    prod_type  = case_when(
      str_detect(raw_prod, regex("RBC|Red Blood|PRBC", ignore_case = TRUE)) ~ "RBC",
      str_detect(raw_prod, regex("Plasma|FFP", ignore_case = TRUE))         ~ "PLASMA",
      TRUE                                                                  ~ "OTHER"
    )
  ) %>%
  filter(!pat_id %in% excluded_trans_recipients & prod_type %in% c("RBC", "PLASMA")) %>%
  inner_join(select(df_master_base, pat_id, entry_date, end_day), by = "pat_id") %>%
  mutate(raw_trans_day = as.numeric(difftime(trans_date, entry_date, units = "days"))) %>%
  filter(!is.na(raw_trans_day) & raw_trans_day >= 0 & raw_trans_day < end_day) %>%
  distinct(pat_id, trans_day = raw_trans_day, prod_type) %>%
  mutate(exposed_flag = 1) %>%
  pivot_wider(names_from = prod_type, values_from = exposed_flag, values_fill = 0) %>%
  arrange(pat_id, trans_day)

for (p in c("RBC", "PLASMA")) {
  if (!p %in% colnames(df_trans_products)) df_trans_products[[p]] <- 0
}

df_trans_products <- df_trans_products %>%
  group_by(pat_id) %>%
  mutate(
    has_rbc    = as.integer(cummax(RBC) >= 1),
    has_plasma = as.integer(cummax(PLASMA) >= 1)
  ) %>%
  ungroup() %>%
  select(pat_id, trans_day, has_rbc, has_plasma)

message("[00_build_cohort] Writing data artifacts to ", output_dir, "...")
saveRDS(df_master_base,    file.path(output_dir, "df_master_base.rds"))
saveRDS(df_trans_clean,    file.path(output_dir, "df_trans_clean.rds"))
saveRDS(df_trans_products, file.path(output_dir, "df_trans_products.rds"))
saveRDS(df_screen_clean,   file.path(output_dir, "df_screen_clean.rds"))

n_hsc_dropped   <- df_demo_raw %>% filter(str_trim(description) %in% c("HSC Donor", "HSC Recipient", "HSC/Kidney Recipient")) %>% nrow()
n_other_dropped <- df_demo_raw %>% filter(!str_trim(description) %in% c("HSC Donor", "HSC Recipient", "HSC/Kidney Recipient")) %>% filter(str_trim(description) %in% descriptions_to_exclude | !sex %in% c("M", "F")) %>% nrow()

supp_table_1_attrition <- tibble(
  `Analytical Stage / Step` = c(
    "1. Registry Population (All Rows Ingested)",
    "2. Hematopoietic Stem Cell Exclusions (HSC Donors & Recipients)",
    "3. Other Clinical Exclusions (Living/Altruistic Donors, Non-Coded Sex)",
    "4. Eligible Solid Organ Candidates (Analytical Spine)",
    "5. Pre-Existing Sensitization Exclusion (Confirmed Positive <= Listing + 7 Days)",
    "6. Unconfirmed Positive Screen Exclusion (POS Screen without SAB in 30 Days)",
    "7. Risk Pool with Qualified Baseline Evaluation",
    "8. Dropped Due to Zero Longitudinal Surveillance Past Baseline Episode",
    "9. Pre-Transfusion Exclusion Solid Organ Cohort",
    "10. Platelet & Cryoprecipitate Recipient Exclusion",
    "11. Final Analyzed Cohort (Approach B: Verified Negative or Screen Negative)",
    "   - Event: Confirmed Incident Sensitization",
    "   - Event-Free: Right-Censored at Removal / Last Screen"
  ),
  `Candidate Count (N)` = c(
    nrow(df_demo_raw),
    n_hsc_dropped,
    n_other_dropped,
    nrow(df_patients_clean),
    nrow(df_baseline_sensitized),
    nrow(df_unconfirmed_pos_screens),
    nrow(df_baseline_qualified),
    nrow(df_baseline_qualified) - nrow(df_master_prelim),
    nrow(df_master_prelim),
    length(excluded_trans_recipients),
    nrow(df_master_base),
    sum(df_master_base$status == 1),
    sum(df_master_base$status == 0)
  )
) %>%
  mutate(
    `% of Eligible Spine` = case_when(
      row_number() == 4 ~ "100.0%",
      row_number() >= 5 ~ sprintf("%.1f%%", (`Candidate Count (N)` / `Candidate Count (N)`[4]) * 100),
      TRUE ~ "-"
    )
  )

write_csv(supp_table_1_attrition, file.path(table_dir, "Supplemental_Table_1_Cohort_Attrition.csv"))
saveRDS(supp_table_1_attrition,   file.path(output_dir, "consort_metrics.rds"))