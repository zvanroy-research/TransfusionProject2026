# ==============================================================================
# SUPPLEMENTAL TABLE S1: COHORT ORGAN/RECIPIENT DESCRIPTION BREAKDOWN
# Cohort: Dynamic Analysis Cohort (N = 5,435)
# ==============================================================================

library(gtsummary)
library(tidyverse)

tbl_supp1 <- df_table1_data %>%
  select(description, trans_tier) %>%
  tbl_summary(
    by = trans_tier,
    label = list(description ~ "Candidate Category"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = all_categorical() ~ c(0, 1),
    missing = "no"
  ) %>%
  add_overall(last = FALSE, col_label = "**Overall Cohort**\n(N = 5,435)") %>%
  bold_labels() %>%
  modify_header(
    label  = "**Organ / Recipient Category**",
    stat_0 = "**Overall**\nN = {N}",
    stat_1 = "**0 Units**\nN = {n}",
    stat_2 = "**1–4 Units**\nN = {n}",
    stat_3 = "**≥ 5 Units**\nN = {n}"
  ) %>%
  modify_footnote(
    all_stat_cols() ~ "Data presented as n (%)"
  )

# Export directly to Word .docx for Supplemental File
tbl_supp1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_S1_Candidate_Categories.docx")

# ==============================================================================
# SUPPLEMENTAL TABLE S2: COHORT ETHNICITY DEMOGRAPHICS (NUMBERS ONLY)
# Cohort: Dynamic Analysis Cohort (N = 5,435)
# ==============================================================================

library(readxl)
library(tidyverse)
library(gtsummary)
library(scales)

# ------------------------------------------------------------------------------
# 1. CLEAN ETHNICITY DATA FROM DEMOGRAPHICS
# ------------------------------------------------------------------------------
df_patients_clean_eth <- df_demo_raw %>%
  filter(!str_trim(description) %in% str_trim(descriptions_to_exclude)) %>%
  filter(sex %in% c("M", "F")) %>%
  mutate(
    pat_id         = as.character(patient_number),
    listing_date   = as.Date(as.numeric(create_date), origin = "1899-12-30"),
    birth_date     = as.Date(as.numeric(dob), origin = "1899-12-30"),
    age_at_listing = as.numeric(difftime(listing_date, birth_date, units = "weeks")) / 52.25,
    sex            = factor(sex, levels = c("F", "M")),
    
    # Clean race/ethnicity variable using `ethnic_code`
    race_ethnicity = if_else(
      is.na(ethnic_code) | str_trim(as.character(ethnic_code)) == "", 
      "Unknown / Unreported", 
      str_trim(as.character(ethnic_code))
    )
  ) %>%
  select(pat_id, listing_date, age_at_listing, sex, race_ethnicity)

# ------------------------------------------------------------------------------
# 2. LINK TO DYNAMIC COHORT (N = 5,435) & TRANSFUSION TIERS
# ------------------------------------------------------------------------------
df_cohort_eth_binned <- df_master_base %>%
  select(pat_id) %>%
  left_join(df_patients_clean_eth, by = "pat_id") %>%
  left_join(df_patient_max_trans, by = "pat_id") %>%
  mutate(
    total_units = coalesce(total_units, 0),
    trans_tier  = case_when(
      total_units == 0 ~ "0 Units",
      total_units >= 1 & total_units <= 4 ~ "1-4 Units",
      TRUE ~ "≥ 5 Units"
    ),
    trans_tier = factor(trans_tier, levels = c("0 Units", "1-4 Units", "≥ 5 Units"))
  )

# ------------------------------------------------------------------------------
# 3. GENERATE GTSUMMARY TABLE S2 (NO P-VALUES)
# ------------------------------------------------------------------------------
tbl_supp2_ethnicity <- df_cohort_eth_binned %>%
  select(race_ethnicity, trans_tier) %>%
  tbl_summary(
    by = trans_tier,
    label = list(race_ethnicity ~ "Race / Ethnicity"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = all_categorical() ~ c(0, 1),
    missing = "no"
  ) %>%
  add_overall(last = FALSE, col_label = "**Overall Cohort**\n(N = 5,435)") %>%
  bold_labels() %>%
  modify_header(
    label  = "**Race / Ethnicity Category**",
    stat_0 = "**Overall**\nN = {N}",
    stat_1 = "**0 Units**\nN = {n}",
    stat_2 = "**1–4 Units**\nN = {n}",
    stat_3 = "**≥ 5 Units**\nN = {n}"
  ) %>%
  modify_footnote(
    all_stat_cols() ~ "Data presented as n (%)"
  )

# Print in RStudio Viewer
print(tbl_supp2_ethnicity)

# ------------------------------------------------------------------------------
# 4. EXPORT WORD DOC & FLAT CSV
# ------------------------------------------------------------------------------
tbl_supp2_ethnicity %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_S2_Ethnicity_Demographics.docx")

df_eth_flat <- df_cohort_eth_binned %>%
  group_by(race_ethnicity, trans_tier) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(trans_tier) %>%
  mutate(pct = (n / sum(n)) * 100) %>%
  mutate(cell_label = paste0(comma(n), " (", sprintf("%.1f", pct), "%)")) %>%
  select(-n, -pct) %>%
  pivot_wider(names_from = trans_tier, values_from = cell_label, values_fill = "0 (0.0%)")

write.csv(df_eth_flat, "Table_S2_Race_Ethnicity_Flat.csv", row.names = FALSE)

# Export Word doc
tbl_supp2_ethnicity %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_S2_Ethnicity.docx")