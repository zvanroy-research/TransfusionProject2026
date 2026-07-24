# ==============================================================================
# SUPPLEMENTAL TABLE S3: BLOOD PRODUCT ADMINISTRATION & DIVERSITY SUMMARY
# Cohort: Dynamic Analysis Cohort (N = 5,435)
# Transfused Cohort (N = 1,992)
# Outputs: Publication-ready Word Document (.docx) & Flat CSV Files
# ==============================================================================

library(tidyverse)
library(flextable)
library(officer)
library(scales)

# ------------------------------------------------------------------------------
# 1. FILTER TRANSFUSION LOG TO DYNAMIC COHORT (N = 5,435)
# ------------------------------------------------------------------------------
# Restrict transfusions to the waitlist observation window (trans_day >= 0 & <= end_day)
df_trans_cohort <- df_trans_raw %>%
  mutate(
    pat_id        = as.character(`Pat ID`),
    trans_date    = as.Date(BLOOD_START_INSTANT),
    blood_product = str_trim(as.character(BLOOD_PRODUCT))
  ) %>%
  inner_join(
    select(df_master_base, pat_id, listing_date, end_day), 
    by = "pat_id"
  ) %>%
  mutate(trans_day = as.numeric(difftime(trans_date, listing_date, units = "days"))) %>%
  filter(!is.na(trans_day) & trans_day >= 0 & trans_day <= end_day)

# Store total cohort denominators
n_cohort_total     <- nrow(df_master_base)             # N = 5,435
n_transfused_total <- n_distinct(df_trans_cohort$pat_id) # N = 1,992

# ------------------------------------------------------------------------------
# 2. TABLE S3A: OVERALL PRODUCT VOLUMES & CANDIDATE EXPOSURE
# ------------------------------------------------------------------------------
df_product_totals <- df_trans_cohort %>%
  group_by(blood_product) %>%
  summarize(
    total_units_administered  = n(),
    unique_patients_receiving = n_distinct(pat_id),
    .groups = "drop"
  ) %>%
  arrange(desc(total_units_administered)) %>%
  mutate(
    pct_of_transfused = (unique_patients_receiving / n_transfused_total) * 100,
    pct_of_cohort     = (unique_patients_receiving / n_cohort_total) * 100,
    patients_display  = paste0(comma(unique_patients_receiving), " (", sprintf("%.1f", pct_of_cohort), "%)")
  ) %>%
  select(
    `Blood Product Type`                 = blood_product,
    `Total Units Administered`           = total_units_administered,
    `Candidates Receiving, n (% Cohort)` = patients_display
  )

# ------------------------------------------------------------------------------
# 3. TABLE S3B: MUTUALLY EXCLUSIVE PRODUCT DIVERSITY BREAKDOWN
# ------------------------------------------------------------------------------
# Calculate distinct product types per patient
df_distinct_types_per_patient <- df_trans_cohort %>%
  group_by(pat_id) %>%
  summarize(distinct_types = n_distinct(blood_product), .groups = "drop")

# Merge back with master cohort (includes 0 distinct types for unexposed)
df_all_patient_type_counts <- df_master_base %>%
  select(pat_id) %>%
  left_join(df_distinct_types_per_patient, by = "pat_id") %>%
  mutate(distinct_types = coalesce(distinct_types, 0))

# Construct mutually exclusive categories summing cleanly to 1,992
df_multi_product_summary <- tibble(
  `Product Diversity Category` = c(
    "1 Product Type Only",
    "2 Distinct Product Types",
    "3 Distinct Product Types",
    "All 4 Distinct Product Types",
    "Total Transfused Candidates"
  ),
  `Candidate Count (N)` = c(
    sum(df_all_patient_type_counts$distinct_types == 1),
    sum(df_all_patient_type_counts$distinct_types == 2),
    sum(df_all_patient_type_counts$distinct_types == 3),
    sum(df_all_patient_type_counts$distinct_types == 4),
    n_transfused_total
  )
) %>%
  mutate(
    `% of Total Cohort (N = 5,435)`         = sprintf("%.2f%%", (`Candidate Count (N)` / n_cohort_total) * 100),
    `% of Transfused Candidates (N = 1,992)` = sprintf("%.2f%%", (`Candidate Count (N)` / n_transfused_total) * 100)
  )

# ------------------------------------------------------------------------------
# 4. FORMAT WITH FLEXTABLE (PUBLICATION STYLE)
# ------------------------------------------------------------------------------
# Style Table S3A
ft_product_totals <- flextable(df_product_totals) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bg(bg = "#F7F9FA", part = "header") %>%
  align(j = 2:3, align = "center", part = "all") %>%
  colformat_double(j = 2, digits = 0, big.mark = ",")

# Style Table S3B
ft_multi_product <- flextable(df_multi_product_summary) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = 5, part = "body") %>% # Bold total row
  bg(bg = "#F7F9FA", part = "header") %>%
  align(j = 2:4, align = "center", part = "all") %>%
  colformat_double(j = 2, digits = 0, big.mark = ",")

# ------------------------------------------------------------------------------
# 5. EXPORT TO WORD DOCUMENT (.DOCX) AND CSV
# ------------------------------------------------------------------------------
doc <- read_docx() %>%
  body_add_par("Supplemental Summary: Blood Product Administration Dynamics", style = "heading 1") %>%
  body_add_par("Dynamic Analysis Cohort (N = 5,435)", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_par("Table S3A. Overall Volume and Candidate Exposure by Blood Product Type", style = "heading 3") %>%
  body_add_flextable(ft_product_totals) %>%
  body_add_par("") %>%
  body_add_par("Table S3B. Mutually Exclusive Candidate Diversity of Transfused Blood Product Types", style = "heading 3") %>%
  body_add_flextable(ft_multi_product)

# Save Word Doc
print(doc, target = "Table_S3_Blood_Product_Summary.docx")

# Save CSV Files
write.csv(df_product_totals, "Table_S3A_Product_Totals.csv", row.names = FALSE)
write.csv(df_multi_product_summary, "Table_S3B_Multi_Product_Thresholds.csv", row.names = FALSE)

cat("\n--- SUCCESS: Table_S3_Blood_Product_Summary.docx generated successfully! ---\n")
