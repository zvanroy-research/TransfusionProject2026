# ==============================================================================
# RE-BINNED TABLE 1: CLINICAL CUTOFF MODEL (OPTION A)
# Re-binned strictly by active waitlist transfusion exposure
# Institutional Floors: 3000 MFI Class I / 1500 MFI Class II
# Outputs: Table1_Clinical_Model_OptionA.docx & Flat CSV
# ==============================================================================

library(gtsummary)
library(flextable)
library(officer)

# ------------------------------------------------------------------------------
# 1. RE-BIN TRANSFUSION TIERS BASED ON ACTIVE WAITLIST LOGS
# ------------------------------------------------------------------------------
# Filter transfusion log to active waitlist observation window
df_trans_cohort_active <- df_trans_raw %>%
  mutate(
    pat_id     = as.character(`Pat ID`),
    trans_date = as.Date(BLOOD_START_INSTANT)
  ) %>%
  inner_join(
    select(df_master_base, pat_id, listing_date, end_day), 
    by = "pat_id"
  ) %>%
  mutate(trans_day = as.numeric(difftime(trans_date, listing_date, units = "days"))) %>%
  filter(!is.na(trans_day) & trans_day >= 0 & trans_day <= end_day)

# Count active waitlist units per patient
df_active_waitlist_units <- df_trans_cohort_active %>%
  group_by(pat_id) %>%
  summarize(active_units = n(), .groups = "drop")

# Merge with master clinical cohort (df_master_base) to re-bin exposure tiers
df_table1_clinical_prep <- df_master_base %>%
  left_join(df_active_waitlist_units, by = "pat_id") %>%
  mutate(
    active_units = coalesce(active_units, 0),
    trans_tier_optionA = case_when(
      active_units == 0 ~ "0 Units",
      active_units >= 1 & active_units <= 4 ~ "1-4 Units",
      active_units >= 5 ~ "≥ 5 Units"
    ),
    trans_tier_optionA = factor(trans_tier_optionA, levels = c("0 Units", "1-4 Units", "≥ 5 Units"))
  )

# Print breakdown in console to verify
cat("\n--- CLINICAL MODEL TABLE 1 COLUMN COUNTS (OPTION A) ---\n")
print(table(df_table1_clinical_prep$trans_tier_optionA, useNA = "always"))

# ------------------------------------------------------------------------------
# 2. GENERATE GTSUMMARY TABLE MATCHING MANUSCRIPT FORMATTING
# ------------------------------------------------------------------------------
table1_clinical_gt <- df_table1_clinical_prep %>%
  select(
    trans_tier_optionA,
    `Age at Listing (years)` = age_at_listing,
    `Sex`                    = sex,
    `History of Pregnancy`   = has_prior_pregnancy,
    `History of Transplant`  = has_prior_transplant
  ) %>%
  tbl_summary(
    by = trans_tier_optionA,
    statistic = list(
      all_continuous()  ~ "{median} ({p25},{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ c(1, 1, 1),
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%
  add_overall(last = FALSE, col_label = "**Overall**\nN = {N}") %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = function(x) style_pvalue(x, digits = 3)
  ) %>%
  bold_labels() %>%
  modify_header(
    label = "**Characteristic**",
    stat_1 = "**0 Units**\nN = {n}",
    stat_2 = "**1-4 Units**\nN = {n}",
    stat_3 = "**≥ 5 Units**\nN = {n}",
    p.value = "**p-value**"
  )

# ------------------------------------------------------------------------------
# 3. CONVERT TO FLEXTABLE & EXPORT TO WORD (.DOCX)
# ------------------------------------------------------------------------------
ft_table1_clinical <- as_flex_table(table1_clinical_gt) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bg(bg = "#F7F9FA", part = "header") %>%
  align(j = 2:6, align = "center", part = "all")

doc_table1_clinical <- read_docx() %>%
  body_add_par("Table 1. Baseline Demographics and Clinical Characteristics Stratified by Transfusion Burden (Clinical Floor Model)", style = "heading 1") %>%
  body_add_par("Institutional Static Threshold Cohort (3,000 MFI MHC-I / 1,500 MFI MHC-II)", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft_table1_clinical) %>%
  body_add_par("") %>%
  body_add_par("1 Data presented as Median (IQR).", style = "Normal") %>%
  body_add_par("2 Data presented as n (%).", style = "Normal") %>%
  body_add_par("3 Evaluated across bins via Kruskal-Wallis rank sum test.", style = "Normal") %>%
  body_add_par("4 Evaluated across bins via Pearson's Chi-squared test.", style = "Normal") %>%
  body_add_par("*Note: Tiers reflect cumulative transfusion events accrued strictly during active waitlist follow-up (trans_day >= 0 to end_day). Candidates with recorded transfusion exposures occurring strictly prior to initial waitlist placement or post-delisting accrued 0 active waitlist units and are categorized under the 0 Units baseline tier. Primary cohort excludes baseline pre-existing sensitizations under institutional clinical floors.", style = "Normal")

# Save Word Doc
print(doc_table1_clinical, target = "Table1_Clinical_Model_OptionA.docx")

# Save flat CSV version
write.csv(as_tibble(table1_clinical_gt), "Table1_Clinical_Model_OptionA.csv", row.names = FALSE)

cat("\n--- SUCCESS: Table1_Clinical_Model_OptionA.docx generated successfully! ---\n")