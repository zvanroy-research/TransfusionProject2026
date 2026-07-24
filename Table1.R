# ==============================================================================
# RE-BINNED TABLE 1: BASELINE DEMOGRAPHICS & CLINICAL CHARACTERISTICS
# Option A: Active Waitlist Transfusions (0 Units N = 3,956)
# Outputs: Table1_Baseline_Characteristics_OptionA.docx
# ==============================================================================

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)

# ------------------------------------------------------------------------------
# 1. RE-BIN TRANSFUSION TIERS BASED ON ACTIVE WAITLIST LOGS
# ------------------------------------------------------------------------------
df_active_waitlist_units <- df_trans_cohort %>%
  group_by(pat_id) %>%
  summarize(active_units = n(), .groups = "drop")

df_table1_prep <- df_master_base %>%
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

# ------------------------------------------------------------------------------
# 2. BUILD GTSUMMARY TABLE WITH EXACT COLUMN NAMES
# ------------------------------------------------------------------------------
table1_gt <- df_table1_prep %>%
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
# 3. EXPORT TO WORD DOCUMENT (.DOCX)
# ------------------------------------------------------------------------------
ft_table1 <- as_flex_table(table1_gt) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bg(bg = "#F7F9FA", part = "header") %>%
  align(j = 2:6, align = "center", part = "all")

doc_table1 <- read_docx() %>%
  body_add_par("Table 1. Baseline Demographics and Clinical Characteristics Stratified by Transfusion Burden", style = "heading 1") %>%
  body_add_par("Dynamic Analysis Cohort (N = 5,435)", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft_table1) %>%
  body_add_par("") %>%
  body_add_par("1 Data presented as Median (IQR).", style = "Normal") %>%
  body_add_par("2 Data presented as n (%).", style = "Normal") %>%
  body_add_par("3 Evaluated across bins via Kruskal-Wallis rank sum test.", style = "Normal") %>%
  body_add_par("4 Evaluated across bins via Pearson's Chi-squared test.", style = "Normal") %>%
  body_add_par("*Note: Tiers reflect cumulative transfusion events accrued strictly during active waitlist follow-up (trans_day >= 0 to end_day). Candidates with recorded transfusion exposures occurring strictly prior to initial waitlist placement or post-delisting (n = 513) accrued 0 active waitlist units and are categorized under the 0 Units baseline tier. Primary cohort excludes baseline pre-existing sensitizations under dynamic model MFI thresholds.", style = "Normal")

print(doc_table1, target = "Table1_Baseline_Characteristics_OptionA.docx")

cat("\n--- SUCCESS: Table 1 generated cleanly without column errors! ---\n")