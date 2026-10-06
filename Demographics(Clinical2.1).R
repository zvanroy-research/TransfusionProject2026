# ==============================================================================
# script 1: demographics
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

root_dir   <- "/Users/zvr/Downloads/Clinical2.0"
output_dir <- file.path(root_dir, "processed_data")
table_dir  <- file.path(root_dir, "tables_output")

df_master_base    <- readRDS(file.path(output_dir, "df_master_base.rds"))
df_trans_clean    <- readRDS(file.path(output_dir, "df_trans_clean.rds"))
df_trans_products <- readRDS(file.path(output_dir, "df_trans_products.rds"))

if (!"raw_desc" %in% names(df_master_base))  df_master_base$raw_desc  <- as.character(df_master_base$organ_group)
if (!"ethnicity" %in% names(df_master_base)) df_master_base$ethnicity <- "Unknown/Unspecified"

transfused_ids_all <- unique(df_trans_clean$pat_id)

df_cohort_profile <- df_master_base %>%
  mutate(
    trans_status = if_else(pat_id %in% transfused_ids_all, "Transfused (1+ Units)", "Never Transfused"),
    trans_status = factor(trans_status, levels = c("Never Transfused", "Transfused (1+ Units)")),
    # Explicitly enforce the new 4-tier taxonomy levels to guarantee correct display order in tables
    organ_group  = factor(organ_group, levels = c("Isolated Kidney", "Isolated Heart", "Isolated Lung", "Multi-Organ & Other"))
  )

format_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", n, (n / total) * 100)
}

format_cont_median <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f [%.1f–%.1f]", q[2], q[1], q[3])
}

# --- table 1---

df_cohort_profile <- df_cohort_profile %>%
  mutate(followup_days = end_day)

df_cohort_profile <- df_cohort_profile %>%
  mutate(
    eth_code = str_to_upper(str_trim(coalesce(as.character(ethnicity), "UNK"))),
    clean_eth = case_when(
      eth_code == "CAU" ~ "White / Caucasian",
      eth_code == "BLK" ~ "Black / African American",
      eth_code == "HIS" ~ "Hispanic / Latino",
      eth_code == "AIN" ~ "American Indian / Alaska Native",
      eth_code == "ASI" ~ "Asian",
      TRUE              ~ "Other / Unknown"
    ),
    clean_eth = factor(clean_eth, levels = c(
      "White / Caucasian", "Black / African American", "Hispanic / Latino", 
      "Asian", "American Indian / Alaska Native", "Other / Unknown"
    )),
    clean_eth = droplevels(clean_eth)
  )

n_tot <- nrow(df_cohort_profile)
d_un  <- filter(df_cohort_profile, trans_status == "Never Transfused")
d_tx  <- filter(df_cohort_profile, trans_status == "Transfused (1+ Units)")
n_un  <- nrow(d_un)
n_tx  <- nrow(d_tx)

format_p <- function(p) {
  if (is.null(p) || is.na(p) || is.nan(p)) return("-")
  if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
}

rows <- list()

# Total Candidates
rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Total Candidates",
  `Overall (N = 2,442)`           = as.character(n_tot),
  `Never Transfused (n = 2,184)`   = as.character(n_un),
  `Transfused (1+ Units) (n = 258)` = as.character(n_tx),
  `p-value`                       = "-"
)

# Follow-up Time
p_fu <- wilcox.test(d_un$followup_days, d_tx$followup_days)$p.value
rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Follow-up Surveillance, Median [IQR], days",
  `Overall (N = 2,442)`           = format_cont_median(df_cohort_profile$followup_days),
  `Never Transfused (n = 2,184)`   = format_cont_median(d_un$followup_days),
  `Transfused (1+ Units) (n = 258)` = format_cont_median(d_tx$followup_days),
  `p-value`                       = format_p(p_fu)
)

# Age at Listing
p_age <- wilcox.test(d_un$age_at_listing, d_tx$age_at_listing)$p.value
rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Age at Listing, Median [IQR], years",
  `Overall (N = 2,442)`           = format_cont_median(df_cohort_profile$age_at_listing),
  `Never Transfused (n = 2,184)`   = format_cont_median(d_un$age_at_listing),
  `Transfused (1+ Units) (n = 258)` = format_cont_median(d_tx$age_at_listing),
  `p-value`                       = format_p(p_age)
)

# Sex
p_sex <- chisq.test(table(df_cohort_profile$sex, df_cohort_profile$trans_status))$p.value
for (sx in levels(df_cohort_profile$sex)) {
  rows[[length(rows) + 1]] <- tibble(
    Characteristic                  = paste0("  Sex: ", if_else(sx == "M", "Male", "Female")),
    `Overall (N = 2,442)`           = format_n_pct(sum(df_cohort_profile$sex == sx), n_tot),
    `Never Transfused (n = 2,184)`   = format_n_pct(sum(d_un$sex == sx), n_un),
    `Transfused (1+ Units) (n = 258)` = format_n_pct(sum(d_tx$sex == sx), n_tx),
    `p-value`                       = if_else(sx == "M", format_p(p_sex), "")
  )
}

# Pregnancy Status
n_fem_tot <- sum(df_cohort_profile$sex == "F")
n_fem_un  <- sum(d_un$sex == "F")
n_fem_tx  <- sum(d_tx$sex == "F")
p_prg     <- chisq.test(table(filter(df_cohort_profile, sex == "F")$has_prior_pregnancy, 
                              filter(df_cohort_profile, sex == "F")$trans_status))$p.value

rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "  Prior Pregnancy (% of female candidates)",
  `Overall (N = 2,442)`           = format_n_pct(sum(df_cohort_profile$sex == "F" & df_cohort_profile$has_prior_pregnancy == 1), n_fem_tot),
  `Never Transfused (n = 2,184)`   = format_n_pct(sum(d_un$sex == "F" & d_un$has_prior_pregnancy == 1), n_fem_un),
  `Transfused (1+ Units) (n = 258)` = format_n_pct(sum(d_tx$sex == "F" & d_tx$has_prior_pregnancy == 1), n_fem_tx),
  `p-value`                       = format_p(p_prg)
)

# Prior Transplant
p_txp <- chisq.test(table(df_cohort_profile$has_prior_transplant, df_cohort_profile$trans_status))$p.value
rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Prior Solid Organ Transplant",
  `Overall (N = 2,442)`           = format_n_pct(sum(df_cohort_profile$has_prior_transplant == 1), n_tot),
  `Never Transfused (n = 2,184)`   = format_n_pct(sum(d_un$has_prior_transplant == 1), n_un),
  `Transfused (1+ Units) (n = 258)` = format_n_pct(sum(d_tx$has_prior_transplant == 1), n_tx),
  `p-value`                       = format_p(p_txp)
)

# Race / Ethnicity
tab_eth <- table(df_cohort_profile$clean_eth, df_cohort_profile$trans_status)
p_eth   <- suppressWarnings(chisq.test(tab_eth, simulate.p.value = TRUE, B = 2000)$p.value)

rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Race / Ethnicity", 
  `Overall (N = 2,442)`           = "", 
  `Never Transfused (n = 2,184)`   = "", 
  `Transfused (1+ Units) (n = 258)` = "", 
  `p-value`                       = format_p(p_eth)
)
for (eth in levels(df_cohort_profile$clean_eth)) {
  n_sub_tot <- sum(df_cohort_profile$clean_eth == eth)
  if (n_sub_tot > 0) {
    rows[[length(rows) + 1]] <- tibble(
      Characteristic                  = paste0("  ", eth),
      `Overall (N = 2,442)`           = format_n_pct(n_sub_tot, n_tot),
      `Never Transfused (n = 2,184)`   = format_n_pct(sum(d_un$clean_eth == eth), n_un),
      `Transfused (1+ Units) (n = 258)` = format_n_pct(sum(d_tx$clean_eth == eth), n_tx),
      `p-value`                       = ""
    )
  }
}

# Organ Candidacy
tab_og <- table(df_cohort_profile$organ_group, df_cohort_profile$trans_status)
p_og   <- suppressWarnings(chisq.test(tab_og, simulate.p.value = TRUE, B = 2000)$p.value)

rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Organ Candidacy Category", 
  `Overall (N = 2,442)`           = "", 
  `Never Transfused (n = 2,184)`   = "", 
  `Transfused (1+ Units) (n = 258)` = "", 
  `p-value`                       = format_p(p_og)
)
for (og in levels(df_cohort_profile$organ_group)) {
  rows[[length(rows) + 1]] <- tibble(
    Characteristic                  = paste0("  ", og),
    `Overall (N = 2,442)`           = format_n_pct(sum(df_cohort_profile$organ_group == og), n_tot),
    `Never Transfused (n = 2,184)`   = format_n_pct(sum(d_un$organ_group == og), n_un),
    `Transfused (1+ Units) (n = 258)` = format_n_pct(sum(d_tx$organ_group == og), n_tx),
    `p-value`                       = ""
  )
}

# Incident Sensitization
p_stat <- chisq.test(table(df_cohort_profile$status, df_cohort_profile$trans_status))$p.value
rows[[length(rows) + 1]] <- tibble(
  Characteristic                  = "Incident Sensitization (Confirmed Event)",
  `Overall (N = 2,442)`           = format_n_pct(sum(df_cohort_profile$status == 1), n_tot),
  `Never Transfused (n = 2,184)`   = format_n_pct(sum(d_un$status == 1), n_un),
  `Transfused (1+ Units) (n = 258)` = format_n_pct(sum(d_tx$status == 1), n_tx),
  `p-value`                       = format_p(p_stat)
)

table_1 <- bind_rows(rows)
write_csv(table_1, file.path(table_dir, "Table_1_Demographics.csv"))
cat("\n--- UPDATED TABLE 1: DEMOGRAPHICS ---\n")
print(as.data.frame(table_1), row.names = FALSE)

# --- Sup table 2

df_tx_summary <- df_trans_clean %>%
  group_by(pat_id) %>%
  summarize(
    total_units  = max(cum_units), 
    first_tx_day = min(trans_day), 
    .groups      = "drop"
  )

df_prod_summary <- df_trans_products %>%
  group_by(pat_id) %>%
  summarize(
    received_rbc    = as.integer(any(has_rbc == 1)), 
    received_plasma = as.integer(any(has_plasma == 1)), 
    .groups         = "drop"
  )

df_tx_full <- df_tx_summary %>%
  left_join(df_prod_summary, by = "pat_id") %>%
  inner_join(
    df_master_base %>% select(pat_id, status, organ_group),
    by = "pat_id"
  ) %>%
  mutate(
    organ_group = factor(
      organ_group, 
      levels = c("Isolated Kidney", "Isolated Heart", "Isolated Lung", "Multi-Organ & Other")
    )
  )

fmt_med_iqr <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return("—")
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f [%.1f–%.1f]", q[2], q[1], q[3])
}

fmt_mean_sd <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return("—")
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_range <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return("—")
  sprintf("%d – %d", as.integer(min(x, na.rm = TRUE)), as.integer(max(x, na.rm = TRUE)))
}

fmt_count_pct <- function(n, total) {
  if (total == 0) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", n, (n / total) * 100)
}

build_tx_column <- function(df_sub) {
  n_sub <- nrow(df_sub)
  if (n_sub == 0) return(rep("—", 24))
  
  c(
    as.character(n_sub),
    "", # Cumulative Transfusion Burden header
    fmt_med_iqr(df_sub$total_units),
    fmt_mean_sd(df_sub$total_units),
    fmt_range(df_sub$total_units),
    "", # Dose header
    fmt_count_pct(sum(df_sub$total_units %in% 1:2), n_sub),
    fmt_count_pct(sum(df_sub$total_units %in% 3:5), n_sub),
    fmt_count_pct(sum(df_sub$total_units %in% 6:9), n_sub),
    fmt_count_pct(sum(df_sub$total_units >= 10), n_sub),
    "", # Transfusion Burden by Organ Group header
    fmt_med_iqr(df_sub$total_units[df_sub$organ_group == "Isolated Kidney"]),
    fmt_med_iqr(df_sub$total_units[df_sub$organ_group == "Isolated Heart"]),
    fmt_med_iqr(df_sub$total_units[df_sub$organ_group == "Isolated Lung"]),
    fmt_med_iqr(df_sub$total_units[df_sub$organ_group == "Multi-Organ & Other"]),
    "", # Product Types Transfused header
    fmt_count_pct(sum(df_sub$received_rbc == 1, na.rm = TRUE), n_sub),
    fmt_count_pct(sum(df_sub$received_plasma == 1, na.rm = TRUE), n_sub),
    fmt_count_pct(sum(df_sub$received_rbc == 1 & df_sub$received_plasma == 1, na.rm = TRUE), n_sub),
    "", # Temporal Exposure Distribution header
    fmt_med_iqr(df_sub$first_tx_day),
    fmt_count_pct(sum(df_sub$first_tx_day <= 30), n_sub),
    fmt_count_pct(sum(df_sub$first_tx_day > 30), n_sub)
  )
}

d_all  <- df_tx_full
d_sens <- filter(df_tx_full, status == 1)
d_uns  <- filter(df_tx_full, status == 0)

n_kidney <- sum(df_tx_full$organ_group == "Isolated Kidney")
n_heart  <- sum(df_tx_full$organ_group == "Isolated Heart")
n_lung   <- sum(df_tx_full$organ_group == "Isolated Lung")
n_multi  <- sum(df_tx_full$organ_group == "Multi-Organ & Other")

supp_table_2 <- tibble(
  `Transfusion Characteristic` = c(
    "Total Transfused Candidates",
    "Cumulative Transfusion Burden",
    "  Median Units [IQR]",
    "  Mean Units (SD)",
    "  Units Range (Min – Max)",
    "Dose",
    "  1–2 Units",
    "  3–5 Units",
    "  6–9 Units",
    "  ≥ 10 Units",
    "Transfusion Burden by Organ Group",
    sprintf("  Isolated Kidney (n = %d)", n_kidney),
    sprintf("  Isolated Heart (n = %d)", n_heart),
    sprintf("  Isolated Lung (n = %d)", n_lung),
    sprintf("  Multi-Organ & Other (n = %d)", n_multi),
    "Product Types Transfused",
    "  Red Blood Cells (RBC)",
    "  Plasma (FFP)",
    "  Both RBC and Plasma",
    "Temporal Exposure Distribution",
    "  Days to First Transfusion, Median [IQR]",
    "  First Transfusion ≤ 30 Days",
    "  First Transfusion > 30 Days"
  ),
  `Overall`      = build_tx_column(d_all),
  `Sensitized`   = build_tx_column(d_sens),
  `Unsensitized` = build_tx_column(d_uns)
)

write_csv(supp_table_2, file.path(table_dir, "Supplementary_Table_2_Transfusion_Characteristics.csv"))
cat("\n--- REVISED SUPPLEMENTARY TABLE 2 (4-TIER ORGAN TAXONOMY) ---\n")
print(as.data.frame(supp_table_2), row.names = FALSE)

# sup table 3

cat("\n--- RAW ETHNIC CODE DISTRIBUTION IN FINAL COHORT ---\n")
print(table(df_cohort_profile$ethnicity, useNA = "always"))

supp_table_3_ethnicity <- df_cohort_profile %>%
  mutate(
    code = str_to_upper(str_trim(coalesce(as.character(ethnicity), "UNK"))),
    clean_eth = case_when(
      code == "CAU" ~ "White / Caucasian",
      code == "BLK" ~ "Black / African American",
      code == "HIS" ~ "Hispanic / Latino",
      code == "AIN" ~ "American Indian / Alaska Native",
      code == "ASI" ~ "Asian",
      code %in% c("OTH", "UNK", "") ~ "Other / Unknown",
      TRUE ~ "Other / Unknown"
    ),
    clean_eth = factor(
      clean_eth, 
      levels = c(
        "White / Caucasian", 
        "Black / African American", 
        "Hispanic / Latino", 
        "Asian", 
        "American Indian / Alaska Native", 
        "Other / Unknown"
      )
    )
  ) %>%
  group_by(clean_eth) %>%
  summarize(
    `Total (N)`            = n(),
    `Never Transfused (n)` = sum(trans_status == "Never Transfused"),
    `Transfused (n)`       = sum(trans_status == "Transfused (1+ Units)"),
    `Sensitized (n)`       = sum(status == 1),
    `Sensitization Rate`   = sprintf("%.1f%%", 100 * mean(status == 1)),
    .groups                = "drop"
  ) %>%
  complete(
    clean_eth, 
    fill = list(
      `Total (N)` = 0, 
      `Never Transfused (n)` = 0, 
      `Transfused (n)` = 0, 
      `Sensitized (n)` = 0, 
      `Sensitization Rate` = "0.0%"
    )
  ) %>%
  filter(`Total (N)` > 0) %>%
  arrange(desc(`Total (N)`)) %>%
  mutate(
    `% of Cohort` = sprintf("%.1f%%", 100 * (`Total (N)` / sum(`Total (N)`)))
  ) %>%
  select(
    `Race / Ethnicity` = clean_eth,
    `Total (N)`,
    `% of Cohort`,
    `Never Transfused (n)`,
    `Transfused (n)`,
    `Sensitized (n)`,
    `Sensitization Rate`
  )

write_csv(supp_table_3_ethnicity, file.path(table_dir, "Supplementary_Table_3_Ethnicity.csv"))

cat("\n--- SUPPLEMENTARY TABLE 3: ETHNICITY IN COHORTS ---\n")
print(as.data.frame(supp_table_3_ethnicity), row.names = FALSE)

# --- Sup table 4

supp_table_4_organs <- df_cohort_profile %>%
  group_by(organ_group, raw_desc) %>%
  summarize(Candidates = n(), Sensitized = sum(status == 1), .groups = "drop") %>%
  group_by(organ_group) %>%
  mutate(subtotal_group = sum(Candidates), pct_within_group = sprintf("%.1f%%", 100 * (Candidates / subtotal_group))) %>%
  ungroup() %>%
  arrange(organ_group, desc(Candidates)) %>%
  mutate(
    pct_of_total_cohort = sprintf("%.1f%%", 100 * (Candidates / nrow(df_cohort_profile))),
    incidence_pct       = sprintf("%.1f%%", 100 * (Sensitized / Candidates))
  ) %>%
  select(
    `Analytical Organ Group`        = organ_group,
    `Granular Registry Description` = raw_desc,
    `Candidate Count (N)`           = Candidates,
    `Group Subtotal`                = subtotal_group,
    `% of Organ Category`           = pct_within_group,
    `% of Total Cohort`             = pct_of_total_cohort,
    `Incident Sensitization (n)`    = Sensitized,
    `Sensitization Incidence (%)`   = incidence_pct
  )

write_csv(supp_table_4_organs, file.path(table_dir, "Supplementary_Table_4_Expanded_Organs.csv"))