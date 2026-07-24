# ==============================================================================
# CONTINUOUS NON-LINEAR TRANSFUSION-ASSOCIATED RISK ENGINE
# Model: Time-Varying Cox Regression with Natural Cubic Splines (df = 3)
# Output: Adjusted Hazard Ratio Curve Relative to 0 Transfused Units
# ==============================================================================

library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(splines)
library(ggplot2)
library(scales)
library(cowplot)

# ------------------------------------------------------------------------------
# 1. FILE PATHS & DATA INGESTION
# ------------------------------------------------------------------------------
file_path <- "/Users/zvr/Downloads/Transfusion_Data.xlsx"
cut_path  <- "/Users/zvr/Downloads/Cuttoffs.xlsx"

df_demo_raw   <- read_excel(file_path, sheet = "Demographics")
df_trans_raw  <- read_excel(file_path, sheet = "Transfusion")
df_events_raw <- read_excel(file_path, sheet = "Patient Events")
df_mhc1_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-I")
df_mhc2_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-II")
df_cutoffs    <- read_excel(cut_path)

# ------------------------------------------------------------------------------
# 2. DEMOGRAPHIC CLEANING & EXCLUSIONS
# ------------------------------------------------------------------------------
descriptions_to_exclude <- c(
  "Kidney Living Donor", "HSC Donor", "Disease Association",
  "Altruistic Donor", "OTHER/RECIPIENT", "Liver Living Donor", "OTHER/DONOR"
)

df_patients_clean <- df_demo_raw %>%
  filter(!str_trim(description) %in% str_trim(descriptions_to_exclude)) %>%
  filter(sex %in% c("M", "F")) %>%
  mutate(
    pat_id         = as.character(patient_number),
    listing_date   = as.Date(as.numeric(create_date), origin = "1899-12-30"),
    birth_date     = as.Date(as.numeric(dob), origin = "1899-12-30"),
    age_at_listing = as.numeric(difftime(listing_date, birth_date, units = "weeks")) / 52.25,
    sex            = factor(sex, levels = c("F", "M"))
  ) %>%
  select(pat_id, listing_date, age_at_listing, sex)

# ------------------------------------------------------------------------------
# 3. LAB MATRIX & EVENT DEFINITIONS
# ------------------------------------------------------------------------------
df_first_positive_mfi <- bind_rows(
  df_mhc1_raw %>% 
    mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>%
    select(pat_id, draw_date_clean, any_of(df_cutoffs$Target)) %>%
    pivot_longer(cols = -c(pat_id, draw_date_clean), names_to = "Target", values_to = "MFI"),
  
  df_mhc2_raw %>% 
    mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>%
    select(pat_id, draw_date_clean, any_of(df_cutoffs$Target)) %>%
    pivot_longer(cols = -c(pat_id, draw_date_clean), names_to = "Target", values_to = "MFI")
) %>%
  left_join(df_cutoffs, by = "Target") %>%
  filter(MFI >= Cuttoff) %>%
  group_by(pat_id) %>%
  summarize(first_pos_date = min(draw_date_clean, na.rm = TRUE)) %>%
  ungroup()

df_last_lab_mfi <- bind_rows(
  df_mhc1_raw %>% mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>% select(pat_id, draw_date_clean),
  df_mhc2_raw %>% mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>% select(pat_id, draw_date_clean)
) %>%
  group_by(pat_id) %>%
  summarize(last_draw_date = max(draw_date_clean, na.rm = TRUE)) %>%
  ungroup()

df_baseline_history <- df_events_raw %>%
  mutate(pat_id = as.character(patient_number), event_date_clean = as.Date(event_date)) %>%
  left_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(event_date_clean < listing_date) %>%
  group_by(pat_id) %>%
  summarize(
    has_prior_transplant = if_else(any(event_code == "OTX"), 1, 0), 
    has_prior_pregnancy  = if_else(any(event_code == "PRG"), 1, 0)
  ) %>%
  ungroup()

df_waitlist_endpoints <- df_events_raw %>%
  mutate(pat_id = as.character(patient_number), event_date_clean = as.Date(event_date)) %>%
  left_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  filter(!is.na(event_date_clean) & event_date_clean >= listing_date & event_code %in% c("OTX", "DEA")) %>%
  group_by(pat_id) %>%
  summarize(removal_date = min(event_date_clean, na.rm = TRUE)) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 4. TRANSFUSION COUNTER & TIME-VARYING MERGE
# ------------------------------------------------------------------------------
df_trans_clean <- df_trans_raw %>%
  mutate(pat_id = as.character(`Pat ID`), trans_date = as.Date(BLOOD_START_INSTANT), units_transfused = 1) %>%
  left_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
  mutate(trans_day = as.numeric(difftime(trans_date, listing_date, units = "days"))) %>%
  filter(!is.na(trans_day) & trans_day >= 0) %>%
  arrange(pat_id, trans_day) %>%
  group_by(pat_id) %>%
  mutate(cum_units = cumsum(units_transfused)) %>%
  ungroup() %>%
  select(pat_id, trans_day, cum_units)

study_freeze_date <- max(c(df_patients_clean$listing_date, df_waitlist_endpoints$removal_date, df_last_lab_mfi$last_draw_date), na.rm = TRUE)

df_master_base <- df_patients_clean %>%
  left_join(df_baseline_history, by = "pat_id") %>%
  left_join(df_waitlist_endpoints, by = "pat_id") %>%
  left_join(df_first_positive_mfi, by = "pat_id") %>%
  left_join(df_last_lab_mfi, by = "pat_id") %>%
  mutate(
    has_prior_transplant = coalesce(has_prior_transplant, 0),
    has_prior_pregnancy  = coalesce(has_prior_pregnancy, 0),
    status               = if_else(!is.na(first_pos_date), 1, 0),
    last_draw_date_post  = if_else(last_draw_date > listing_date, last_draw_date, as.Date(NA)),
    final_date = case_when(
      !is.na(first_pos_date)      ~ first_pos_date, 
      !is.na(removal_date)        ~ removal_date, 
      !is.na(last_draw_date_post) ~ last_draw_date_post, 
      TRUE                        ~ study_freeze_date
    ),
    end_day = as.numeric(difftime(final_date, listing_date, units = "days"))
  ) %>%
  filter(is.na(first_pos_date) | as.numeric(difftime(first_pos_date, listing_date, units="days")) > 14) %>%
  mutate(end_day = if_else(end_day <= 0 & is.na(first_pos_date), as.numeric(difftime(coalesce(removal_date, study_freeze_date), listing_date, units="days")), end_day)) %>%
  filter(!is.na(end_day) & end_day > 0)

df_long_base <- tmerge(data1 = df_master_base, data2 = df_master_base, id = pat_id, tstart = 0, tstop = end_day)
df_long_base <- tmerge(data1 = df_long_base, data2 = df_master_base, id = pat_id, event = event(end_day, status))

df_long_spline <- tmerge(
  data1 = df_long_base,
  data2 = df_trans_clean,
  id = pat_id,
  cum_units = tdc(trans_day, cum_units)
) %>%
  mutate(cum_units = coalesce(cum_units, 0))

# Sample dataset of exposure values for the density rug plot (subsampled for speed/clarity)
df_rug_data <- df_long_spline %>% 
  filter(event == 1 | cum_units > 0) %>% 
  sample_n(min(n(), 1500))

# ------------------------------------------------------------------------------
# 5. FIT NATURAL CUBIC SPLINE MODEL (df = 3)
# ------------------------------------------------------------------------------
message("Fitting Natural Cubic Spline Cox Proportional Hazards Model...")

fit_spline <- coxph(
  Surv(tstart, tstop, event) ~ ns(cum_units, df = 3) + age_at_listing + sex + 
    has_prior_transplant + has_prior_pregnancy, 
  data = df_long_spline
)

# ------------------------------------------------------------------------------
# 5B. STATISTICAL TESTS FOR SPLINE FIT & NON-LINEARITY (FAIL-SAFE)
# ------------------------------------------------------------------------------
# Fit linear Cox model to compare against natural cubic spline
fit_linear <- coxph(
  Surv(tstart, tstop, event) ~ cum_units + age_at_listing + sex + 
    has_prior_transplant + has_prior_pregnancy, 
  data = df_long_spline
)

# Perform Likelihood Ratio Test comparing Linear vs. Spline
lrt_res <- anova(fit_linear, fit_spline, test = "Chisq")

# Extract p-value safely across survival package versions
p_nonlinearity <- na.omit(lrt_res[[grep("P|p", colnames(lrt_res))[1]]])[1]

# Format helper
format_pval_safe <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("N/A")
  if (p < 0.001) return("p < 0.001")
  return(sprintf("p = %.3f", p))
}

cat("\n======================================================================\n")
cat("                SPLINE MODEL STATISTICAL METRICS                      \n")
cat("======================================================================\n")
cat(sprintf("LRT p-value for Non-Linearity: %s\n", format_pval_safe(p_nonlinearity)))
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# 6. EXTRACT PREDICTED RELATIVE HAZARD RATIOS & 95% CIs
# ------------------------------------------------------------------------------
pred_grid <- data.frame(
  cum_units            = seq(0, 25, by = 0.1),
  age_at_listing       = 50,
  sex                  = factor("F", levels = c("F", "M")),
  has_prior_transplant = 0,
  has_prior_pregnancy  = 0
)

# Baseline reference prediction at 0 units
ref_pred  <- predict(fit_spline, newdata = transform(pred_grid[1,], cum_units = 0), type = "lp", se.fit = TRUE)
fit_preds <- predict(fit_spline, newdata = pred_grid, type = "lp", se.fit = TRUE)

df_spline_plot <- pred_grid %>%
  mutate(
    log_hr   = fit_preds$fit - ref_pred$fit,
    se       = fit_preds$se.fit,
    hr       = exp(log_hr),
    lower_ci = exp(log_hr - 1.96 * se),
    upper_ci = exp(log_hr + 1.96 * se)
  )

# ------------------------------------------------------------------------------
# 6B. EXTRACT PREDICTED HRS AT LANDMARK EXPOSURE THRESHOLDS
# ------------------------------------------------------------------------------
landmark_units <- c(1, 2, 3, 4, 5, 8, 10, 15, 20)

df_landmark_table <- df_spline_plot %>%
  filter(cum_units %in% landmark_units) %>%
  mutate(
    `Transfused Units` = sprintf("%d units", cum_units),
    `Adjusted HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", hr, lower_ci, upper_ci)
  ) %>%
  select(`Transfused Units`, `Adjusted HR (95% CI)`)

print(df_landmark_table)

# Export as CSV for Supplementary Material (Table S4)
write.csv(df_landmark_table, "Table_S4_Spline_Landmark_HRs.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. GENERATE PUBLICATION-READY SPLINE VISUALIZATION
# ------------------------------------------------------------------------------
p_spline <- ggplot(df_spline_plot, aes(x = cum_units, y = hr)) +
  # Reference line at HR = 1.0
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.5) +
  
  # 95% Confidence Interval Ribbon
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "#1F618D", alpha = 0.18) +
  
  # Main Spline Curve
  geom_line(color = "#1F618D", linewidth = 1.1) +
  
  # Baseline Reference Marker at (0, 1)
  annotate("point", x = 0, y = 1.0, color = "#2C3E50", size = 2.0) +
  
  # Exposure Density Rug Plot along X-axis
  geom_rug(data = df_rug_data, aes(x = cum_units), inherit.aes = FALSE, 
           sides = "b", color = "#2C3E50", alpha = 0.25, length = unit(0.02, "npc")) +
  
  # Axis Scales & Limits
  scale_x_continuous(
    breaks = seq(0, 25, by = 5), 
    limits = c(0, 25), 
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0.50, 2.00, by = 0.25), 
    limits = c(0.40, 2.10),
    expand = c(0.01, 0)
  ) +
  
  # Publication Labels
  labs(
    title = "Dose-Dependent Sensitization Risk by Transfused Units",
    subtitle = "Utilizing Dynamic Analysis, Relative to Unexposed",
    x = "Cumulative Blood Units Transfused",
    y = "Adjusted Hazard Ratio"
  ) +
  
  # Publication Theme
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(face = "plain", size = 9.5, hjust = 0.5, color = "#333333", margin = margin(b = 12)),
    axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 10, margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 9),
    axis.line = element_line(color = "black", linewidth = 0.6),
    plot.margin = margin(t = 12, r = 12, b = 10, l = 10)
  )

print(p_spline)

# ------------------------------------------------------------------------------
# 8. EXPORT HIGH-RESOLUTION SUBMISSION FILES
# ------------------------------------------------------------------------------
ggsave("Figure_3_Spline.pdf", p_spline, width = 7.0, height = 5.0, units = "in", device = cairo_pdf)
ggsave("Figure_3_Spline.tiff", p_spline, width = 7.0, height = 5.0, units = "in", dpi = 300, compression = "lzw")

message("Success! Figure 3 exported in high-resolution PDF & 300 DPI TIFF formats.")

# ==============================================================================
# EXPORT SUPPLEMENTARY TABLE S4: LANDMARK SPLINE HAZARD RATIOS
# Outputs: Table_S4_Spline_Landmark_HRs.docx & Table_S4_Spline_Landmark_HRs.csv
# ==============================================================================

library(flextable)
library(officer)

# ------------------------------------------------------------------------------
# 1. PREPARE LANDMARK SUMMARY DATA FRAME
# ------------------------------------------------------------------------------
landmark_units <- c(1, 2, 3, 4, 5, 8, 10, 15, 20)

df_table_s4_data <- df_spline_plot %>%
  filter(cum_units %in% landmark_units) %>%
  mutate(
    `Cumulative Blood Exposure` = sprintf("%d %s", cum_units, if_else(cum_units == 1, "Unit", "Units")),
    `Adjusted Hazard Ratio (95% CI)*` = sprintf("%.2f (%.2f–%.2f)", hr, lower_ci, upper_ci)
  ) %>%
  select(`Cumulative Blood Exposure`, `Adjusted Hazard Ratio (95% CI)*`)

# ------------------------------------------------------------------------------
# 2. BUILD FLEXTABLE FOR PUBLICATION
# ------------------------------------------------------------------------------
ft_s4 <- flextable(df_table_s4_data) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bg(bg = "#F7F9FA", part = "header") %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2, align = "center", part = "all")

# ------------------------------------------------------------------------------
# 3. EXPORT TO WORD DOCUMENT (.DOCX) AND CSV
# ------------------------------------------------------------------------------
doc_s4 <- read_docx() %>%
  body_add_par("Supplementary Table S4. Predicted Hazard Ratios of de novo HLA Alloimmunization at Landmark Cumulative Transfusion Volumes", style = "heading 1") %>%
  body_add_par("Continuous Natural Cubic Spline Model (df = 3)", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(ft_s4) %>%
  body_add_par("") %>%
  body_add_par("*Adjusted for candidate age at listing, sex, prior organ transplantation history, and prior pregnancy history. Reference group represents unexposed candidates (0 transfused units; HR = 1.00).", style = "Normal")

# Save Word document
print(doc_s4, target = "Table_S4_Spline_Landmark_HRs.docx")

# Save flat CSV file
write.csv(df_table_s4_data, "Table_S4_Spline_Landmark_HRs.csv", row.names = FALSE)

cat("\n--- SUCCESS: Table_S4_Spline_Landmark_HRs.docx and .csv generated! ---\n")