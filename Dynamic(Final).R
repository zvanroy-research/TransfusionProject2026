# ==============================================================================
# MASTER TRANSFUSION VOLUMETRIC PIPELINE: ALLELE-SPECIFIC CUTOFF MATRIX
# Core Architecture: Time-Varying Cox Proportional Hazards Regression
# Primary Model: Reconciled 3-Stage Volumetric Dashboard (0 vs 1-4 vs 5+ Units)
# Target Processing: Allele-Specific Threshold Engine (Ingests Cuttoffs.xlsx)
# Output Architecture: 2 Models (Unadjusted & Clinically Adjusted)
# ==============================================================================

library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)

# ==============================================================================
# 0. MASTER CONFIGURATION & GLOBAL PARAMETERS
# ==============================================================================
max_x_days <- 730   # 2-Year observation window cap for visualization
file_path  <- "/Users/zvr/Downloads/Transfusion_Data.xlsx"
cut_path   <- "/Users/zvr/Downloads/Cuttoffs.xlsx"

# ==============================================================================
# 1. DATABASE SHEET INGESTION
# ==============================================================================
message("Ingesting raw database registry sheets & custom cutoff keys...")
df_demo_raw   <- read_excel(file_path, sheet = "Demographics")
df_trans_raw  <- read_excel(file_path, sheet = "Transfusion")
df_events_raw <- read_excel(file_path, sheet = "Patient Events")
df_mhc1_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-I")
df_mhc2_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-II")
df_cutoffs    <- read_excel(cut_path)

# ==============================================================================
# 2. DEMOGRAPHIC CLEANING & DISCIPLINED EXCLUSIONS
# ==============================================================================
message("Applying clinical cohort exclusion parameters...")
descriptions_to_exclude <- c(
  "Kidney Living Donor", "HSC Donor", "Disease Association",
  "Altruistic Donor", "OTHER/RECIPIENT", "Liver Living Donor", "OTHER/DONOR"
)

df_patients_clean <- df_demo_raw %>%
  filter(!str_trim(description) %in% str_trim(descriptions_to_exclude)) %>%
  filter(sex %in% c("M", "F")) %>% # FILTER OUT UNKNOWN / UNCODED SEX (sex U)
  mutate(
    pat_id         = as.character(patient_number),
    listing_date   = as.Date(as.numeric(create_date), origin = "1899-12-30"),
    birth_date     = as.Date(as.numeric(dob), origin = "1899-12-30"),
    age_at_listing = as.numeric(difftime(listing_date, birth_date, units = "weeks")) / 52.25,
    sex            = factor(sex, levels = c("F", "M")) # Sets Female as reference level
  ) %>%
  select(pat_id, listing_date, age_at_listing, sex)

# ==============================================================================
# 3. LAB MATRIX ENGINE: ALLELE-SPECIFIC CRITICAL CUTOFF MATCHING
# ==============================================================================
message("Processing Luminex matrix using target probe-specific thresholds...")

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
  filter(MFI >= Cuttoff) %>% # Locks directly onto your spreadsheet's column spelling
  group_by(pat_id) %>%
  summarize(first_pos_date = min(draw_date_clean, na.rm = TRUE), .groups = "drop")

df_last_lab_mfi <- bind_rows(
  df_mhc1_raw %>% mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>% select(pat_id, draw_date_clean),
  df_mhc2_raw %>% mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>% select(pat_id, draw_date_clean)
) %>%
  group_by(pat_id) %>%
  summarize(last_draw_date = max(draw_date_clean, na.rm = TRUE), .groups = "drop")

# ==============================================================================
# 4. RETROSPECTIVE IMMUNOLOGICAL PRIMING HISTORY & COMPETITORS
# ==============================================================================
message("Reconstructing pre-listing historical events...")
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

# ==============================================================================
# 5. TIME-VARYING GRANULAR VOLUMETRIC TRANSFUSION COUNTER
# ==============================================================================
message("Compiling longitudinal transfusion record arrays...")
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

# ==============================================================================
# 6. ASSEMBLE MODEL INTEGRATION MATRIX WITH WATERFALL RECONCILIATION
# ==============================================================================
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
    
    # Strip pre-listing draws from censoring thresholds
    last_draw_date_post  = if_else(last_draw_date > listing_date, last_draw_date, as.Date(NA)),
    
    final_date = case_when(
      !is.na(first_pos_date)       ~ first_pos_date, 
      !is.na(removal_date)         ~ removal_date, 
      !is.na(last_draw_date_post)  ~ last_draw_date_post, 
      TRUE                         ~ study_freeze_date
    ),
    end_day   = as.numeric(difftime(final_date, listing_date, units = "days"))
  ) %>%
  filter(is.na(first_pos_date) | as.numeric(difftime(first_pos_date, listing_date, units="days")) > 14) %>%
  # Retain patients without post-listing draws as event-free controls
  mutate(end_day = if_else(end_day <= 0 & is.na(first_pos_date), as.numeric(difftime(coalesce(removal_date, study_freeze_date), listing_date, units="days")), end_day)) %>%
  filter(!is.na(end_day) & end_day > 0)

# ==============================================================================
# 7. PARSE TIME-DEPENDENT COVARIATE MATRICES (PRIMARY ARCHITECTURE)
# ==============================================================================
df_long_base <- tmerge(data1 = df_master_base, data2 = df_master_base, id = pat_id, tstart = 0, tstop = end_day)
df_long_base <- tmerge(data1 = df_long_base, data2 = df_master_base, id = pat_id, event = event(end_day, status))

df_long_binned <- tmerge(
  data1 = df_long_base,
  data2 = df_trans_clean,
  id = pat_id,
  cum_units = tdc(trans_day, cum_units)
) %>%
  mutate(
    cum_units = coalesce(cum_units, 0),
    vol_tier = case_when(cum_units == 0 ~ "0", cum_units >= 1 & cum_units <= 4 ~ "1-4", TRUE ~ "5+"),
    vol_tier = factor(vol_tier, levels = c("0", "1-4", "5+"))
  )

# ==============================================================================
# 8. RUN STAGE-ACCELERATING REGRESSION ARC (UNADJUSTED vs ADJUSTED)
# ==============================================================================
cat("\n======================================================================\n")
cat("       VOLUMETRIC TRANSFUSION SURVIVAL MODELS (ALLELE CUTOFFS)        \n")
cat("======================================================================\n")

message("\n--- MODEL 1: UNADJUSTED MODEL ---")
fit_stage1d <- coxph(Surv(tstart, tstop, event) ~ vol_tier, data = df_long_binned)
print(summary(fit_stage1d))

message("\n--- MODEL 2: CLINICALLY ADJUSTED MODEL (DEMOGRAPHICS, PREGNANCY, TRANSPLANT) ---")
fit_stage2d <- coxph(
  Surv(tstart, tstop, event) ~ vol_tier + age_at_listing + sex + 
    has_prior_transplant + has_prior_pregnancy, 
  data = df_long_binned
)
print(summary(fit_stage2d))

# ==============================================================================
# 9. DUAL-PANEL PUBLICATION DASHBOARD (UNEXPOSED LINE FORCED ON TOP)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)

message("\nRendering final publication dashboard and exporting high-res file...")

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS & DATA PREPARATION
# ------------------------------------------------------------------------------
extract_incidence_df <- function(s_fit) {
  surv_matrix <- s_fit$surv
  d0    <- data.frame(time = s_fit$time, incidence = 1 - surv_matrix[,1], Group = "Unexposed (0 Units)")
  dlow  <- data.frame(time = s_fit$time, incidence = 1 - surv_matrix[,2], Group = "Low Exposure (1-4 Units)")
  dhigh <- data.frame(time = s_fit$time, incidence = 1 - surv_matrix[,3], Group = "High Exposure (≥ 5 Units)")
  
  bind_rows(dlow, dhigh, d0) %>% 
    filter(time <= max_x_days) %>%
    # Group: controls legend order (Unexposed first)
    mutate(Group = factor(Group, levels = c("Unexposed (0 Units)", "Low Exposure (1-4 Units)", "High Exposure (≥ 5 Units)"))) %>%
    # group_draw: controls rendering layer order (Unexposed last = drawn on top)
    mutate(group_draw = factor(Group, levels = c("Low Exposure (1-4 Units)", "High Exposure (≥ 5 Units)", "Unexposed (0 Units)")))
}

format_pval <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste0("p = ", sprintf("%.3f", p)))
  }
}

group_colors <- c(
  "Unexposed (0 Units)"     = "#2B3A42", # Deep Slate / Charcoal
  "Low Exposure (1-4 Units)" = "#E67E22", # Warm Amber / Ochre
  "High Exposure (≥ 5 Units)" = "#C0392B"  # Deep Crimson / Red
)

# ------------------------------------------------------------------------------
# 2. PANEL A: UNADJUSTED MODEL
# ------------------------------------------------------------------------------
sum_s1    <- summary(fit_stage1d)
p_s1_low  <- sum_s1$coefficients["vol_tier1-4", "Pr(>|z|)"]
p_s1_high <- sum_s1$coefficients["vol_tier5+", "Pr(>|z|)"]

box_s1    <- paste0(
  "Low Exposure vs. Unexposed: ", format_pval(p_s1_low), "\n",
  "High Exposure vs. Unexposed: ", format_pval(p_s1_high)
)

fit_c_s1   <- survival::survfit(fit_stage1d, newdata = data.frame(vol_tier = factor(c("0", "1-4", "5+"), levels = c("0", "1-4", "5+"))))
df_plot_s1 <- extract_incidence_df(fit_c_s1)

plot_s1 <- ggplot(df_plot_s1, aes(x = time, y = incidence, color = Group, group = group_draw)) +
  geom_step(linewidth = 0.85, alpha = 0.95) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), 
    limits = c(0, 0.50), 
    breaks = seq(0, 0.50, by = 0.10), 
    expand = c(0.02, 0)
  ) +
  scale_x_continuous(breaks = seq(0, max_x_days, by = 200), expand = c(0.01, 0)) +
  labs(
    tag = "A",
    title = "Unadjusted Model", 
    x = "Days on Waitlist", 
    y = "Probability of Sensitization"
  ) +
  theme_classic(base_size = 11, base_family = "sans") + 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10)
  ) +
  annotate(
    "label", x = max_x_days * 0.40, y = 0.45, label = box_s1, 
    size = 3.0, fontface = "plain", family = "sans",
    fill = alpha("white", 0.9), color = "black", label.size = 0.3, hjust = 0.5
  )

# ------------------------------------------------------------------------------
# 3. PANEL B: ADJUSTED MODEL
# ------------------------------------------------------------------------------
sum_s2    <- summary(fit_stage2d)
p_s2_low  <- sum_s2$coefficients["vol_tier1-4", "Pr(>|z|)"]
p_s2_high <- sum_s2$coefficients["vol_tier5+", "Pr(>|z|)"]

box_s2    <- paste0(
  "Low Exposure vs. Unexposed: ", format_pval(p_s2_low), "\n",
  "High Exposure vs. Unexposed: ", format_pval(p_s2_high)
)

prof_s2  <- data.frame(vol_tier = factor(c("0", "1-4", "5+"), levels = c("0", "1-4", "5+")), age_at_listing = 50, sex = "M", has_prior_transplant = 0, has_prior_pregnancy = 0)
fit_c_s2 <- survival::survfit(fit_stage2d, newdata = prof_s2)
df_plot_s2 <- extract_incidence_df(fit_c_s2)

plot_s2 <- ggplot(df_plot_s2, aes(x = time, y = incidence, color = Group, group = group_draw)) +
  geom_step(linewidth = 0.85, alpha = 0.95) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), 
    limits = c(0, 0.50), 
    breaks = seq(0, 0.50, by = 0.10),
    expand = c(0.02, 0)
  ) +
  scale_x_continuous(breaks = seq(0, max_x_days, by = 200), expand = c(0.01, 0)) +
  labs(
    tag = "B",
    title = "Adjusted Model", 
    x = "Days on Waitlist", 
    y = ""
  ) +
  theme_classic(base_size = 11, base_family = "sans") + 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10)
  ) +
  annotate(
    "label", x = max_x_days * 0.40, y = 0.45, label = box_s2, 
    size = 3.0, fontface = "plain", family = "sans",
    fill = alpha("white", 0.9), color = "black", label.size = 0.3, hjust = 0.5
  )

# ------------------------------------------------------------------------------
# 4. COMPILE VIEWPORT & SAVE PUBLICATION FILES
# ------------------------------------------------------------------------------
plot_for_legend <- ggplot(df_plot_s1, aes(x = time, y = incidence, color = Group)) +
  geom_step() +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5, face = "bold"),
    legend.key.width = unit(1.2, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

raw_legend <- cowplot::get_legend(plot_for_legend)
centered_legend <- ggdraw(raw_legend)

cohort_header <- textGrob(
  "Dynamic Analysis", 
  gp = gpar(fontface = "bold", fontsize = 12, fontfamily = "sans")
)

final_grob <- arrangeGrob(
  cohort_header,
  arrangeGrob(plot_s1, plot_s2, ncol = 2),
  centered_legend,
  nrow = 3, heights = c(0.8, 10, 1)
)

# Display in active session
grid.newpage()
grid.draw(final_grob)

# Save high-resolution PDF & TIFF files for journal submission
ggsave("Figure_2_Dynamic.pdf", final_grob, width = 8.5, height = 5.2, units = "in", device = cairo_pdf)
ggsave("Figure_2_Dynamic.tiff", final_grob, width = 8.5, height = 5.2, units = "in", dpi = 300, compression = "lzw")

message("Success! Unexposed line now rendered on top. Figures exported.")