# ==============================================================================
# SECONDARY ANALYSIS PIPELINE: BINARY EXPOSURE (EXPOSED VS UNEXPOSED)
# Targets: Red Blood Cells, Platelets, Plasma, Cryoprecipitate
# Core Architecture: Time-Varying Cox Models & 2x2 Publication Viewport
# ==============================================================================

library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)

# ------------------------------------------------------------------------------
# 0. CONFIGURATION & FILE INGESTION
# ------------------------------------------------------------------------------
max_x_days <- 730
file_path  <- "/Users/zvr/Downloads/Transfusion_Data.xlsx"
cut_path   <- "/Users/zvr/Downloads/Cuttoffs.xlsx"

message("Ingesting raw database registry sheets...")
df_demo_raw   <- read_excel(file_path, sheet = "Demographics")
df_trans_raw  <- read_excel(file_path, sheet = "Transfusion")
df_events_raw <- read_excel(file_path, sheet = "Patient Events")
df_mhc1_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-I")
df_mhc2_raw   <- read_excel(file_path, sheet = "Luminex Confirmation MHC-II")
df_cutoffs    <- read_excel(cut_path)

# ------------------------------------------------------------------------------
# 1. BASELINE DEMOGRAPHICS & COHORT CLEANING
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
# 2. ALLELE-SPECIFIC MFI MATRIX & ENDPOINTS
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
  summarize(first_pos_date = min(draw_date_clean, na.rm = TRUE), .groups = "drop")

df_last_lab_mfi <- bind_rows(
  df_mhc1_raw %>% mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>% select(pat_id, draw_date_clean),
  df_mhc2_raw %>% mutate(pat_id = as.character(patient_number), draw_date_clean = as.Date(draw_date)) %>% select(pat_id, draw_date_clean)
) %>%
  group_by(pat_id) %>%
  summarize(last_draw_date = max(draw_date_clean, na.rm = TRUE), .groups = "drop")

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
    end_day   = as.numeric(difftime(final_date, listing_date, units = "days"))
  ) %>%
  filter(is.na(first_pos_date) | as.numeric(difftime(first_pos_date, listing_date, units="days")) > 14) %>%
  mutate(end_day = if_else(end_day <= 0 & is.na(first_pos_date), as.numeric(difftime(coalesce(removal_date, study_freeze_date), listing_date, units="days")), end_day)) %>%
  filter(!is.na(end_day) & end_day > 0)

df_long_base <- tmerge(data1 = df_master_base, data2 = df_master_base, id = pat_id, tstart = 0, tstop = end_day)
df_long_base <- tmerge(data1 = df_long_base, data2 = df_master_base, id = pat_id, event = event(end_day, status))

# ------------------------------------------------------------------------------
# 3. HELPER & PLOTTING FUNCTIONS FOR BINARY SECONDARY ANALYSIS
# ------------------------------------------------------------------------------
format_pval <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste0("p = ", sprintf("%.3f", p)))
  }
}

group_colors_binary <- c(
  "Unexposed (0 Units)" = "#2B3A42", # Deep Slate / Charcoal
  "Exposed (≥1 Units)"  = "#C0392B"  # Deep Crimson / Red
)

extract_incidence_binary_df <- function(s_fit) {
  surv_matrix <- s_fit$surv
  d0    <- data.frame(time = s_fit$time, incidence = 1 - surv_matrix[,1], Group = "Unexposed (0 Units)")
  dexp  <- data.frame(time = s_fit$time, incidence = 1 - surv_matrix[,2], Group = "Exposed (≥1 Units)")
  
  bind_rows(dexp, d0) %>% 
    filter(time <= max_x_days) %>%
    mutate(Group = factor(Group, levels = c("Unexposed (0 Units)", "Exposed (≥1 Units)"))) %>%
    mutate(group_draw = factor(Group, levels = c("Exposed (≥1 Units)", "Unexposed (0 Units)")))
}

# Main iterative model-fitting and plotting function (Binary Exposure)
run_product_binary_analysis <- function(product_name, tag_label) {
  message(paste0("\nProcessing binary analysis for product: ", product_name, "..."))
  
  # Filter transfusion dataset strictly for target BLOOD_PRODUCT
  df_trans_prod <- df_trans_raw %>%
    filter(str_trim(BLOOD_PRODUCT) == product_name) %>%
    mutate(pat_id = as.character(`Pat ID`), trans_date = as.Date(BLOOD_START_INSTANT), units_transfused = 1) %>%
    left_join(select(df_patients_clean, pat_id, listing_date), by = "pat_id") %>%
    mutate(trans_day = as.numeric(difftime(trans_date, listing_date, units = "days"))) %>%
    filter(!is.na(trans_day) & trans_day >= 0) %>%
    arrange(pat_id, trans_day) %>%
    group_by(pat_id) %>%
    mutate(cum_units = cumsum(units_transfused)) %>%
    ungroup() %>%
    select(pat_id, trans_day, cum_units)
  
  # Build time-dependent covariate dataset with binary exposure factor
  df_long_binary <- tmerge(
    data1 = df_long_base,
    data2 = df_trans_prod,
    id = pat_id,
    cum_units = tdc(trans_day, cum_units)
  ) %>%
    mutate(
      cum_units = coalesce(cum_units, 0),
      is_exposed = if_else(cum_units >= 1, "Exposed", "Unexposed"),
      is_exposed = factor(is_exposed, levels = c("Unexposed", "Exposed"))
    )
  
  # Fit Clinically Adjusted Model
  fit_adj <- coxph(
    Surv(tstart, tstop, event) ~ is_exposed + age_at_listing + sex + 
      has_prior_transplant + has_prior_pregnancy, 
    data = df_long_binary
  )
  
  sum_adj <- summary(fit_adj)
  
  p_exp <- if("is_exposedExposed" %in% rownames(sum_adj$coefficients)) {
    sum_adj$coefficients["is_exposedExposed", "Pr(>|z|)"]
  } else {
    NA
  }
  
  box_label <- paste0("Exposed vs. Unexposed: ", format_pval(p_exp))
  
  prof_adj <- data.frame(
    is_exposed = factor(c("Unexposed", "Exposed"), levels = c("Unexposed", "Exposed")), 
    age_at_listing = 50, sex = "M", has_prior_transplant = 0, has_prior_pregnancy = 0
  )
  
  fit_c_adj <- survival::survfit(fit_adj, newdata = prof_adj)
  df_plot   <- extract_incidence_binary_df(fit_c_adj)
  
  # Generate Panel Plot
  plot_obj <- ggplot(df_plot, aes(x = time, y = incidence, color = Group, group = group_draw)) +
    geom_step(linewidth = 0.85, alpha = 0.95) +
    scale_color_manual(values = group_colors_binary) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1), 
      limits = c(0, 0.50), 
      breaks = seq(0, 0.50, by = 0.10), 
      expand = c(0.02, 0)
    ) +
    scale_x_continuous(breaks = seq(0, max_x_days, by = 200), expand = c(0.01, 0)) +
    labs(
      tag = tag_label,
      title = paste0(product_name, " (Adjusted)"), 
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
      "label", x = max_x_days * 0.40, y = 0.45, label = box_label, 
      size = 3.0, fontface = "plain", family = "sans",
      fill = alpha("white", 0.9), color = "black", label.size = 0.3, hjust = 0.5
    )
  
  return(plot_obj)
}

# ------------------------------------------------------------------------------
# 4. EXECUTE ITERATIVE ANALYSES ACROSS ALL BLOOD PRODUCTS
# ------------------------------------------------------------------------------
plot_rbc     <- run_product_binary_analysis("Red Blood Cells", "A")
plot_plt     <- run_product_binary_analysis("Platelets", "B")
plot_plasma  <- run_product_binary_analysis("Plasma", "C")
plot_cryo    <- run_product_binary_analysis("Cryoprecipitate", "D")

# Remove redundant y-axis titles for right-column panels
plot_plt  <- plot_plt + labs(y = "")
plot_cryo <- plot_cryo + labs(y = "")

# ------------------------------------------------------------------------------
# 5. COMPILE 2x2 PUBLICATION DASHBOARD
# ------------------------------------------------------------------------------
raw_legend <- cowplot::get_legend(
  plot_rbc + theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5, face = "bold"),
    legend.key.width = unit(1.2, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
)

centered_legend <- ggdraw(raw_legend)

cohort_header <- textGrob(
  "Secondary Analysis: Product-Specific Sensitization (Exposed vs. Unexposed)", 
  gp = gpar(fontface = "bold", fontsize = 13, fontfamily = "sans")
)

final_2x2_grob <- arrangeGrob(
  cohort_header,
  arrangeGrob(plot_rbc, plot_plt, plot_plasma, plot_cryo, ncol = 2, nrow = 2),
  centered_legend,
  nrow = 3, heights = c(0.6, 10, 0.8)
)

# Render in active viewport
grid.newpage()
grid.draw(final_2x2_grob)

# ------------------------------------------------------------------------------
# 6. EXPORT HIGH-RESOLUTION SUBMISSION FILES
# ------------------------------------------------------------------------------
ggsave("Figure_S1_Product_Binary_Kinetics.pdf", final_2x2_grob, width = 10.0, height = 8.5, units = "in", device = cairo_pdf)
ggsave("Figure_S1_Product_Binary_Kinetics.tiff", final_2x2_grob, width = 10.0, height = 8.5, units = "in", dpi = 300, compression = "lzw")

message("\nSuccess! Binary 2x2 secondary analysis dashboard generated and exported to PDF & 300 DPI TIFF.")