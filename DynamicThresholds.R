# ==============================================================================
# MULTIMODAL-ROBUST ITERATIVE ALLELE CUTOFF & DIAGNOSTIC PLOTTING ENGINE
# ==============================================================================
library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(scales)

# Step 1: Load the raw histogram matrix
message("Loading MFI Histogram Table...")
df_hist <- fread("/Users/zvr/Downloads/Histogram Table.csv")
allele_names <- names(df_hist)

# CLINICAL SETUP: Global baseline fallback cushion
fallback_val <- 3000 

# ==============================================================================
# 2. RUN REGULAR DATA-DRIVEN INFLECTION HUNTING (MULTIMODAL ARCHITECTURE)
# ==============================================================================
find_column_valley_multimodal <- function(column_name, mfi_vector, global_fallback = fallback_val) {
  valid_mfis <- as.numeric(mfi_vector)
  valid_mfis <- valid_mfis[!is.na(valid_mfis) & valid_mfis > 0]
  
  if (length(valid_mfis) < 500) {
    return(data.frame(Allele = column_name, Cutoff_MFI = global_fallback, Method = "Fallback (Barren Cohort)"))
  }
  
  log_mfis <- log10(valid_mfis)
  
  # MULTIMODAL FIX: adjust = 1.4 smooths out micro-ripples within multi-tiered peaks
  dens <- density(log_mfis, n = 512, adjust = 1.4)
  
  # Search window: 500 MFI (log10=2.70) up to 15,000 MFI (log10=4.17)
  search_indices <- which(dens$x >= 2.70 & dens$x <= 4.17)
  
  if (length(search_indices) == 0) {
    return(data.frame(Allele = column_name, Cutoff_MFI = global_fallback, Method = "Fallback (Empty Window)"))
  }
  
  y_vals <- dens$y
  local_minima <- c()
  
  for (i in search_indices) {
    if (y_vals[i] < y_vals[i - 1] && y_vals[i] < y_vals[i + 1]) {
      local_minima <- c(local_minima, i)
    }
  }
  
  if (length(local_minima) > 0) {
    # MULTIMODAL FIX: Pick the EARLIEST valley floor along the X-axis (dens$x)
    # to catch the true separation boundary before low/medium/high positive tiers split
    optimal_idx <- local_minima[which.min(dens$x[local_minima])]
    cutoff_mfi  <- round(10^(dens$x[optimal_idx]))
    
    return(data.frame(
      Allele = column_name, 
      Cutoff_MFI = cutoff_mfi, 
      Method = "Success (First Visual Valley Locked)"
    ))
  } else {
    return(data.frame(Allele = column_name, Cutoff_MFI = global_fallback, Method = "Fallback (No Distinct Minimum)"))
  }
}

message("Calculating allele-specific cutoffs using X-axis tracing...")
df_allele_cutoffs <- purrr::map_df(allele_names, function(col) {
  find_column_valley_multimodal(col, df_hist[[col]])
})

# Save the updated lookup key file
fwrite(df_allele_cutoffs, "/Users/zvr/Downloads/automated_allele_cutoffs_key.csv")

# ==============================================================================
# MULTI-PAGE ANALYTE MFI HISTOGRAM ENGINE (GAP-FREE BINNING)
# ==============================================================================

library(readxl)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

message("\nRendering gap-free multi-page analyte histogram report...")

file_path <- "/Users/zvr/Downloads/Transfusion_Data.xlsx"
cut_path  <- "/Users/zvr/Downloads/Cuttoffs.xlsx"

df_mhc1_raw <- read_excel(file_path, sheet = "Luminex Confirmation MHC-I")
df_mhc2_raw <- read_excel(file_path, sheet = "Luminex Confirmation MHC-II")
df_cutoffs  <- read_excel(cut_path)

df_mhc1_long <- df_mhc1_raw %>% 
  mutate(pat_id = as.character(patient_number)) %>%
  select(pat_id, draw_date, any_of(df_cutoffs$Target)) %>%
  pivot_longer(cols = -c(pat_id, draw_date), names_to = "Target", values_to = "MFI") %>%
  mutate(global_cutoff = 3000, Class = "MHC-I")

df_mhc2_long <- df_mhc2_raw %>% 
  mutate(pat_id = as.character(patient_number)) %>%
  select(pat_id, draw_date, any_of(df_cutoffs$Target)) %>%
  pivot_longer(cols = -c(pat_id, draw_date), names_to = "Target", values_to = "MFI") %>%
  mutate(global_cutoff = 1500, Class = "MHC-II")

df_mfi_long <- bind_rows(df_mhc1_long, df_mhc2_long) %>%
  filter(!is.na(MFI) & MFI > 0) %>%
  left_join(df_cutoffs, by = "Target") %>%
  rename(custom_cutoff = Cuttoff)

df_analyte_stats <- df_mfi_long %>%
  group_by(Target, Class, global_cutoff, custom_cutoff) %>%
  summarize(.groups = "drop")

create_analyte_histogram <- function(target_name, df_data, df_stats) {
  
  df_sub <- df_data %>% filter(Target == target_name)
  stats  <- df_stats %>% filter(Target == target_name)
  
  dyn_val   <- stats$custom_cutoff
  clin_val  <- stats$global_cutoff
  class_val <- stats$Class
  
  dyn_label  <- "Dynamic Analysis Threshold"
  clin_label <- paste0("Clinical Analysis Threshold (", class_val, ")")
  
  # Adjusted bin count (bins = 22) eliminates discrete integer sampling gaps
  n_bins <- 22
  hist_counts <- hist(log10(df_sub$MFI), breaks = n_bins, plot = FALSE)$counts
  local_max   <- max(hist_counts, na.rm = TRUE)
  y_ceiling   <- local_max * 1.28
  
  p <- ggplot(df_sub, aes(x = MFI)) +
    geom_histogram(fill = "#2B3A42", color = "white", bins = n_bins, alpha = 0.85) +
    
    geom_vline(aes(xintercept = dyn_val, color = dyn_label, linetype = dyn_label), linewidth = 0.8) +
    geom_vline(aes(xintercept = clin_val, color = clin_label, linetype = clin_label), linewidth = 0.8) +
    
    annotate(
      "label", x = dyn_val, y = local_max * 1.16, 
      label = paste0(comma(dyn_val)), 
      size = 2.4, fontface = "bold", color = "#E67E22", fill = alpha("white", 0.9),
      label.size = 0.2, label.padding = unit(0.12, "lines")
    ) +
    annotate(
      "label", x = clin_val, y = local_max * 1.04, 
      label = paste0(comma(clin_val)), 
      size = 2.4, fontface = "bold", color = "#C0392B", fill = alpha("white", 0.9),
      label.size = 0.2, label.padding = unit(0.12, "lines")
    ) +
    
    scale_color_manual(name = "", values = setNames(c("#E67E22", "#C0392B"), c(dyn_label, clin_label))) +
    scale_linetype_manual(name = "", values = setNames(c("dashed", "dotted"), c(dyn_label, clin_label))) +
    
    scale_x_log10(
      labels = scales::label_comma(),
      breaks = c(100, 500, 1000, 1500, 3000, 5000, 10000, 20000)
    ) +
    
    scale_y_continuous(
      limits = c(0, y_ceiling), 
      expand = c(0, 0),
      labels = scales::label_comma()
    ) +
    
    labs(
      title = paste0("Target: ", target_name, " (", class_val, ")"),
      x = "MFI (Log10 Scale)",
      y = "Count"
    ) +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 8),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      legend.position = "top",
      legend.key.height = unit(0.3, "cm"),
      legend.text = element_text(size = 7.0, face = "bold"),
      legend.margin = margin(0, 0, -2, 0),
      plot.margin = margin(t = 6, r = 8, b = 6, l = 8)
    )
  
  return(p)
}

unique_targets <- sort(unique(df_mfi_long$Target))
plot_list <- list()

for (tgt in unique_targets) {
  plot_list[[tgt]] <- create_analyte_histogram(tgt, df_mfi_long, df_analyte_stats)
}

pdf_filename <- "Analyte_MFI_Histograms_Dynamic_vs_Clinical.pdf"

multipage_grob <- marrangeGrob(grobs = plot_list, nrow = 3, ncol = 2, top = NULL)

ggsave(
  filename = pdf_filename,
  plot = multipage_grob,
  width = 8.5,
  height = 11.0,
  units = "in",
  device = "pdf"
)

message(paste0("Success! Smooth gap-free PDF report saved as: ", pdf_filename))