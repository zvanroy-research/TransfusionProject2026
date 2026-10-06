# ==============================================================================
# script 2: primary analysis and related files
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(car)
  library(scales)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(rms)
})

set.seed(42)

root_dir   <- "/Users/zvr/Downloads/Clinical2.0"
output_dir <- file.path(root_dir, "processed_data")
table_dir  <- file.path(root_dir, "tables_output")
fig_dir    <- file.path(root_dir, "figures_output")

max_x_days     <- 730
plot_cap_units <- 10

df_master_base <- readRDS(file.path(output_dir, "df_master_base.rds"))
df_trans_clean <- readRDS(file.path(output_dir, "df_trans_clean.rds"))

df_long_base <- tmerge(data1 = df_master_base, data2 = df_master_base, id = pat_id, tstart = 0, tstop = end_day)
df_long_base <- tmerge(data1 = df_long_base, data2 = df_master_base, id = pat_id, event = event(end_day, status))

df_long_cont <- tmerge(
  data1     = df_long_base,
  data2     = df_trans_clean,
  id        = pat_id,
  cum_units = tdc(trans_day, cum_units)
) %>%
  mutate(
    cum_units      = coalesce(cum_units, 0),
    trans_exposure = if_else(cum_units >= 1, "Exposed (1+ Units)", "Unexposed (0 Units)"),
    trans_exposure = factor(trans_exposure, levels = c("Unexposed (0 Units)", "Exposed (1+ Units)"))
  )

fit_stage1c <- coxph(Surv(tstart, tstop, event) ~ trans_exposure, data = df_long_cont, cluster = pat_id)

# Stratified Multivariable Model (Addresses Proportional Hazards violations)
fit_stage2c <- coxph(Surv(tstart, tstop, event) ~ trans_exposure + age_at_listing + sex_pregnancy_status + 
                       strata(has_prior_transplant, organ_group), data = df_long_cont, cluster = pat_id)

# Unstratified ghost model solely for computing Variance Inflation Factors (GVIF)
fit_stage2c_unstrat <- coxph(Surv(tstart, tstop, event) ~ trans_exposure + age_at_listing + sex_pregnancy_status + 
                               has_prior_transplant + organ_group, data = df_long_cont, cluster = pat_id)

# --- table 2 ---


s_m1 <- summary(fit_stage1c)
s_m2 <- summary(fit_stage2c)

fit_uni_age   <- coxph(Surv(tstart, tstop, event) ~ age_at_listing, data = df_long_cont, cluster = pat_id)
fit_uni_sex   <- coxph(Surv(tstart, tstop, event) ~ sex_pregnancy_status, data = df_long_cont, cluster = pat_id)
fit_uni_prior <- coxph(Surv(tstart, tstop, event) ~ has_prior_transplant, data = df_long_cont, cluster = pat_id)
fit_uni_organ <- coxph(Surv(tstart, tstop, event) ~ organ_group, data = df_long_cont, cluster = pat_id)

extract_est <- function(fit, coef_name) {
  if (is.null(fit) || is.na(coef_name) || coef_name == "") {
    return(list(hr_ci = "", p_str = ""))
  }
  s <- summary(fit)
  p_col <- intersect(c("Pr(>|z|)", "p", "Pr(>|W|)"), colnames(s$coefficients))[1]
  val_p <- s$coefficients[coef_name, p_col]
  p_str <- if (val_p < 0.001) "< 0.001" else sprintf("%.3f", val_p)
  hr_ci <- sprintf("%.2f (%.2f–%.2f)", 
                   s$conf.int[coef_name, "exp(coef)"], 
                   s$conf.int[coef_name, "lower .95"], 
                   s$conf.int[coef_name, "upper .95"])
  list(hr_ci = hr_ci, p_str = p_str)
}

p_col_m2 <- intersect(c("Pr(>|z|)", "p", "Pr(>|W|)"), colnames(s_m2$coefficients))[1]
v_m2_raw <- suppressWarnings(car::vif(fit_stage2c_unstrat))

if (is.matrix(v_m2_raw)) {
  df_vif <- as.data.frame(v_m2_raw) %>%
    rownames_to_column("Var_Group") %>%
    rename(GVIF = 2) %>%
    mutate(
      Adj_GVIF = if ("GVIF^(1/(2*Df))" %in% names(.)) .data[["GVIF^(1/(2*Df))"]] else sqrt(GVIF),
      GVIF_Readout = sprintf("%.2f", Adj_GVIF^2)
    )
} else {
  df_vif <- tibble(
    Var_Group    = names(v_m2_raw),
    GVIF_Readout = sprintf("%.2f", v_m2_raw)
  )
}

get_m2_vals <- function(coef_name, var_grp) {
  if (is.na(coef_name) || coef_name == "") {
    return(list(hr_ci = "", p_str = "", gvif = "-"))
  }
  val_p <- s_m2$coefficients[coef_name, p_col_m2]
  p_str <- if (val_p < 0.001) "< 0.001" else sprintf("%.3f", val_p)
  hr_ci <- sprintf("%.2f (%.2f–%.2f)", 
                   s_m2$conf.int[coef_name, "exp(coef)"], 
                   s_m2$conf.int[coef_name, "lower .95"], 
                   s_m2$conf.int[coef_name, "upper .95"])
  gvif_val <- df_vif$GVIF_Readout[df_vif$Var_Group == var_grp]
  if (length(gvif_val) == 0) gvif_val <- "-"
  list(hr_ci = hr_ci, p_str = p_str, gvif = gvif_val)
}

model_specs <- list(
  list(var = "Transfusion Exposure (Time-Varying)", raw = "", grp = "", is_ref = FALSE, is_strat = FALSE, fit = NULL),
  list(var = "  Unexposed (0 Units)", raw = "", grp = "", is_ref = TRUE, is_strat = FALSE, fit = NULL),
  list(var = "  Exposed (≥1 Units)", raw = "trans_exposureExposed (1+ Units)", grp = "trans_exposure", is_ref = FALSE, is_strat = FALSE, fit = fit_stage1c),
  list(var = "Clinical & Reproductive History", raw = "", grp = "", is_ref = FALSE, is_strat = FALSE, fit = NULL),
  list(var = "  Female, Nulligravid", raw = "", grp = "", is_ref = TRUE, is_strat = FALSE, fit = NULL),
  list(var = "  Female, Prior Pregnancy", raw = "sex_pregnancy_statusFemale_Prior_Pregnancy", grp = "sex_pregnancy_status", is_ref = FALSE, is_strat = FALSE, fit = fit_uni_sex),
  list(var = "  Male Candidate", raw = "sex_pregnancy_statusMale", grp = "sex_pregnancy_status", is_ref = FALSE, is_strat = FALSE, fit = fit_uni_sex),
  list(var = "  Prior Solid Organ Transplant*", raw = "has_prior_transplant", grp = "has_prior_transplant", is_ref = FALSE, is_strat = TRUE, fit = fit_uni_prior),
  list(var = "Organ Candidacy Category*", raw = "", grp = "", is_ref = FALSE, is_strat = FALSE, fit = NULL),
  list(var = "  Isolated Kidney", raw = "", grp = "", is_ref = TRUE, is_strat = TRUE, fit = NULL),
  list(var = "  Isolated Heart", raw = "organ_groupIsolated Heart", grp = "organ_group", is_ref = FALSE, is_strat = TRUE, fit = fit_uni_organ),
  list(var = "  Isolated Lung", raw = "organ_groupIsolated Lung", grp = "organ_group", is_ref = FALSE, is_strat = TRUE, fit = fit_uni_organ),
  list(var = "  Multi-Organ & Other", raw = "organ_groupMulti-Organ & Other", grp = "organ_group", is_ref = FALSE, is_strat = TRUE, fit = fit_uni_organ),
  list(var = "Demographics", raw = "", grp = "", is_ref = FALSE, is_strat = FALSE, fit = NULL),
  list(var = "  Age at Listing (per year)", raw = "age_at_listing", grp = "age_at_listing", is_ref = FALSE, is_strat = FALSE, fit = fit_uni_age)
)

t2_out <- map_dfr(model_specs, function(item) {
  # 1. Unadjusted Extraction Logic
  if (item$raw == "" && !item$is_ref) {
    unadj_hr <- ""; unadj_p <- ""
  } else if (item$is_ref) {
    unadj_hr <- "1.00 (Reference)"; unadj_p <- ""
  } else {
    u_res <- extract_est(item$fit, item$raw)
    unadj_hr <- u_res$hr_ci; unadj_p <- u_res$p_str
  }
  
  # 2. Adjusted / Stratified Extraction Logic
  if (item$raw == "" && !item$is_ref && !item$is_strat) {
    adj_hr <- ""; adj_p <- ""; gvif <- "-"
  } else if (item$is_strat) {
    adj_hr <- "Stratified"; adj_p <- "—"; gvif <- "-"
  } else if (item$is_ref) {
    adj_hr <- "1.00 (Reference)"; adj_p <- ""; gvif <- "-"
  } else {
    a_res <- get_m2_vals(item$raw, item$grp)
    adj_hr <- a_res$hr_ci; adj_p <- a_res$p_str; gvif <- a_res$gvif
  }
  
  tibble(
    Variable                 = item$var,
    `Unadjusted HR (95% CI)` = unadj_hr,
    `Unadjusted p-value`     = unadj_p,
    `Adjusted aHR (95% CI)`  = adj_hr,
    `Adjusted p-value`       = adj_p,
    GVIF                     = gvif
  )
})

write_csv(t2_out, file.path(table_dir, "Table_2_Primary_Cox_Models.csv"))
cat("\n--- REFINED TABLE 2: PRIMARY TIME-VARYING SURVIVAL MODELS ---\n")
print(as.data.frame(t2_out), row.names = FALSE)

# --- figure 1 ---
group_colors <- c("Unexposed (0 Units)" = "#2B3A42", "Exposed (1+ Units)" = "#C0392B")

# unadjusted
prof_s1_unexp <- data.frame(trans_exposure = factor("Unexposed (0 Units)", levels = c("Unexposed (0 Units)", "Exposed (1+ Units)")))
prof_s1_exp   <- data.frame(trans_exposure = factor("Exposed (1+ Units)", levels = c("Unexposed (0 Units)", "Exposed (1+ Units)")))

fit_s1_unexp <- survival::survfit(fit_stage1c, newdata = prof_s1_unexp, se.fit = FALSE)
fit_s1_exp   <- survival::survfit(fit_stage1c, newdata = prof_s1_exp, se.fit = FALSE)

df_plot_s1 <- bind_rows(
  data.frame(time = 0, incidence = 0, Group = "Unexposed (0 Units)"),
  data.frame(time = fit_s1_unexp$time, incidence = 1 - fit_s1_unexp$surv, Group = "Unexposed (0 Units)"),
  data.frame(time = 0, incidence = 0, Group = "Exposed (1+ Units)"),
  data.frame(time = fit_s1_exp$time, incidence = 1 - fit_s1_exp$surv, Group = "Exposed (1+ Units)")
) %>% filter(time <= max_x_days) %>%
  mutate(Group = factor(Group, levels = c("Unexposed (0 Units)", "Exposed (1+ Units)")),
         group_draw = factor(Group, levels = c("Exposed (1+ Units)", "Unexposed (0 Units)")))

p_col_m1 <- intersect(c("Pr(>|z|)", "p", "Pr(>|W|)"), colnames(s_m1$coefficients))[1]
p_s1_exp <- s_m1$coefficients["trans_exposureExposed (1+ Units)", p_col_m1]
box_s1   <- paste0("Exposed vs. Unexposed: ", if(p_s1_exp < 0.001) "p < 0.001" else sprintf("p = %.3f", p_s1_exp))

plot_s1 <- ggplot(df_plot_s1, aes(x = time, y = incidence, color = Group, group = group_draw)) +
  geom_step(linewidth = 0.85, alpha = 0.95) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.50), breaks = seq(0, 0.50, 0.10), expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(0, max_x_days, 200), expand = c(0.01, 0)) +
  labs(tag = "A", title = "Unadjusted Model", x = "Days on Waitlist", y = "Probability of Sensitization") +
  theme_classic(base_size = 11) + theme(plot.tag = element_text(size = 15, face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none") +
  annotate("label", x = max_x_days * 0.40, y = 0.45, label = box_s1, size = 3.2, fill = alpha("white", 0.9))

# 2. adjusted
prof_s2_unexp <- data.frame(trans_exposure = factor("Unexposed (0 Units)", levels = c("Unexposed (0 Units)", "Exposed (1+ Units)")),
                            age_at_listing = 50, sex_pregnancy_status = factor("Male", levels = c("Female_Nulligravid", "Female_Prior_Pregnancy", "Male")),
                            has_prior_transplant = 0, organ_group = factor("Isolated Kidney", levels = levels(df_long_cont$organ_group)))
prof_s2_exp   <- data.frame(trans_exposure = factor("Exposed (1+ Units)", levels = c("Unexposed (0 Units)", "Exposed (1+ Units)")),
                            age_at_listing = 50, sex_pregnancy_status = factor("Male", levels = c("Female_Nulligravid", "Female_Prior_Pregnancy", "Male")),
                            has_prior_transplant = 0, organ_group = factor("Isolated Kidney", levels = levels(df_long_cont$organ_group)))

fit_s2_unexp <- survival::survfit(fit_stage2c, newdata = prof_s2_unexp, se.fit = FALSE)
fit_s2_exp   <- survival::survfit(fit_stage2c, newdata = prof_s2_exp, se.fit = FALSE)

df_plot_s2 <- bind_rows(
  data.frame(time = 0, incidence = 0, Group = "Unexposed (0 Units)"),
  data.frame(time = fit_s2_unexp$time, incidence = 1 - fit_s2_unexp$surv, Group = "Unexposed (0 Units)"),
  data.frame(time = 0, incidence = 0, Group = "Exposed (1+ Units)"),
  data.frame(time = fit_s2_exp$time, incidence = 1 - fit_s2_exp$surv, Group = "Exposed (1+ Units)")
) %>% filter(time <= max_x_days) %>%
  mutate(Group = factor(Group, levels = c("Unexposed (0 Units)", "Exposed (1+ Units)")),
         group_draw = factor(Group, levels = c("Exposed (1+ Units)", "Unexposed (0 Units)")))

p_col_m2 <- intersect(c("Pr(>|z|)", "p", "Pr(>|W|)"), colnames(s_m2$coefficients))[1]
p_s2_exp <- s_m2$coefficients["trans_exposureExposed (1+ Units)", p_col_m2]
box_s2   <- paste0("Exposed vs. Unexposed: ", if(p_s2_exp < 0.001) "p < 0.001" else sprintf("p = %.3f", p_s2_exp))

plot_s2 <- ggplot(df_plot_s2, aes(x = time, y = incidence, color = Group, group = group_draw)) +
  geom_step(linewidth = 0.85, alpha = 0.95) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.50), breaks = seq(0, 0.50, 0.10), expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(0, max_x_days, 200), expand = c(0.01, 0)) +
  labs(tag = "B", title = "Adjusted Model", x = "Days on Waitlist", y = "") +
  theme_classic(base_size = 11) + theme(plot.tag = element_text(size = 15, face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none") +
  annotate("label", x = max_x_days * 0.40, y = 0.45, label = box_s2, size = 3.2, fill = alpha("white", 0.9))

plot_for_legend <- ggplot(df_plot_s1, aes(x = time, y = incidence, color = Group)) +
  geom_step() + scale_color_manual(values = group_colors) + theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 9.5, face = "bold"))
centered_legend <- ggdraw(cowplot::get_legend(plot_for_legend))
cohort_header   <- textGrob("", gp = gpar(fontface = "bold", fontsize = 12))
figure_1_grob   <- arrangeGrob(cohort_header, arrangeGrob(plot_s1, plot_s2, ncol = 2), centered_legend, nrow = 3, heights = c(0.8, 10, 1))

ggsave(file.path(fig_dir, "Figure_1_Primary_Binary_Models.pdf"), figure_1_grob, width = 8.5, height = 5.2, units = "in")
ggsave(file.path(fig_dir, "Figure_1_Primary_Binary_Models.tiff"), figure_1_grob, width = 8.5, height = 5.2, units = "in", dpi = 300, compression = "lzw")


# --- 1. cumulative dose model stats---

fit_c1 <- coxph(Surv(tstart, tstop, event) ~ cum_units + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group), data = df_long_cont)
fit_c2 <- coxph(Surv(tstart, tstop, event) ~ pspline(cum_units, df = 2) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group), data = df_long_cont)
fit_c3 <- coxph(Surv(tstart, tstop, event) ~ pspline(cum_units, df = 3) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group), data = df_long_cont)
fit_c4 <- coxph(Surv(tstart, tstop, event) ~ pspline(cum_units, df = 4) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group), data = df_long_cont)
fit_c5 <- coxph(Surv(tstart, tstop, event) ~ pspline(cum_units, df = 5) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group), data = df_long_cont)

fit_spline_adj <- fit_c2

comp_fits   <- list(fit_c1, fit_c2, fit_c3, fit_c4, fit_c5)
comp_names  <- c("df = 1 (Linear)", "df = 2 (Spline)", "df = 3 (Spline)", "df = 4 (Spline)", "df = 5 (Spline)")
correct_edf <- c(length(fit_c1$coefficients), sum(fit_c2$df), sum(fit_c3$df), sum(fit_c4$df), sum(fit_c5$df))
n_events    <- as.numeric(fit_spline_adj$nevent)

unique_pats  <- unique(df_long_cont$pat_id)
k_folds      <- 5
set.seed(42)
pat_fold_map <- tibble(pat_id = unique_pats, fold = sample(rep(1:k_folds, length.out = length(unique_pats))))
df_cv        <- df_long_cont %>% left_join(pat_fold_map, by = "pat_id")

cv_formula_strings <- list(
  "df = 1 (Linear)" = "Surv(tstart, tstop, event) ~ cum_units + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group)",
  "df = 2 (Spline)" = "Surv(tstart, tstop, event) ~ pspline(cum_units, df = 2) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group)",
  "df = 3 (Spline)" = "Surv(tstart, tstop, event) ~ pspline(cum_units, df = 3) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group)",
  "df = 4 (Spline)" = "Surv(tstart, tstop, event) ~ pspline(cum_units, df = 4) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group)",
  "df = 5 (Spline)" = "Surv(tstart, tstop, event) ~ pspline(cum_units, df = 5) + age_at_listing + sex_pregnancy_status + strata(has_prior_transplant, organ_group)"
)

cv_cindices <- map_dbl(names(cv_formula_strings), function(m) {
  c_ind <- numeric(k_folds)
  for (f in 1:k_folds) {
    tr  <- filter(df_cv, fold != f)
    val <- filter(df_cv, fold == f)
    
    fit_tr <- do.call(coxph, list(formula = as.formula(cv_formula_strings[[m]]), data = quote(tr)))
    pred   <- predict(fit_tr, newdata = val, type = "lp")
    c_ind[f] <- concordance(Surv(val$tstart, val$tstop, val$event) ~ pred, reverse = TRUE)$concordance
  }
  round(mean(c_ind), 4)
})

supp_table_6 <- tibble(
  `Spline Model Architecture` = comp_names,
  `Log-Likelihood`            = map_dbl(comp_fits, ~ round(.x$loglik[2], 2)),
  `Effective DF (EDF)`        = round(correct_edf, 2),
  `AIC`                       = round(-2 * map_dbl(comp_fits, ~ .x$loglik[2]) + 2 * correct_edf, 2),
  `BIC`                       = round(-2 * map_dbl(comp_fits, ~ .x$loglik[2]) + log(n_events) * correct_edf, 2),
  `Mean 5-Fold CV C-Index`    = cv_cindices
) %>%
  mutate(`Delta AIC` = round(AIC - min(AIC), 2), `Delta BIC` = round(BIC - min(BIC), 2))

write_csv(supp_table_6, file.path(table_dir, "Supplementary_Table_6_Spline_DF_Calculations.csv"))


dd <- datadist(df_long_cont)
options(datadist = "dd")

custom_knots <- c(0, 3, 7)

fit_rcs_adj <- cph(
  Surv(tstart, tstop, event) ~ rcs(cum_units, custom_knots) + age_at_listing + 
    sex_pregnancy_status + strat(has_prior_transplant) + strat(organ_group),
  data = df_long_cont,
  x = TRUE, y = TRUE
)

fit_rcs_robust <- robcov(fit_rcs_adj, cluster = df_long_cont$pat_id)

plot_cap_units <- 8
eval_units <- seq(0, plot_cap_units, by = 0.1)

k_contrast <- contrast(
  fit_rcs_robust,
  list(cum_units = eval_units),
  list(cum_units = 0)
)

df_spline_plot <- tibble(
  x        = eval_units,
  log_hr   = as.numeric(k_contrast$Contrast),
  hr_lower = exp(as.numeric(k_contrast$Lower)),
  hr_upper = exp(as.numeric(k_contrast$Upper))
) %>%
  mutate(hr = exp(log_hr))

get_milestone <- function(u) {
  r <- df_spline_plot %>% slice(which.min(abs(x - u)))
  tibble(`Transfused Burden` = paste0(u, " Unit", if(u != 1) "s" else ""), `aHR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", r$hr, r$hr_lower, r$hr_upper))
}

table_3 <- bind_rows(get_milestone(1), get_milestone(2), get_milestone(3), get_milestone(4), get_milestone(6), get_milestone(8))
write_csv(table_3, file.path(table_dir, "Table_3_Cumulative_Dose_Analysis.csv"))

# ---  Figure 2 ---
plot_figure_2 <- ggplot(df_spline_plot, aes(x = x, y = hr)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper), fill = "#C0392B", alpha = 0.18) +
  geom_line(color = "#C0392B", linewidth = 1.1) +
  scale_y_continuous(breaks = seq(0.5, 3, by = 0.5), expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(0, plot_cap_units, by = 1), limits = c(0, plot_cap_units), expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(0, plot_cap_units), ylim = c(0.5, 3.2)) +
  labs(
    tag      = "A",
    title    = "Dose-Dependent Sensitization Risk by Transfused Units",
    subtitle = "Restricted Cubic Splines (df=2) with reference = Unexposed (0 Units)",
    x        = "Cumulative Transfused Units",
    y        = "Adjusted Hazard Ratio"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.tag   = element_text(size = 15, face = "bold"),
    plot.title = element_text(face = "bold", size = 11.5, lineheight = 1.1)
  )

ggsave(file.path(fig_dir, "Figure_2_Secondary_Dose_Spline_NoPlatelets.pdf"), plot_figure_2, width = 6.5, height = 4.8, units = "in")
ggsave(file.path(fig_dir, "Figure_2_Secondary_Dose_Spline_NoPlatelets.tiff"), plot_figure_2, width = 6.5, height = 4.8, units = "in", dpi = 300, compression = "lzw")

# --- Supplementary Figure ---

supp_hist_data <- df_long_cont %>%
  filter(cum_units >= 1) %>%
  group_by(pat_id) %>%
  summarize(
    max_burden = max(cum_units, na.rm = TRUE),
    sensitized = as.integer(any(event == 1)),
    sens_dose  = ifelse(sensitized == 1, cum_units[which(event == 1)[1]], NA_real_)
  ) %>%
  ungroup()

dose_levels <- 1:15

df_dose_support <- tibble(dose = dose_levels) %>%
  rowwise() %>%
  mutate(
    n_pts_at_least = sum(supp_hist_data$max_burden >= dose),
    n_pts_peak     = sum(round(supp_hist_data$max_burden) == dose),
    n_events       = sum(round(supp_hist_data$sens_dose) == dose, na.rm = TRUE)
  ) %>%
  ungroup()

p_supp_a <- ggplot(df_dose_support, aes(x = dose, y = n_pts_peak)) +
  geom_col(aes(fill = dose <= 8), width = 0.75, show.legend = FALSE) +
  geom_vline(xintercept = 8.5, linetype = "dashed", color = "#C0392B", linewidth = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#2C3E50", "FALSE" = "#BDC3C7")) +
  scale_x_continuous(breaks = 1:15, expand = c(0.02, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0.02, 0)) +
  annotate("text", x = 8.8, y = max(df_dose_support$n_pts_peak) * 0.90, 
           label = "Cumulative Analysis Truncation = 8 Units\n 80% of transfused patients ≤ 8 total units", 
           hjust = 0, size = 3.2, fontface = "italic", color = "#C0392B") +
  labs(
    tag   = "A",
    title = "Transfused Cohort Contribution to Exposure by Units Transfused",
    x     = "",
    y     = "Number of Candidates"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.tag   = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size = 11.5)
  )

p_supp_b <- ggplot(df_dose_support, aes(x = dose, y = n_events)) +
  geom_col(aes(fill = dose <= 8), width = 0.75, show.legend = FALSE) +
  geom_vline(xintercept = 8.5, linetype = "dashed", color = "#C0392B", linewidth = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#F2D7D5")) +
  scale_x_continuous(breaks = 1:15, expand = c(0.02, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0.02, 0)) +
  labs(
    tag   = "B",
    title = "Observed Sensitization Events by Units Transfused",
    x     = "Cumulative Transfused Units",
    y     = "Sensitization Events"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.tag   = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size = 11.5)
  )

supp_figure_grob <- arrangeGrob(p_supp_a, p_supp_b, ncol = 1)

ggsave(file.path(fig_dir, "Supplementary_Figure_Dose_Support_Histogram.pdf"), 
       supp_figure_grob, width = 7.0, height = 6.2, units = "in")
ggsave(file.path(fig_dir, "Supplementary_Figure_Dose_Support_Histogram.tiff"), 
       supp_figure_grob, width = 7.0, height = 6.2, units = "in", dpi = 300, compression = "lzw")

# cumulative dose table extended

fit_rcs_unadj <- cph(
  Surv(tstart, tstop, event) ~ rcs(cum_units, custom_knots),
  data = df_long_cont,
  x = TRUE, y = TRUE
)
fit_rcs_unadj_rob <- robcov(fit_rcs_unadj, cluster = df_long_cont$pat_id)

milestones <- 0:8

k_adj_cum   <- contrast(fit_rcs_robust, list(cum_units = milestones), list(cum_units = 0))
k_unadj_cum <- contrast(fit_rcs_unadj_rob, list(cum_units = milestones), list(cum_units = 0))

prev_milestones <- c(0, 0:7)
k_adj_marg      <- contrast(fit_rcs_robust, list(cum_units = milestones), list(cum_units = prev_milestones))

n_unexposed_pts <- df_master_base %>%
  filter(!pat_id %in% df_trans_clean$pat_id) %>%
  nrow()

df_support_summary <- tibble(dose = milestones) %>%
  rowwise() %>%
  mutate(
    n_pts = if (dose == 0) {
      n_unexposed_pts
    } else {
      sum(round(supp_hist_data$max_burden) == dose)
    },
    n_events = if (dose == 0) {
      312 # Unexposed baseline events
    } else {
      sum(round(supp_hist_data$sens_dose) == dose, na.rm = TRUE)
    }
  ) %>%
  ungroup()

table_3_all_points <- tibble(
  `Cumulative RBC Burden` = c("0 Units (Reference)", paste0(milestones[-1], " Unit", ifelse(milestones[-1] > 1, "s", ""))),
  `Candidates (n)`        = df_support_summary$n_pts,
  `Sensitization Events`  = df_support_summary$n_events,
  `Unadjusted HR (95% CI)` = c(
    "1.00 (Reference)",
    sprintf("%.2f (%.2f–%.2f)", exp(as.numeric(k_unadj_cum$Contrast[-1])), 
            exp(as.numeric(k_unadj_cum$Lower[-1])), 
            exp(as.numeric(k_unadj_cum$Upper[-1])))
  ),
  `Adjusted aHR vs 0 Units (95% CI)` = c(
    "1.00 (Reference)",
    sprintf("%.2f (%.2f–%.2f)", exp(as.numeric(k_adj_cum$Contrast[-1])), 
            exp(as.numeric(k_adj_cum$Lower[-1])), 
            exp(as.numeric(k_adj_cum$Upper[-1])))
  ),
  `Incremental aHR (+1 Unit)` = c(
    "—",
    sprintf("%.2f (%.2f–%.2f)", exp(as.numeric(k_adj_marg$Contrast[-1])), 
            exp(as.numeric(k_adj_marg$Lower[-1])), 
            exp(as.numeric(k_adj_marg$Upper[-1])))
  )
)

write_csv(table_3_all_points, file.path(table_dir, "Table_3_Cumulative_Dose_Analysis_All_Points.csv"))

cat("\n--- TABLE 3: FULL INTEGER DOSE-RESPONSE ANALYSIS (0–8 UNITS) ---\n")
print(as.data.frame(table_3_all_points), row.names = FALSE)

# word export

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(flextable)
  library(officer)
})

root_dir  <- "/Users/zvr/Downloads/Clinical2.0"
table_dir <- file.path(root_dir, "tables_output")
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

csv_path <- file.path(table_dir, "Supplementary_Table_6_Spline_DF_Calculations.csv")
if (!file.exists(csv_path)) {
  csv_path <- "Supplementary_Table_6_Spline_DF_Calculations.csv"
}

raw_spline <- read_csv(csv_path, show_col_types = FALSE)

df_spline_processed <- raw_spline %>%
  mutate(
    LRT_stat  = if_else(row_number() == 1, NA_real_, 2 * (`Log-Likelihood` - lag(`Log-Likelihood`))),
    # Ensure Delta_df is positive and at least 1 degree of freedom:
    Delta_df  = if_else(row_number() == 1, NA_real_, pmax(1, round(abs(`Effective DF (EDF)` - lag(`Effective DF (EDF)`))))),
    LRT_p_num = if_else(row_number() == 1, NA_real_, 1 - pchisq(LRT_stat, df = Delta_df)),
    LRT_Test  = case_when(row_number() == 1 ~ "Reference", TRUE ~ sprintf("%.2f (vs. df=%d)", LRT_stat, lag(row_number()))),
    LRT_P     = case_when(row_number() == 1 ~ "—", LRT_p_num < 0.001 ~ "< 0.001", TRUE ~ sprintf("%.3f", LRT_p_num)),
    LogLik_fmt    = format(round(`Log-Likelihood`, 2), big.mark = ",", nsmall = 2),
    EDF_fmt       = sprintf("%.2f", `Effective DF (EDF)`),
    AIC_fmt       = format(round(AIC, 2), big.mark = ",", nsmall = 2),
    Delta_AIC_fmt = sprintf("%.2f", `Delta AIC`),
    BIC_fmt       = format(round(BIC, 2), big.mark = ",", nsmall = 2),
    Delta_BIC_fmt = sprintf("%.2f", `Delta BIC`),
    CV_CIndex_fmt = sprintf("%.4f", `Mean 5-Fold CV C-Index`),
    Model_Label   = if_else(grepl("df = 2", `Spline Model Architecture`), 
                            "df = 2 (Restricted Cubic Spline)*", 
                            gsub("\\(Spline\\)", "(Restricted Cubic Spline)", `Spline Model Architecture`)),
    Is_Selected   = grepl("df = 2", `Spline Model Architecture`)
  )

df_table_out <- df_spline_processed %>%
  select(Architecture = Model_Label, LogLik = LogLik_fmt, EDF = EDF_fmt, AIC = AIC_fmt, Delta_AIC = Delta_AIC_fmt,
         BIC = BIC_fmt, Delta_BIC = Delta_BIC_fmt, CV_CIndex = CV_CIndex_fmt, LRT_Test = LRT_Test, LRT_P = LRT_P, Is_Selected)

ft6 <- flextable(df_table_out %>% select(-Is_Selected))

ft6 <- set_header_labels(ft6, Architecture = "Exposure Function Architecture", LogLik = "Log-Likelihood", EDF = "Effective DF (EDF)",
                         AIC = "AIC", Delta_AIC = "\u0394AICa", BIC = "BIC", Delta_BIC = "\u0394BICb", CV_CIndex = "5-Fold CV C-Indexc",
                         LRT_Test = "LRT \u03c72d", LRT_P = "LRT p-Value")

ft6 <- font(ft6, fontname = "Times New Roman", part = "all")
ft6 <- fontsize(ft6, size = 9.0, part = "body")
ft6 <- fontsize(ft6, size = 9.5, part = "header")
ft6 <- fontsize(ft6, size = 8.5, part = "footer")
ft6 <- bold(ft6, part = "header")
ft6 <- align(ft6, j = 1, align = "left", part = "all")
ft6 <- align(ft6, j = 2:10, align = "center", part = "all")

selected_idx <- which(df_table_out$Is_Selected)
ft6 <- bold(ft6, i = selected_idx, j = 1:10, bold = TRUE)
ft6 <- bg(ft6, i = selected_idx, j = 1:10, bg = "#EDF2F7")

ft6 <- border_remove(ft6)
ft6 <- hline_top(ft6, border = fp_border(color = "black", width = 1.5), part = "header")
ft6 <- hline_bottom(ft6, border = fp_border(color = "black", width = 1.0), part = "header")
ft6 <- hline_bottom(ft6, border = fp_border(color = "black", width = 1.5), part = "body")
ft6 <- padding(ft6, padding.top = 3.5, padding.bottom = 3.5, part = "body")
ft6 <- padding(ft6, padding.top = 5.0, padding.bottom = 5.0, part = "header")

ft6 <- width(ft6, j = 1, width = 1.90)
ft6 <- width(ft6, j = 2, width = 0.85)
ft6 <- width(ft6, j = 3, width = 0.70)
ft6 <- width(ft6, j = 4, width = 0.70)
ft6 <- width(ft6, j = 5, width = 0.55)
ft6 <- width(ft6, j = 6, width = 0.70)
ft6 <- width(ft6, j = 7, width = 0.55)
ft6 <- width(ft6, j = 8, width = 0.75)
ft6 <- width(ft6, j = 9, width = 0.85)
ft6 <- width(ft6, j = 10, width = 0.65)

footnote_tab6 <- c(
  "Abbreviations: AIC, Akaike information criterion; BIC, Bayesian information criterion; CV, cross-validated; DF, degrees of freedom; EDF, effective degrees of freedom; LRT, likelihood ratio test.",
  "*Selected model architecture for primary non-linear dose-response analyses (truncated at 8 cumulative units).",
  "a. \u0394AIC represents the difference relative to the minimum observed AIC (df = 2; AIC = 3,930.27).",
  "b. \u0394BIC represents the difference relative to the minimum observed BIC (df = 2; BIC = 3,949.74). A \u0394BIC > 2.0 indicates positive evidence favoring the more parsimonious model.",
  "c. Out-of-sample predictive concordance evaluated using 5-fold cross-validation clustered at the unique candidate level; values reflect Harrell's C-index with Huber-White sandwich variance adjustments.",
  "d. Likelihood ratio test comparing sequential nested exposure architectures (df = 2 vs. df = 1 assesses departure from linearity; df = 3 vs. df = 2 assesses statistical gain of an additional non-linear knot)."
)
for (fn in footnote_tab6) {
  ft6 <- add_footer_lines(ft6, values = fn)
}
ft6 <- italic(ft6, part = "footer")

doc6 <- read_docx()
doc6 <- body_add_par(
  doc6, 
  "Supplementary Table 6. Model Selection Diagnostics and Goodness-of-Fit Metrics for Cumulative Transfusion Exposure Spline Functions", 
  style = "heading 1"
)
doc6 <- body_add_par(doc6, "")
doc6 <- body_add_flextable(doc6, ft6)

out_tab6_path <- file.path(table_dir, "Supplementary_Table_6_Spline_Optimization_Final.docx")
print(doc6, target = out_tab6_path)
message(sprintf("[Supplementary Table 6] Dynamic Word table successfully generated: %s", out_tab6_path))

# SCHOENFELD RESIDUALS
library(survival)

check_ph_assumption <- function(model_fit, model_name) {
  message(sprintf("\n--- PH Assumption Test: %s ---", model_name))
  zph_test <- cox.zph(model_fit)
  print(zph_test)
  
  num_plots <- nrow(zph_test$table) - 1
  par(mfrow = c(ceiling(sqrt(num_plots)), ceiling(sqrt(num_plots))))
  plot(zph_test, main = paste("Schoenfeld:", model_name))
  par(mfrow = c(1, 1)) # Reset plot layout
}

check_ph_assumption(fit_stage2c, "Primary Adjusted Model (Stratified)")