# ==============================================================================
# script 3: sub-analyses 
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(broom)
})

root_dir   <- "/Users/zvr/Downloads/Clinical2.0"
output_dir <- file.path(root_dir, "processed_data")
table_dir  <- file.path(root_dir, "tables_output")
fig_dir    <- file.path(root_dir, "figures_output")

df_master_base    <- readRDS(file.path(output_dir, "df_master_base.rds"))
df_trans_products <- readRDS(file.path(output_dir, "df_trans_products.rds"))
df_trans_clean    <- readRDS(file.path(output_dir, "df_trans_clean.rds"))
df_screen_clean   <- readRDS(file.path(output_dir, "df_screen_clean.rds"))

# ==============================================================================
# 1. blood products
# ==============================================================================

df_long_base <- tmerge(data1 = df_master_base, data2 = df_master_base, id = pat_id, tstart = 0, tstop = end_day)
df_long_base <- tmerge(data1 = df_long_base, data2 = df_master_base, id = pat_id, event = event(end_day, status))

df_long_prod <- tmerge(
  data1      = df_long_base,
  data2      = df_trans_products,
  id         = pat_id,
  exp_rbc    = tdc(trans_day, has_rbc),
  exp_plasma = tdc(trans_day, has_plasma)
) %>%
  mutate(
    exp_rbc    = coalesce(exp_rbc, 0),
    exp_plasma = coalesce(exp_plasma, 0)
  )


fit_prod_un_rbc <- coxph(Surv(tstart, tstop, event) ~ exp_rbc, data = df_long_prod, cluster = pat_id)
fit_prod_un_pla <- coxph(Surv(tstart, tstop, event) ~ exp_plasma, data = df_long_prod, cluster = pat_id)

# Stratified Multivariable Model
fit_prod_adj <- coxph(
  Surv(tstart, tstop, event) ~ exp_rbc + exp_plasma + age_at_listing + 
    sex_pregnancy_status + strata(has_prior_transplant, organ_group),
  data = df_long_prod, 
  cluster = pat_id
)

res_rbc_un  <- tidy(fit_prod_un_rbc, exponentiate = TRUE, conf.int = TRUE)
res_pla_un  <- tidy(fit_prod_un_pla, exponentiate = TRUE, conf.int = TRUE)
res_adj     <- tidy(fit_prod_adj, exponentiate = TRUE, conf.int = TRUE)

cat("\n--- BLOOD PRODUCT SUBANALYSIS (RBC vs PLASMA) ---\n")
cat(sprintf("RBC Unadjusted HR: %.2f (%.2f–%.2f), p = %.3f\n", res_rbc_un$estimate[1], res_rbc_un$conf.low[1], res_rbc_un$conf.high[1], res_rbc_un$p.value[1]))
cat(sprintf("Plasma Unadjusted HR: %.2f (%.2f–%.2f), p = %.3f\n", res_pla_un$estimate[1], res_pla_un$conf.low[1], res_pla_un$conf.high[1], res_pla_un$p.value[1]))

adj_rbc <- res_adj %>% filter(term == "exp_rbc")
adj_pla <- res_adj %>% filter(term == "exp_plasma")
cat(sprintf("RBC Adjusted HR: %.2f (%.2f–%.2f), p = %.3f\n", adj_rbc$estimate, adj_rbc$conf.low, adj_rbc$conf.high, adj_rbc$p.value))
cat(sprintf("Plasma Adjusted HR: %.2f (%.2f–%.2f), p = %.3f\n", adj_pla$estimate, adj_pla$conf.low, adj_pla$conf.high, adj_pla$p.value))

supp_tab_components <- tibble(
  `Blood Product Component` = c("Red Blood Cells (RBC)", "Plasma (FFP)"),
  `Unadjusted HR (95% CI)`  = c(
    sprintf("%.2f (%.2f–%.2f)", res_rbc_un$estimate[1], res_rbc_un$conf.low[1], res_rbc_un$conf.high[1]),
    sprintf("%.2f (%.2f–%.2f)", res_pla_un$estimate[1], res_pla_un$conf.low[1], res_pla_un$conf.high[1])
  ),
  `Unadjusted p-Value`      = c(
    if (res_rbc_un$p.value[1] < 0.001) "< 0.001" else sprintf("%.3f", res_rbc_un$p.value[1]),
    if (res_pla_un$p.value[1] < 0.001) "< 0.001" else sprintf("%.3f", res_pla_un$p.value[1])
  ),
  `Mutually Adjusted aHR (95% CI)*` = c(
    sprintf("%.2f (%.2f–%.2f)", adj_rbc$estimate, adj_rbc$conf.low, adj_rbc$conf.high),
    sprintf("%.2f (%.2f–%.2f)", adj_pla$estimate, adj_pla$conf.low, adj_pla$conf.high)
  ),
  `Adjusted p-Value`        = c(
    if (adj_rbc$p.value < 0.001) "< 0.001" else sprintf("%.3f", adj_rbc$p.value),
    if (adj_pla$p.value < 0.001) "< 0.001" else sprintf("%.3f", adj_pla$p.value)
  )
)

write_csv(supp_tab_components, file.path(table_dir, "Supplementary_Table_Blood_Components.csv"))
cat("\n--- BLOOD PRODUCT COMPONENT TABLE ---\n")
print(as.data.frame(supp_tab_components), row.names = FALSE)

# ==============================================================================
# 2. antibody durability
# ==============================================================================
df_sensitized <- df_master_base %>% filter(status == 1) %>%
  mutate(trans_prior = if_else(pat_id %in% df_trans_clean$pat_id, 1, 0))

df_post_screens <- df_sensitized %>%
  select(pat_id, first_pos_date) %>%
  left_join(df_screen_clean, by = "pat_id", relationship = "many-to-many") %>%
  filter(!is.na(draw_date_clean) & draw_date_clean > first_pos_date)

df_fu_agg <- df_post_screens %>%
  group_by(pat_id) %>%
  summarize(
    has_followup = n() > 0,
    fu_duration  = as.numeric(difftime(max(draw_date_clean), first_pos_date[1], units = "days")),
    clear_date   = suppressWarnings(min(draw_date_clean[screen_pos == 0], na.rm = TRUE)),
    .groups      = "drop"
  ) %>%
  mutate(
    clear_date = if_else(is.infinite(clear_date), NA_Date_, clear_date)
  )

df_revert_status <- df_sensitized %>%
  left_join(df_fu_agg, by = "pat_id") %>%
  mutate(
    has_followup = coalesce(has_followup, FALSE),
    cleared_flag = if_else(!is.na(clear_date), 1, 0),
    clear_days   = if_else(cleared_flag == 1, 
                           as.numeric(difftime(clear_date, first_pos_date, units = "days")), 
                           fu_duration),
    pos_at_90 = case_when(
      !has_followup | coalesce(fu_duration, 0) < 90 ~ NA_integer_,
      cleared_flag == 1 & clear_days <= 90 ~ 0,
      TRUE ~ 1
    ),
    pos_at_180 = case_when(
      !has_followup | coalesce(fu_duration, 0) < 180 ~ NA_integer_,
      cleared_flag == 1 & clear_days <= 180 ~ 0,
      TRUE ~ 1
    ),
    pos_at_365 = case_when(
      !has_followup | coalesce(fu_duration, 0) < 365 ~ NA_integer_,
      cleared_flag == 1 & clear_days <= 365 ~ 0,
      TRUE ~ 1
    )
  )

df_kinetics <- df_revert_status %>% filter(has_followup)

fit_km_persist <- survfit(Surv(clear_days, cleared_flag) ~ trans_prior, data = df_kinetics)

safe_km <- function(fit, grp_val) {
  s <- summary(fit, times = 365, extend = TRUE)
  idx <- which(as.character(s$strata) == paste0("trans_prior=", grp_val))
  if (length(idx) == 1 && length(s$surv) >= idx && !is.na(s$surv[idx])) {
    return(sprintf("%.1f%% (%.1f%%–%.1f%%)", s$surv[idx]*100, s$lower[idx]*100, s$upper[idx]*100))
  } else {
    return("—")
  }
}

km_exp_str   <- safe_km(fit_km_persist, 1)
km_unexp_str <- safe_km(fit_km_persist, 0)

sdiff <- survdiff(Surv(clear_days, cleared_flag) ~ trans_prior, data = df_kinetics)
p_logrank <- pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
p_logrank_str <- if(length(p_logrank) > 0 && !is.na(p_logrank)) {
  if(p_logrank < 0.001) "< 0.001" else sprintf("%.3f", p_logrank)
} else { "—" }

fit_cox_revert <- coxph(
  Surv(clear_days, cleared_flag) ~ trans_prior + age_at_listing + 
    sex_pregnancy_status + strata(has_prior_transplant, organ_group), 
  data = df_kinetics
)

s_rev <- summary(fit_cox_revert)
if (!is.null(s_rev$conf.int) && nrow(s_rev$conf.int) >= 1) {
  hr_est <- s_rev$conf.int[1, "exp(coef)"]
  hr_low <- s_rev$conf.int[1, "lower .95"]
  hr_upp <- s_rev$conf.int[1, "upper .95"]
  p_val  <- s_rev$coefficients[1, 5]
  revert_hr_str <- sprintf("%.2f (%.2f–%.2f)", hr_est, hr_low, hr_upp)
  revert_p_str  <- if(!is.na(p_val) && p_val < 0.001) "< 0.001" else sprintf("%.3f", p_val)
} else {
  revert_hr_str <- "Not estimable (0 events)"
  revert_p_str  <- "—"
}



fmt_pct <- function(n, total) {
  if (is.na(total) || total == 0) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", n, 100 * (n / total))
}

safe_fisher <- function(var1, var2) {
  tbl <- table(var1, var2)
  if (nrow(tbl) >= 2 && ncol(tbl) >= 2) {
    return(sprintf("%.3f", fisher.test(tbl)$p.value))
  } else {
    return("—") 
  }
}

d_all_e <- filter(df_revert_status, trans_prior == 1)
d_kin_e <- filter(df_kinetics, trans_prior == 1)
n_tot_e <- nrow(d_all_e)

d_all_u <- filter(df_revert_status, trans_prior == 0)
d_kin_u <- filter(df_kinetics, trans_prior == 0)
n_tot_u <- nrow(d_all_u)

col_metric <- c(
  "Post-Sensitization Follow-up Duration Analysis 1",
  "  Total Sensitization Events",
  "  Excluded: No post-sensitization follow-up testing 2",
  "  Filtered Analysis Cohort",
  "  Post-Sensitization Follow-up Duration 3",
  "Periodic Analysis of Antibody Persistence 4",
  "  Excluded: post-sensitization follow-up testing < 90 days 2",
  "  Excluded: post-sensitization follow-up testing < 180 days 2",
  "  Excluded: post-sensitization follow-up testing <= 1 year 2",
  "  Post-Sensitization Seropositivity",
  "    Retention at >= 90 Days 5",
  "    Retention at >= 180 Days 5",
  "    Retention at >= 1 Year 5",
  "Longitudinal Antibody Clearance Analyses",
  "  1-Year Persistent Retention Probability (KM) 6",
  "  Adjusted Antibody Clearance Hazard (aHR) 7"
)

col_exp <- c(
  "",
  as.character(n_tot_e),
  fmt_pct(sum(!d_all_e$has_followup), n_tot_e),
  fmt_pct(sum(d_all_e$has_followup), n_tot_e),
  sprintf("%.0f [%.0f–%.0f]", median(d_kin_e$fu_duration, na.rm=TRUE), quantile(d_kin_e$fu_duration, 0.25, na.rm=TRUE), quantile(d_kin_e$fu_duration, 0.75, na.rm=TRUE)),
  "",
  fmt_pct(sum(coalesce(d_all_e$fu_duration, 0) < 90), n_tot_e),
  fmt_pct(sum(coalesce(d_all_e$fu_duration, 0) < 180), n_tot_e),
  fmt_pct(sum(coalesce(d_all_e$fu_duration, 0) < 365), n_tot_e),
  "",
  fmt_pct(sum(d_kin_e$pos_at_90 == 1, na.rm=TRUE), sum(!is.na(d_kin_e$pos_at_90))),
  fmt_pct(sum(d_kin_e$pos_at_180 == 1, na.rm=TRUE), sum(!is.na(d_kin_e$pos_at_180))),
  fmt_pct(sum(d_kin_e$pos_at_365 == 1, na.rm=TRUE), sum(!is.na(d_kin_e$pos_at_365))),
  "",
  km_exp_str,
  revert_hr_str
)

col_unexp <- c(
  "",
  as.character(n_tot_u),
  fmt_pct(sum(!d_all_u$has_followup), n_tot_u),
  fmt_pct(sum(d_all_u$has_followup), n_tot_u),
  sprintf("%.0f [%.0f–%.0f]", median(d_kin_u$fu_duration, na.rm=TRUE), quantile(d_kin_u$fu_duration, 0.25, na.rm=TRUE), quantile(d_kin_u$fu_duration, 0.75, na.rm=TRUE)),
  "",
  fmt_pct(sum(coalesce(d_all_u$fu_duration, 0) < 90), n_tot_u),
  fmt_pct(sum(coalesce(d_all_u$fu_duration, 0) < 180), n_tot_u),
  fmt_pct(sum(coalesce(d_all_u$fu_duration, 0) < 365), n_tot_u),
  "",
  fmt_pct(sum(d_kin_u$pos_at_90 == 1, na.rm=TRUE), sum(!is.na(d_kin_u$pos_at_90))),
  fmt_pct(sum(d_kin_u$pos_at_180 == 1, na.rm=TRUE), sum(!is.na(d_kin_u$pos_at_180))),
  fmt_pct(sum(d_kin_u$pos_at_365 == 1, na.rm=TRUE), sum(!is.na(d_kin_u$pos_at_365))),
  "",
  km_unexp_str,
  "1.00 (Reference)"
)

col_pval <- c(
  "", "—", "—", "—",
  sprintf("%.3f", wilcox.test(fu_duration ~ trans_prior, data = df_kinetics)$p.value),
  "", "—", "—", "—", "—",
  safe_fisher(df_kinetics$pos_at_90, df_kinetics$trans_prior),
  safe_fisher(df_kinetics$pos_at_180, df_kinetics$trans_prior),
  safe_fisher(df_kinetics$pos_at_365, df_kinetics$trans_prior),
  "",
  p_logrank_str,
  revert_p_str
)

supp_tab_5 <- tibble(
  `Surveillance Metric` = col_metric,
  `Exposed` = col_exp,
  `Unexposed` = col_unexp,
  `p-Value` = col_pval
)

write_csv(supp_tab_5, file.path(table_dir, "Supplementary_Table_5_Antibody_Kinetics_2.csv"))

# --- proportional hazards check ---
check_ph_assumption <- function(model_fit, model_name) {
  message(sprintf("\n--- PH Assumption Test: %s ---", model_name))
  zph_test <- cox.zph(model_fit)
  print(zph_test)
  
  num_plots <- nrow(zph_test$table) - 1
  par(mfrow = c(ceiling(sqrt(num_plots)), ceiling(sqrt(num_plots))))
  plot(zph_test, main = paste("Schoenfeld:", model_name))
  par(mfrow = c(1, 1)) # Reset plot layout
}

check_ph_assumption(fit_prod_adj, "Blood Product Subanalysis (Stratified)")
check_ph_assumption(fit_cox_revert, "Antibody Clearance Kinetics (Stratified)")
