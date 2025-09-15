## ---- Regression-ready panel ----
load("Data/df_weekly.RData")

# 1) Compute pre-COVID stats on log1p scale
prewin_stats <- df_weekly %>%
  filter(weeks_seq >  as.Date("2019-03-01"),
         weeks_seq <  as.Date("2020-03-01")) %>%
  mutate(cons_log = log1p(weekly_total_cons)) %>%   # safe for zeros
  group_by(id) %>%
  summarise(
    n_finite = sum(is.finite(cons_log)),
    pre_mean = ifelse(n_finite >= 1, mean(cons_log[is.finite(cons_log)]), NA_real_),
    pre_var  = ifelse(n_finite >= 2, var(cons_log[is.finite(cons_log)]), 0),  # 0 if <2 obs
    .groups  = "drop"
  )

# 2) Assign terciles (1 = Low, 2 = Medium, 3 = High)
prewin_flags <- prewin_stats %>%
  mutate(
    con_group = ntile(pre_mean, 3),
    var_group = ntile(pre_var,  3)
  ) %>%
  mutate(
    con_group = factor(con_group, labels = c("Low", "Medium", "High")),
    var_group = factor(var_group, labels = c("Low", "Medium", "High"))
  ) %>%
  select(id, con_group, var_group)

# 3) Build df_reg and join groups
df_reg <- df_weekly %>%
  left_join(prewin_flags, by = "id") %>%
  select(
    id, usage, weeks_seq, week, power_grid_location,
    weekly_total_cons, weekly_mid_cons, weekly_peak_cons, weekly_low_cons,
    weekly_cf_total, weekly_cf_mid, weekly_cf_peak, weekly_cf_low,
    con_group, var_group
  ) %>%
  mutate(
    across(starts_with("weekly_"), log1p)   # logs for outcomes
  ) %>%
  filter(
    weeks_seq >= as.Date("2019-03-01"),
    weeks_seq <  as.Date("2021-01-01")
  ) %>%
  mutate(
    month = month(weeks_seq),
    year  = year(weeks_seq),
    post  = if_else(weeks_seq >= as.Date("2020-03-01"), 1L, 0L)
  )


save(df_reg, file = "Data/df_reg.RData")

## ---- DiD models ----

## ---- Main DiD ----
df_did <- df_reg %>%
  filter(weeks_seq >= as.Date("2020-01-11"), weeks_seq < as.Date("2020-06-21")) %>%
  select(
    weeks_seq, id, power_grid_location, month, year,
    total_treated = weekly_total_cons, total_control = weekly_cf_total,
    mid_treated   = weekly_mid_cons,   mid_control   = weekly_cf_mid,
    peak_treated  = weekly_peak_cons,  peak_control  = weekly_cf_peak,
    low_treated   = weekly_low_cons,   low_control   = weekly_cf_low,
    con_group, var_group
  ) %>%
  mutate(quarter = ceiling(month / 3)) %>%
  pivot_longer(-c(weeks_seq, id, power_grid_location, month, quarter, year, con_group, var_group),
               names_to = c(".value","group"), names_sep = "_") %>%
  mutate(
    treated = as.integer(group == "treated"),
    period  = as.integer(difftime(weeks_seq, as.Date("2020-03-01"), units = "weeks")),
    year    = if_else(group == "control", 2019L, year),
    week    = isoweek(weeks_seq),
    post    = if_else(weeks_seq >= as.Date("2020-03-01"), 1L, 0L)
  ) %>%
  filter(!(is.na(total) & is.na(mid) & is.na(peak) & is.na(low)))

# 1) Run models (you already did)
did_total <- run_did_model(df_did, "total")
did_mid   <- run_did_model(df_did, "mid")
did_peak  <- run_did_model(df_did, "peak")
did_low   <- run_did_model(df_did, "low")

models <- list(
  "Total"    = did_total,
  "Mid-Peak" = did_mid,
  "Peak"     = did_peak,
  "Off-Peak" = did_low
)

# show only Num. Obs. and R2
gof_map <- data.frame(
  raw   = c("nobs", "r.squared"),
  clean = c("Num.Obs.", "R2"),
  fmt   = c(0, 3),
  omit  = FALSE,
  stringsAsFactors = FALSE
)

latex_main <- modelsummary(
  models,
  output = "latex",
  stars = TRUE,
  coef_rename = c(
    "treated"       = "Treated",
    "post"          = "Post",
    "treated:post"  = "Treated × Post"
  ),
  statistic = "({std.error})",
  statistic_vertical = TRUE,   # keep your current layout
  gof_map = gof_map,           # << only Num.Obs. & R2 will be shown
  title = paste(
    "DiD Estimates of the Impact of COVID-19 on Electricity Consumption",
    "\\label{tab:did_results}"
  ),
  escape = FALSE,
  align = "lcccc",                    # one col for term + 4 models
  add_header_above = c(" " = 1, "Consumption Type" = 4),
  latex_options = c("hold_position")
)

# stretch to full width with tabular*
latex_main <- sub(
  "\\\\begin\\{tabular\\}\\[t\\]\\{lcccc\\}",
  "\\\\begin{tabular*}{\\\\textwidth}{@{\\\\extracolsep{\\\\fill}}lcccc}",
  latex_main
)
latex_main <- sub("\\\\end\\{tabular\\}", "\\\\end{tabular*}", latex_main)

# tighter columns
latex_main <- sub("\\\\centering", "\\\\centering\n\\\\setlength{\\\\tabcolsep}{6pt}", latex_main, fixed = TRUE)

# write out
cat(latex_main, file = "Tables/did_results.tex")


# ---- Consumption Based and Variance Based DiD ----

df_high_cons_did <- df_did %>%
  filter(con_group == "High")

df_med_cons_did <- df_did %>%
  filter(con_group == "Medium")

df_low_cons_did <- df_did %>%
  filter(con_group == "Low")

df_high_var_did <- df_did %>%
  filter(var_group == "High")

df_med_var_did <- df_did %>%
  filter(var_group == "Medium")

df_low_var_did <- df_did %>%
  filter(var_group == "Low")

did_high_cons <- run_did_model(df_high_cons_did, "total")
did_med_cons <- run_did_model(df_med_cons_did, "total")
did_low_cons <- run_did_model(df_low_cons_did, "total")

did_high_var <- run_did_model(df_high_var_did, "total")
did_med_var <- run_did_model(df_med_var_did, "total")
did_low_var <- run_did_model(df_low_var_did, "total")


dir.create("Tables", showWarnings = FALSE, recursive = TRUE)


## =========================
## A) By baseline consumption terciles
## =========================
models_cons <- list(
  "Low"    = did_low_cons,
  "Medium" = did_med_cons,
  "High"   = did_high_cons
)

# Show only Num.Obs. and R2 in GOF
gof_map_cons <- data.frame(
  raw   = c("nobs", "r.squared"),
  clean = c("Num.Obs.", "R2"),
  fmt   = c(0, 3),
  omit  = FALSE,
  stringsAsFactors = FALSE
)

latex_str <- modelsummary(
  models_cons,
  output = "latex",
  stars = TRUE,
  coef_rename = c(
    "treated"      = "Treated",
    "post"         = "Post",
    "treated:post" = "Treated × Post"
  ),
  statistic = "({std.error})",
  statistic_vertical = TRUE,
  gof_map = gof_map_cons,   # << keep only Num.Obs. & R2
  title = paste(
    "DiD Estimates by Baseline Consumption Intensity",
    "\\label{tab:did_by_consumption_terciles}"
  ),
  escape = FALSE,
  align = "lccc",
  add_header_above = c(" " = 1, "Baseline Consumption Intensity" = 3),
  latex_options = c("hold_position")
)

## Make it full-width (tabular* with elastic spacing)
latex_str <- sub(
  "\\\\begin\\{tabular\\}\\[t\\]\\{lccc\\}",
  "\\\\begin{tabular*}{\\\\textwidth}{@{\\\\extracolsep{\\\\fill}}lccc}",
  latex_str
)
latex_str <- sub("\\\\end\\{tabular\\}", "\\\\end{tabular*}", latex_str)

## Slightly tighter columns
latex_str <- sub("\\\\centering", "\\\\centering\n\\\\setlength{\\\\tabcolsep}{6pt}", latex_str, fixed = TRUE)

## Write to file
cat(latex_str, file = "Tables/did_by_consumption_terciles.tex")


## =========================
## B) By variance (stability) terciles
## =========================
models_var <- list(
  "Low"    = did_low_var,
  "Medium" = did_med_var,
  "High"   = did_high_var
)


# Show only Num.Obs. and R2
gof_map_var <- data.frame(
  raw   = c("nobs", "r.squared"),
  clean = c("Num.Obs.", "R2"),
  fmt   = c(0, 3),
  omit  = FALSE,
  stringsAsFactors = FALSE
)

# 1) Generate LaTeX code for variance-terciles table
latex_var <- modelsummary(
  models_var,
  output = "latex",
  stars = TRUE,
  coef_rename = c(
    "treated"      = "Treated",
    "post"         = "Post",
    "treated:post" = "Treated × Post"
  ),
  statistic = "({std.error})",
  statistic_vertical = TRUE,
  gof_map = gof_map_var,   # << only Num.Obs. & R2
  title = paste(
    "DiD Estimates by Pre-COVID Consumption Stability (Variance Terciles)",
    "\\label{tab:did_by_variance_terciles}"
  ),
  escape = FALSE,
  align = "lccc",
  add_header_above = c(" " = 1, "Variance Group" = 3),
  latex_options = c("hold_position")
)

# 2) Stretch to full width (tabular* with elastic spacing)
latex_var <- sub(
  "\\\\begin\\{tabular\\}\\[t\\]\\{lccc\\}",
  "\\\\begin{tabular*}{\\\\textwidth}{@{\\\\extracolsep{\\\\fill}}lccc}",
  latex_var
)
latex_var <- sub("\\\\end\\{tabular\\}", "\\\\end{tabular*}", latex_var)

# Optional: tighten inter-column padding
latex_var <- sub("\\\\centering", "\\\\centering\n\\\\setlength{\\\\tabcolsep}{6pt}", latex_var, fixed = TRUE)

# 3) Write the final LaTeX table to file
cat(latex_var, file = "Tables/did_by_variance_terciles.tex")


## ---- Event study (main) ----
df_es <- df_reg %>%
  filter(weeks_seq >= as.Date("2020-01-11"), weeks_seq < as.Date("2020-06-21")) %>%
  select(
    weeks_seq, id, power_grid_location, month, year,
    total_treated = weekly_total_cons, total_control = weekly_cf_total,
    mid_treated   = weekly_mid_cons,   mid_control   = weekly_cf_mid,
    peak_treated  = weekly_peak_cons,  peak_control  = weekly_cf_peak,
    low_treated   = weekly_low_cons,   low_control   = weekly_cf_low
  ) %>%
  mutate(quarter = ceiling(month / 3)) %>%
  pivot_longer(-c(weeks_seq, id, power_grid_location, month, quarter, year),
               names_to = c(".value","group"), names_sep = "_") %>%
  mutate(
    treated = as.integer(group == "treated"),
    period  = as.integer(difftime(weeks_seq, as.Date("2020-03-01"), units = "weeks")),
    week    = isoweek(weeks_seq)
  ) %>%
  filter(!(is.na(total) & is.na(mid) & is.na(peak) & is.na(low))) %>%
  filter(period >= -20, period < 29)

event_total <- run_event_study(df_es, "total")
event_mid   <- run_event_study(df_es, "mid")
event_peak  <- run_event_study(df_es, "peak")
event_low   <- run_event_study(df_es, "low")

# ----------------- USAGE with your models ----------------------------------

models_es <- list(
  "Total"    = event_total,
  "Mid-Peak" = event_mid,
  "Peak"     = event_peak,
  "Off-Peak" = event_low
)

make_es_table(
  models_named = models_es,
  caption = "Event-Study Coefficients",
  label   = "tab:eventstudy_coefs",
  outfile = "Tables/eventstudy_coefs.tex",
  ref_period = -1,
)



dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
save_event_plot(event_total, -1, "Event Study - Total",    "Figures/event_total.pdf")
save_event_plot(event_mid,   -1, "Event Study - Mid-Peak", "Figures/event_mid.pdf")
save_event_plot(event_peak,  -1, "Event Study - Peak",     "Figures/event_peak.pdf")
save_event_plot(event_low,   -1, "Event Study - Off-Peak", "Figures/event_low.pdf")