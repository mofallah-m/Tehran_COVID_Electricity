## ---- Pre-trend joint F-test ----
lead_coefs <- grep("^period::-(2|3|4|5|6|7):treated$", names(coef(event_total)), value = TRUE)
H <- paste(lead_coefs, "= 0")
pretrend_F <- car::linearHypothesis(event_total, hypothesis.matrix = H, test = "F")
print(pretrend_F)

## ---- Placebo event study ----
load("Data/df_reg.RData")

df_placebo <- df_reg %>%
  filter(weeks_seq >= as.Date("2019-12-01"),
         weeks_seq <  as.Date("2020-06-21")) %>%
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

event_total_placebo <- run_placebo_event_study(df_placebo, "total", placebo_shift_weeks = 6)
event_mid_placebo   <- run_placebo_event_study(df_placebo, "mid",   placebo_shift_weeks = 6)
event_peak_placebo  <- run_placebo_event_study(df_placebo, "peak",  placebo_shift_weeks = 6)
event_low_placebo   <- run_placebo_event_study(df_placebo, "low",   placebo_shift_weeks = 6)

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
p_placebo_total <- plot_placebo(event_total_placebo, ref_period = -1, title = "Placebo Test - Total",     event_line_at = 6)
p_placebo_mid   <- plot_placebo(event_mid_placebo,   ref_period = -1, title = "Placebo Test - Mid-Peak",  event_line_at = 6)
p_placebo_peak  <- plot_placebo(event_peak_placebo,  ref_period = -1, title = "Placebo Test - Peak",      event_line_at = 6)
p_placebo_low   <- plot_placebo(event_low_placebo,   ref_period = -1, title = "Placebo Test - Off-Peak",  event_line_at = 6)
ggsave("Figures/placebo_total.pdf", p_placebo_total, width = 7, height = 5)
ggsave("Figures/placebo_mid.pdf",   p_placebo_mid,   width = 7, height = 5)
ggsave("Figures/placebo_peak.pdf",  p_placebo_peak,  width = 7, height = 5)
ggsave("Figures/placebo_low.pdf",   p_placebo_low,   width = 7, height = 5)

# 1) Bundle the placebo ES models
models_es_placebo <- list(
  "Total"    = event_total_placebo,
  "Mid-Peak" = event_mid_placebo,
  "Peak"     = event_peak_placebo,
  "Off-Peak" = event_low_placebo
)

# 2) Build the placebo ES coefficients table
make_es_table(
  models_named = models_es_placebo,
  caption = "Placebo Event-Study Coefficients by Relative Week (Treatment shifted by 6 weeks)",
  label   = "tab:eventstudy_placebo",
  outfile = "Tables/eventstudy_placebo.tex",
  ref_period = -1,          # adjust if your ref differs
  # window = c(-6, 15)      # optional: trim displayed horizon
)

## ---- Event study (expanded window) ----
df_es_ex <- df_reg %>%
  select(
    weeks_seq, id, power_grid_location, month, year,
    total_treated = weekly_total_cons, total_control = weekly_cf_total,
    mid_treated   = weekly_mid_cons,   mid_control   = weekly_cf_mid,
    peak_treated  = weekly_peak_cons,  peak_control  = weekly_cf_peak,
    low_treated   = weekly_low_cons,   low_control  = weekly_cf_low
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

event_total <- run_event_study(df_es_ex, "total")
event_mid   <- run_event_study(df_es_ex, "mid")
event_peak  <- run_event_study(df_es_ex, "peak")
event_low   <- run_event_study(df_es_ex, "low")

save_event_plot(event_total, -1, "Event Study (Expanded Periods) - Total",    "Figures/event_ex_total.pdf")
save_event_plot(event_mid,   -1, "Event Study (Expanded Periods) - Mid-Peak", "Figures/event_ex_mid.pdf")
save_event_plot(event_peak,  -1, "Event Study (Expanded Periods) - Peak",     "Figures/event_ex_peak.pdf")
save_event_plot(event_low,   -1, "Event Study (Expanded Periods) - Off-Peak", "Figures/event_ex_low.pdf")

# bundle models (expanded sample)
models_es_ex <- list(
  "Total"    = event_total,
  "Mid-Peak" = event_mid,
  "Peak"     = event_peak,
  "Off-Peak" = event_low
)

# table of ES coefficients (full window)
make_es_table(
  models_named = models_es_ex,
  caption = "Event-Study Coefficients by Relative Week (Expanded Sample)",
  label   = "tab:eventstudy_coefs_expanded",
  outfile = "Tables/eventstudy_coefs_expanded.tex",
  ref_period = -1,       # set to your ES reference (e.g., -1)
  window = c(-13, 17)   # optional: display a tighter window
)

## ---- Event study by location ----
p_loc_total <- plot_event_study_by_location(df_es, "total", ref_period = -1)
ggsave("Figures/event_loc_total.pdf", p_loc_total, width = 8.5, height = 11, units = "in", dpi = 300)
p_loc_mid <- plot_event_study_by_location(df_es, "mid", ref_period = -1)
ggsave("Figures/event_loc_mid.pdf",   p_loc_mid,  width = 8.5, height = 11, units = "in", dpi = 300)
p_loc_peak <- plot_event_study_by_location(df_es, "peak", ref_period = -1)
ggsave("Figures/event_loc_peak.pdf",  p_loc_peak, width = 8.5, height = 11, units = "in", dpi = 300)
p_loc_low <- plot_event_study_by_location(df_es, "low", ref_period = -1)
ggsave("Figures/event_loc_low.pdf",   p_loc_low,  width = 8.5, height = 11, units = "in", dpi = 300)

## ---- Mean-control (2018, 2019, 2021) event study ----
load("Data/weekly_coefs.RData")  

df_weekly_alt <- dt %>%
  dplyr::filter(week != 53,
                lubridate::year(current_record) < 2022,
                lubridate::year(current_record) > 2017) %>%
  dplyr::left_join(weekly_coefs %>% dplyr::select(week, estimate), by = "week") %>%
  dplyr::mutate(weight = ifelse(is.na(estimate), 1, exp(estimate))) %>%
  dplyr::group_by(id, previous_record, current_record) %>%
  dplyr::mutate(
    weekly_total_cons = (weight / sum(weight, na.rm = TRUE)) * total_cons,
    weekly_mid_cons   = (weight / sum(weight, na.rm = TRUE)) * mid_cons,
    weekly_peak_cons  = (weight / sum(weight, na.rm = TRUE)) * peak_cons,
    weekly_low_cons   = (weight / sum(weight, na.rm = TRUE)) * low_cons
  ) %>% dplyr::ungroup() %>% dplyr::select(-estimate, -weight) %>%
  dplyr::group_by(id, usage, weeks_seq, week, power_grid_location) %>%
  dplyr::summarise(
    weekly_total_cons = mean(weekly_total_cons, na.rm = TRUE),
    weekly_mid_cons   = mean(weekly_mid_cons,   na.rm = TRUE),
    weekly_peak_cons  = mean(weekly_peak_cons,  na.rm = TRUE),
    weekly_low_cons   = mean(weekly_low_cons,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(year = lubridate::year(weeks_seq),
                week = lubridate::isoweek(weeks_seq))

weekly_means_cf <- df_weekly_alt %>%
  dplyr::filter(year != 2020) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(
    weekly_cf_total = mean(log(weekly_total_cons + 1), na.rm = TRUE),
    weekly_cf_mid   = mean(log(weekly_mid_cons   + 1), na.rm = TRUE),
    weekly_cf_peak  = mean(log(weekly_peak_cons  + 1), na.rm = TRUE),
    weekly_cf_low   = mean(log(weekly_low_cons   + 1), na.rm = TRUE),
    .groups = "drop"
  )

df_weekly_alt <- df_weekly_alt %>%
  dplyr::left_join(weekly_means_cf, by = c("week"))

df_reg_alt <- df_weekly_alt %>%
  dplyr::select(
    id, usage, weeks_seq, week, power_grid_location,
    weekly_total_cons, weekly_mid_cons, weekly_peak_cons, weekly_low_cons,
    weekly_cf_total,  weekly_cf_mid,   weekly_cf_peak,   weekly_cf_low
  ) %>%
  dplyr::mutate(dplyr::across(ends_with("_cons"), log1p)) %>%
  dplyr::filter(weeks_seq >= as.Date("2019-03-01"),
                weeks_seq <  as.Date("2021-01-01")) %>%
  dplyr::mutate(
    month = lubridate::month(weeks_seq),
    year  = lubridate::year(weeks_seq),
    post  = ifelse(weeks_seq >= as.Date("2020-03-01"), 1, 0)
  )

df_es_alt <- df_reg_alt %>%
  dplyr::filter(weeks_seq >= as.Date("2020-01-11"),
                weeks_seq <  as.Date("2020-06-21")) %>%
  dplyr::select(
    weeks_seq, id, power_grid_location, month, year,
    total_treated = weekly_total_cons, total_control = weekly_cf_total,
    mid_treated   = weekly_mid_cons,   mid_control   = weekly_cf_mid,
    peak_treated  = weekly_peak_cons,  peak_control  = weekly_cf_peak,
    low_treated   = weekly_low_cons,   low_control   = weekly_cf_low
  ) %>%
  dplyr::mutate(quarter = ceiling(month / 3)) %>%
  tidyr::pivot_longer(-c(weeks_seq, id, power_grid_location, month, quarter, year),
                      names_to = c(".value","group"), names_sep = "_") %>%
  dplyr::mutate(
    treated = as.integer(group == "treated"),
    period  = as.integer(difftime(weeks_seq, as.Date("2020-03-01"), units = "weeks")),
    week    = lubridate::isoweek(weeks_seq)
  ) %>%
  dplyr::filter(!(is.na(total) & is.na(mid) & is.na(peak) & is.na(low))) %>%
  dplyr::filter(period >= -20, period < 29)

event_total_alt <- run_event_study(df_es_alt, "total")
event_mid_alt   <- run_event_study(df_es_alt, "mid")
event_peak_alt  <- run_event_study(df_es_alt, "peak")
event_low_alt   <- run_event_study(df_es_alt, "low")

save_event_plot(event_total_alt, -1, "Event Study (2018, 2019, and 2021 Mean Control) - Total",    "Figures/event_mean_total.pdf")
save_event_plot(event_mid_alt,   -1, "Event Study (2018, 2019, and 2021 Mean Control) - Mid-Peak", "Figures/event_mean_mid.pdf")
save_event_plot(event_peak_alt,  -1, "Event Study (2018, 2019, and 2021 Mean Control) - Peak",     "Figures/event_mean_peak.pdf")
save_event_plot(event_low_alt,   -1, "Event Study (2018, 2019, and 2021 Mean Control) - Off-Peak", "Figures/event_mean_low.pdf")


models_es <- list(
  "Total"    = event_total_alt,
  "Mid-Peak" = event_mid_alt,
  "Peak"     = event_peak_alt,
  "Off-Peak" = event_low_alt
)

make_es_table(
  models_named = models_es,
  caption = "Event-Study Coefficients by Relative Week (Control: Mean of 2018, 2019, and 2021)",
  label   = "tab:eventstudy_coefs_altcontrol",
  outfile = "Tables/eventstudy_coefs_altcontrol.tex",
  ref_period = -1
)