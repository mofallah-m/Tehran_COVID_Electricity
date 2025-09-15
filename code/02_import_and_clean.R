## ---- Working directory & IO ----
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Data", showWarnings = FALSE, recursive = TRUE)
dir.create("Tables",  showWarnings = FALSE, recursive = TRUE)

## ---- Load raw & combine ----
load("Data/Data97to01.RData")
load("Data/Data94to98.RData")
Data <- bind_rows(Data94to98, Data97to01)

## ---- Pre-COVID weekly FE & smooth ----
PreCovid <- Data %>%
  mutate(
    previous_record = as.Date(previous_record),
    current_record  = as.Date(current_record),
    row_id          = dplyr::row_number()
  ) %>%
  filter(current_record <= as.Date("2020-03-10"))

PreCovid_dt <- as.data.table(PreCovid)
PreCovid_dt[, `:=`(
  week_start = isoweek(as.Date(previous_record)),
  week_end   = isoweek(as.Date(current_record) - 1)
)]
n_rows   <- nrow(PreCovid_dt)
week_mat <- matrix(0L, nrow = n_rows, ncol = 53)
for (i in seq_len(n_rows)) {
  sw <- PreCovid_dt$week_start[i]; ew <- PreCovid_dt$week_end[i]
  if (!is.na(sw) && !is.na(ew) && sw <= ew) week_mat[i, sw:ew] <- 1L
}
week_cols   <- paste0("week_", 1:53)
PreCovid_dt <- cbind(PreCovid_dt, as.data.table(week_mat))
setnames(PreCovid_dt, old = paste0("V", 1:53), new = week_cols)
PreCovid_dt <- PreCovid_dt %>% dplyr::select(-row_id)
PreCovid_dt[, c("week_start", "week_end") := NULL]

PreCovid_weekly <- PreCovid_dt %>%
  dplyr::select(-week_53) %>%
  dplyr::mutate(
    log_cons = log(total_cons / days + 1),
    year     = year(as.Date(previous_record))
  )

weekly_vars <- grep("^week_", names(PreCovid_weekly), value = TRUE)
reg_formula <- as.formula(paste("log_cons ~ 0 +", paste(weekly_vars, collapse = " + "),
                                "| power_grid_location + year"))
model_weekly <- feols(reg_formula, data = PreCovid_weekly, cluster = ~ id)

weekly_coefs <- broom::tidy(model_weekly, conf.int = TRUE) %>%
  dplyr::filter(grepl("^week_", term)) %>%
  dplyr::mutate(week = as.numeric(gsub("week_", "", term))) %>%
  dplyr::arrange(week)

save(weekly_coefs,  file = "Data/weekly_coefs.RData")

gam_fit <- mgcv::gam(estimate ~ s(week, bs = "cc", k = 50), data = weekly_coefs)
weekly_coefs <- dplyr::mutate(weekly_coefs,
                              smoothed = predict(gam_fit, newdata = weekly_coefs))

p <- ggplot(weekly_coefs, aes(x = week)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(aes(y = estimate, linetype = "Estimate"), size = 0.8) +
  geom_line(aes(y = smoothed,  linetype = "Smoothed"),  size = 0.8) +
  labs(title = "Estimated Weekly Effects on Log Consumption",
       x = "Week of Year", y = "Coefficient", linetype = NULL) +
  scale_linetype_manual(values = c("Estimate" = "solid", "Smoothed" = "dashed")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")
ggsave(file.path("Figures", "weekly_effects_commercial.pdf"),
       plot = p, width = 6, height = 5, units = "in", device = cairo_pdf)

## ---- Build weekly panel & CF ----
dt <- as.data.table(Data97to01)
dt[, previous_record := as.Date(previous_record)]
dt[, current_record  := as.Date(current_record)]
dt <- dt[, .(
  weeks_seq = seq(
    floor_date(previous_record, unit = "week", week_start = 1),
    floor_date(current_record,  unit = "week", week_start = 1),
    by = "1 week"
  )
), by = .(id, previous_record, current_record, power_grid_location, usage,
          total_cons, mid_cons, peak_cons, low_cons)]
dt[, week := isoweek(weeks_seq)]

df_weekly <- dt %>%
  filter(week != 53, year(current_record) < 2022, year(current_record) > 2017) %>%
  left_join(weekly_coefs %>% dplyr::select(week, estimate), by = "week") %>%
  mutate(weight = ifelse(is.na(estimate), 1, exp(estimate))) %>%
  group_by(id, previous_record, current_record) %>%
  mutate(
    weekly_total_cons = (weight / sum(weight, na.rm = TRUE)) * total_cons,
    weekly_mid_cons   = (weight / sum(weight, na.rm = TRUE)) * mid_cons,
    weekly_peak_cons  = (weight / sum(weight, na.rm = TRUE)) * peak_cons,
    weekly_low_cons   = (weight / sum(weight, na.rm = TRUE)) * low_cons
  ) %>% ungroup() %>% select(-estimate, -weight)

df_weekly <- df_weekly %>%
  group_by(id, usage, weeks_seq, week, power_grid_location) %>%
  summarise(
    weekly_total_cons = mean(weekly_total_cons, na.rm = TRUE),
    weekly_mid_cons   = mean(weekly_mid_cons,   na.rm = TRUE),
    weekly_peak_cons  = mean(weekly_peak_cons,  na.rm = TRUE),
    weekly_low_cons   = mean(weekly_low_cons,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    year      = year(weeks_seq),
    week      = isoweek(weeks_seq),
    prev_year = year - 1
  )

cf_lookup <- df_weekly %>%
  transmute(
    id, power_grid_location,
    cf_year = year(weeks_seq),
    cf_week = isoweek(weeks_seq),
    weekly_cf_total = weekly_total_cons,
    weekly_cf_mid   = weekly_mid_cons,
    weekly_cf_peak  = weekly_peak_cons,
    weekly_cf_low   = weekly_low_cons
  )

df_weekly <- df_weekly %>%
  left_join(cf_lookup,
            by = c("id","power_grid_location","prev_year"="cf_year","week"="cf_week"))

## ---- Trend dataset & plots + summary table ----
df_trend <- df_weekly %>%
  mutate(
    log_weekly_total_cons = log(weekly_total_cons + 1),
    log_weekly_mid_cons   = log(weekly_mid_cons   + 1),
    log_weekly_peak_cons  = log(weekly_peak_cons  + 1),
    log_weekly_low_cons   = log(weekly_low_cons   + 1)
  ) %>%
  group_by(weeks_seq) %>%
  summarise(
    mean_log_total_cons = mean(log_weekly_total_cons, na.rm = TRUE),
    mean_log_mid_cons   = mean(log_weekly_mid_cons,   na.rm = TRUE),
    mean_log_peak_cons  = mean(log_weekly_peak_cons,  na.rm = TRUE),
    mean_log_low_cons   = mean(log_weekly_low_cons,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!weeks_seq %in% as.Date(c(
    "2018-01-01","2018-01-08","2018-01-15",
    "2021-12-06","2021-12-13","2021-12-20","2021-12-27"
  ))) %>%
  mutate(year = year(weeks_seq), week = isoweek(weeks_seq)) %>%
  mutate(
    year_week      = paste0(year, "_", sprintf("%02d", week)),
    prev_year_week = paste0(year - 1, "_", sprintf("%02d", week))
  )

df_prev_year <- df_trend %>%
  select(
    year_week,
    mean_log_total_cons, mean_log_mid_cons, mean_log_peak_cons, mean_log_low_cons
  ) %>% rename_with(~ paste0("cf_", .x), -year_week)

df_trend <- df_trend %>%
  left_join(df_prev_year, by = c("prev_year_week" = "year_week")) %>%
  rename(
    counterfactual_total = cf_mean_log_total_cons,
    counterfactual_mid   = cf_mean_log_mid_cons,
    counterfactual_peak  = cf_mean_log_peak_cons,
    counterfactual_low   = cf_mean_log_low_cons
  )

save(df_trend,  file = "Data/df_trend.RData")
save(df_weekly, file = "Data/df_weekly.RData")

# Plot consumption trends
important_dates <- as.numeric(as.Date(c("2018-03-01", "2019-03-01", "2021-03-01")))
plot_consumption_with_counterfactual(df_trend, "mean_log_total_cons", "Log Weekly Consumption – Total", "consumption_total")
plot_consumption_with_counterfactual(df_trend, "mean_log_mid_cons",   "Log Weekly Consumption – Mid-Peak", "consumption_mid")
plot_consumption_with_counterfactual(df_trend, "mean_log_peak_cons",  "Log Weekly Consumption – Peak", "consumption_peak")
plot_consumption_with_counterfactual(df_trend, "mean_log_low_cons",   "Log Weekly Consumption – Off-Peak", "consumption_low")

# Summary statistics table
df_filtered <- df_weekly %>%
  filter(
    (weeks_seq >= as.Date("2019-01-11") & weeks_seq < as.Date("2019-06-21")) |
      (weeks_seq >= as.Date("2020-01-11") & weeks_seq < as.Date("2020-06-21"))
  )

df_selected <- df_filtered %>%
  mutate(
    year = year(weeks_seq),
    `Log Total Consumption`    = log1p(weekly_total_cons),
    `Log Mid-Peak Consumption` = log1p(weekly_mid_cons),
    `Log Peak Consumption`     = log1p(weekly_peak_cons),
    `Log Off-Peak Consumption` = log1p(weekly_low_cons),
    treated = as.integer(year == 2020)
  ) %>%
  select(`Log Total Consumption`,`Log Mid-Peak Consumption`,
         `Log Peak Consumption`,`Log Off-Peak Consumption`,
         treated, week, year) %>%
  drop_na()

var_order <- c("Log Total Consumption","Log Mid-Peak Consumption",
               "Log Peak Consumption","Log Off-Peak Consumption","treated","week","year")
var_labels <- c(
  "Log Total Consumption"     = "Log Total Consumption",
  "Log Mid-Peak Consumption"  = "Log Mid-Peak Consumption",
  "Log Peak Consumption"      = "Log Peak Consumption",
  "Log Off-Peak Consumption"  = "Log Off-Peak Consumption",
  treated                     = "Treated (2020 = 1)",
  week                        = "Week of Year",
  year                        = "Year"
)

tab_full <- summarize_vars(df_selected, var_order, var_labels)
tab_2019 <- summarize_vars(filter(df_selected, year == 2019), var_order, var_labels)
tab_2020 <- summarize_vars(filter(df_selected, year == 2020), var_order, var_labels)

ADD <- "__ADDLINE__"
panel_table <- dplyr::bind_rows(
  tibble::tibble(Statistics = "\\textbf{Panel A. Full sample}",
                 N = NA_real_, Mean = NA_real_, `Std. Dev.` = NA_real_,
                 Min = NA_real_, Median = NA_real_, Max = NA_real_),
  tab_full,
  tibble::tibble(Statistics = ADD),
  tibble::tibble(Statistics = "\\textbf{Panel B. 2019 sample}"),
  tab_2019,
  tibble::tibble(Statistics = ADD),
  tibble::tibble(Statistics = "\\textbf{Panel C. 2020 sample}"),
  tab_2020
)

panel_table_fmt <- panel_table %>%
  mutate(across(c(N, Mean, `Std. Dev.`, Min, Median, Max), fmt_num))

latex_tbl <- panel_table_fmt %>%
  knitr::kable(format = "latex", booktabs = TRUE,
               caption = "Summary Statistics for Weekly Electricity Consumption",
               label = "summary_stats",
               align = c("l", rep("r", 6)),
               escape = FALSE,
               col.names = c("Statistics","N","Mean","Std. Dev.","Min","Median","Max"),
               linesep = "") %>%
  kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

latex_str <- as.character(latex_tbl)

# keep your addlinespace replacements
latex_str <- gsub("(?m)^__ADDLINE__\\s*&.*?\\\\\\\\\\s*$", "\\\\addlinespace", latex_str, perl = TRUE)

# 🔧 make it full width with elastic spacing
latex_str <- sub(
  "\\\\begin\\{tabular\\}\\[t\\]\\{lrrrrrr\\}",
  "\\\\begin{tabular*}{\\\\textwidth}{@{\\\\extracolsep{\\\\fill}}lrrrrrr}",
  latex_str
)
latex_str <- sub("\\\\end\\{tabular\\}", "\\\\end{tabular*}", latex_str)

# 🔧 gentle column padding tweak
latex_str <- sub("\\\\centering", "\\\\centering\n\\\\setlength{\\\\tabcolsep}{6pt}", latex_str, fixed = TRUE)

writeLines(latex_str, "Tables/summary_stats.tex")

# 2019 vs 2020 aligned-week plots
plot_aligned_weeks_comparison(df_trend, "mean_log_total_cons", "Total", "total")
plot_aligned_weeks_comparison(df_trend, "mean_log_mid_cons",   "Mid-Peak", "mid")
plot_aligned_weeks_comparison(df_trend, "mean_log_peak_cons",  "Peak", "peak")
plot_aligned_weeks_comparison(df_trend, "mean_log_low_cons",   "Off-Peak", "low")