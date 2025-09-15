## ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(data.table); library(lubridate); library(jalcal)
  library(skimr); library(knitr); library(kableExtra)
  library(ggplot2); library(patchwork); library(viridis)
  library(fixest); library(broom); library(modelsummary); library(stargazer)
  library(purrr)
})

## ---- Utilities ----
fmt_num <- function(x) {
  ifelse(is.na(x), "",
         ifelse(abs(x - round(x)) < 1e-9,
                format(round(x), big.mark = ",", trim = TRUE),
                format(round(x, 2), nsmall = 2, big.mark = ",", trim = TRUE)))
}

summarize_vars <- function(d, vars, labels) {
  purrr::map_dfr(vars, function(v) {
    x <- d[[v]]
    tibble::tibble(
      Statistics  = if (!is.null(labels[[v]])) labels[[v]] else v,
      N           = sum(!is.na(x)),
      Mean        = mean(x, na.rm = TRUE),
      `Std. Dev.` = sd(x, na.rm = TRUE),
      Min         = suppressWarnings(min(x, na.rm = TRUE)),
      Median      = median(x, na.rm = TRUE),
      Max         = suppressWarnings(max(x, na.rm = TRUE))
    )
  })
}

## ---- Estimation helpers ----
run_did_model <- function(df, outcome) {
  fml <- as.formula(paste0(outcome, " ~ treated * post + week | id"))
  feols(fml, data = df, cluster = ~ id)
}

run_event_study <- function(df, outcome) {
  df <- df %>% mutate(outcome_var = .data[[outcome]])
  fml <- as.formula("outcome_var ~ i(period, treated) | week + id")
  feols(fml, cluster = ~ id, data = df)
}

run_placebo_event_study <- function(df, outcome, placebo_shift_weeks = 6) {
  df <- df %>%
    mutate(
      outcome_var    = .data[[outcome]],
      placebo_period = as.integer(
        difftime(weeks_seq, as.Date("2020-03-01") - weeks(placebo_shift_weeks), units = "weeks")
      )
    )
  fml <- as.formula("outcome_var ~ i(placebo_period, treated) | week + id")
  feols(fml, cluster = ~ id, data = df)
}

## ---- Plotting helpers ----
plot_centered_event_study <- function(model, ref_period, title = "Event Study") {
  est <- broom::tidy(model, conf.int = TRUE)
  est <- est[grepl("^period::(-?\\d+):treated$", est$term), ]
  est$period <- as.numeric(gsub("^period::(-?\\d+):treated$", "\\1", est$term))
  ref_row <- est[est$period == ref_period, ]
  if (nrow(ref_row) == 0) stop("Reference period not found in coefficients.")
  ref_est <- ref_row$estimate
  est <- est %>%
    mutate(
      estimate = estimate - ref_est,
      conf.low = conf.low - ref_est,
      conf.high = conf.high - ref_est,
      conf.low  = ifelse(period == ref_period, estimate, conf.low),
      conf.high = ifelse(period == ref_period, estimate, conf.high)
    )
  ggplot(est, aes(x = period, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "solid",  color = "black") +
    geom_point(color = "black", size = 1.75) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, color = "black") +
    ylim(-0.25, 0.07) +
    labs(x = "Weeks to Shock", y = expression(beta[k]), title = title) +
    theme_minimal(base_size = 13)
}

plot_placebo <- function(model, ref_period = -1, title = "Event Study", event_line_at = NULL) {
  est <- broom::tidy(model, conf.int = TRUE)
  term_match <- regmatches(est$term, regexpr("^[^:]+", est$term))
  interaction_var <- unique(term_match[grepl(":treated$", est$term)])
  if (length(interaction_var) != 1) stop("Couldn't identify interaction variable for plot.")
  pattern <- paste0("^", interaction_var, "::(-?\\d+):treated$")
  est <- est[grepl(pattern, est$term), ]
  est$period <- as.numeric(gsub(pattern, "\\1", est$term))
  ref_row <- est[est$period == ref_period, ]
  if (nrow(ref_row) == 0) stop("Reference period not found in coefficients.")
  ref_est <- ref_row$estimate
  est <- est %>%
    mutate(
      estimate = estimate - ref_est,
      conf.low = conf.low - ref_est,
      conf.high = conf.high - ref_est,
      conf.low  = ifelse(period == ref_period, estimate, conf.low),
      conf.high = ifelse(period == ref_period, estimate, conf.high)
    )
  p <- ggplot(est, aes(x = period, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "solid",  color = "black") +
    geom_point(color = "black", size = 1.75) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, color = "black") +
    ylim(-0.23, 0.07) +
    labs(x = "Weeks to (Placebo) Shock", y = expression(beta[k]), title = title) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (!is.null(event_line_at)) {
    p <- p + geom_vline(xintercept = event_line_at, linetype = "dashed", color = "black")
  }
  p
}

plot_aligned_weeks_comparison <- function(df, value_col, title_suffix = "", filename_suffix = "") {
  df <- df %>%
    mutate(
      group = case_when(
        weeks_seq >= as.Date("2019-01-11") & weeks_seq < as.Date("2019-07-04") ~ "2019",
        weeks_seq >= as.Date("2020-01-11") & weeks_seq < as.Date("2020-07-04") ~ "2020",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    arrange(weeks_seq, .by_group = TRUE) %>%
    mutate(week_aligned = row_number() - 7) %>%
    ungroup()
  df1 <- df %>% filter(group %in% c("2019", "2020"))
  p <- ggplot(df1, aes(x = week_aligned, y = .data[[value_col]], linetype = group)) +
    geom_line(color = "black", size = 0.8) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.6) +
    scale_linetype_manual(values = c("2019" = "dashed", "2020" = "solid")) +
    labs(
      title = paste("Electricity Consumption: 2019 vs 2020 —", title_suffix),
      x = "Weeks to Shock", y = "Log of Weekly Consumption (Mean)", linetype = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x  = element_text(hjust = 1),
          plot.title   = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom")
  dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
  filename <- file.path("Figures", paste0("2019Vs2020_", filename_suffix, ".pdf"))
  ggsave(filename, plot = p, width = 6, height = 4.5, units = "in", device = cairo_pdf)
  p
}

plot_consumption_with_counterfactual <- function(data, actual_col, title_text, filename_suffix = "") {
  p <- ggplot(data, aes(x = weeks_seq)) +
    geom_line(aes_string(y = actual_col), color = "black", size = 0.8, linetype = "solid") +
    geom_vline(xintercept = as.numeric(as.Date("2020-03-01")), linetype = "dashed", color = "black", size = 0.7) +
    geom_vline(xintercept = important_dates, linetype = "dotted", color = "black", alpha = 0.5, size = 0.6) +
    labs(title = title_text, x = "Date (Weeks)", y = "Log of Weekly Consumption (Mean)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold", hjust = 0.5))
  dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
  filename <- file.path("Figures", paste0(filename_suffix, ".pdf"))
  ggsave(filename, plot = p, width = 6, height = 4.5, units = "in", device = cairo_pdf)
  p
}

save_event_plot <- function(model, ref_period, title_text, file_name) {
  plot <- plot_centered_event_study(
    model      = model,
    ref_period = ref_period,
    title      = title_text
  ) +
    labs(title = title_text) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  dir.create(dirname(file_name), showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = file_name, plot = plot, width = 7, height = 5)
  plot
}

plot_event_study_by_location <- function(df, outcome_var, ref_period = -1) {
  location_order <- df %>%
    dplyr::group_by(power_grid_location) %>%
    dplyr::summarise(level_sum = sum(.data[[outcome_var]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(level_sum)) %>%
    dplyr::pull(power_grid_location)
  
  df <- df %>% dplyr::mutate(outcome = .data[[outcome_var]])
  
  df <- df %>%
    dplyr::group_split(power_grid_location) %>%
    purrr::map_dfr(function(g) {
      loc <- unique(g$power_grid_location)
      m <- feols(outcome ~ i(period, treated) | week + id + year, cluster = ~ id, data = g)
      td <- broom::tidy(m, conf.int = TRUE) %>%
        dplyr::filter(grepl("^period::(-?\\d+):treated$", term)) %>%
        dplyr::mutate(period = as.numeric(sub("^period::(-?\\d+):treated$", "\\1", term)))
      ref_row <- dplyr::filter(td, period == ref_period)
      if (nrow(ref_row) == 0) return(NULL)
      ref_est <- ref_row$estimate
      td %>%
        dplyr::mutate(
          estimate = estimate - ref_est,
          conf.low = conf.low - ref_est,
          conf.high = conf.high - ref_est,
          conf.low  = ifelse(period == ref_period, 0, conf.low),
          conf.high = ifelse(period == ref_period, 0, conf.high),
          power_grid_location = loc
        ) %>%
        dplyr::select(power_grid_location, period, estimate, conf.low, conf.high)
    })
  
  df <- df %>% dplyr::mutate(power_grid_location = factor(power_grid_location, levels = location_order))
  
  ggplot2::ggplot(df, ggplot2::aes(x = period, y = estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), fill = "black", alpha = 0.2) +
    ggplot2::geom_line(color = "black", linewidth = 0.8) +
    ggplot2::facet_wrap(~ power_grid_location, ncol = 4, scales = "free_y") +
    ggplot2::labs(x = "Weeks to Shock", y = expression(beta[k]),
                  title = paste("Event Study by Location -", toupper(outcome_var))) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   strip.text = ggplot2::element_text(size = 8))
}

# ---- 1) Helper: extract ES rows from a fixest model (robust to naming)
tidy_es <- function(mod, ref_period = -1) {
  tb <- broom::tidy(mod)
  
  # keep only interaction terms involving 'period' and 'treated'
  tb <- tb %>%
    filter(str_detect(term, "(period.*treated|treated.*period|i\\(period.*\\))"))
  
  # pull out the period integer from the term name
  tb <- tb %>%
    mutate(
      rel = as.integer(str_replace(term, ".*?(\\-?\\d+).*", "\\1"))
    ) %>%
    filter(!is.na(rel), rel != ref_period) %>%  # drop reference period
    transmute(rel_week = rel,
              estimate, std.error, p.value)
  tb
}

# ---- 2) Helper: stars and pretty "est (se)" cell
stars_from_p <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01 ,  "**",
                ifelse(p < 0.05 ,   "*",
                       ifelse(p < 0.10 ,   "+", ""))))
}
cell_fmt <- function(est, se, p) {
  paste0(sprintf("%.3f", est), stars_from_p(p), " (", sprintf("%.3f", se), ")")
}

# ---- 3) Build a single ES coefficient table (wide, full-width)
tidy_es <- function(mod, ref_period = -1) {
  tb <- broom::tidy(mod)
  
  # keep terms that clearly belong to the ES (involving 'period' & 'treated' or fixest i() )
  tb <- tb %>%
    filter(str_detect(term, "period") & str_detect(term, "treated") | str_detect(term, "^i\\(period")) %>%
    mutate(
      rel_week = suppressWarnings(as.integer(str_replace(term, ".*?(-?\\d+).*", "\\1")))
    ) %>%
    filter(!is.na(rel_week), rel_week != ref_period) %>%
    transmute(rel_week, estimate, std.error, p.value)
  tb
}

# --- Helper: compute total KWh loss for a scenario
compute_loss_for <- function(annual97, annual98, label) {
  annual_map <- c(`1397` = annual97, `1398` = annual98)
  
  weekly_baseline <- wk_subset %>%
    mutate(weekly_kwh = global_weight * annual_map[as.character(jyear)]) %>%
    select(week, weekly_kwh)
  
  loss_tbl <- coef_tbl %>%
    left_join(weekly_baseline, by = "week") %>%
    mutate(loss_kwh = weekly_kwh * (-ES_coefs))   # ES_coefs as % loss
  
  total_loss <- sum(loss_tbl$loss_kwh, na.rm = TRUE)
  
  cat(
    label, "is",
    format(round(total_loss, 0), big.mark = ","),
    "KWh.\n"
  )
  invisible(total_loss)
}

# --- Helper

stars_from_p <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01 ,  "**",
                ifelse(p < 0.05 ,   "*",
                       ifelse(p < 0.10 ,   "+", ""))))
}
cell_fmt <- function(est, se, p) {
  paste0(sprintf("%.3f", est), stars_from_p(p), " (", sprintf("%.3f", se), ")")
}

# --- main builder: writes a full-width LaTeX table with no notes -----------
make_es_table <- function(models_named, caption, label,
                          outfile = "Tables/eventstudy_coefs.tex",
                          ref_period = -1,
                          window = NULL # e.g., c(-6, 15) if you want to trim
) {
  stopifnot(length(models_named) >= 1)
  
  # Tidy each model & tag with its outcome name
  es_list <- lapply(names(models_named), function(nm) {
    tidy_es(models_named[[nm]], ref_period = ref_period) %>%
      dplyr::mutate(outcome = nm)
  })
  es_all <- dplyr::bind_rows(es_list)
  
  # Optional: restrict horizon
  if (!is.null(window) && length(window) == 2) {
    es_all <- es_all %>% filter(rel_week >= window[1], rel_week <= window[2])
  }
  
  # wide table: rows = rel_week, columns = outcomes, values = "est (se)"
  es_wide <- es_all %>%
    mutate(cell = cell_fmt(estimate, std.error, p.value)) %>%
    select(rel_week, outcome, cell) %>%
    tidyr::pivot_wider(names_from = outcome, values_from = cell) %>%
    arrange(rel_week)
  
  es_wide[is.na(es_wide)] <- ""
  
  # ----- compose LaTeX manually (tabular* to fill \textwidth) -------------
  ncols <- ncol(es_wide)
  stopifnot(ncols >= 2)
  
  colspec <- paste0("@{\\extracolsep{\\fill}}", paste0(c("r", rep("c", ncols - 1)), collapse = ""))
  
  header_line <- paste(c("Rel. Week", colnames(es_wide)[-1]), collapse = " & ")
  
  body_lines <- apply(es_wide, 1, function(row) paste(row, collapse = " & "))
  
  lines <- c(
  "\\begin{table}[H]",
  "\\centering",
  sprintf("\\caption{%s}", caption),   # caption appears first
  sprintf("\\label{%s}", label),       # label after caption
  sprintf("\\begin{tabular*}{\\textwidth}{%s}", colspec),
  "\\toprule",
  header_line, "\\\\",
  "\\midrule",
  paste0(body_lines, collapse = " \\\\\n"),
  "\\bottomrule",
  "\\midrule",
  sprintf("\\multicolumn{%d}{l}{\\rule{0pt}{1em}+ p $<$ 0.1, * p $<$ 0.05, ** p $<$ 0.01, *** p $<$ 0.001}\\\\", ncols),
  "\\end{tabular*}",
  "\\end{table}"
)
  
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, con = outfile)
}