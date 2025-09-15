# --- Load inputs you already have
load("Data/weekly_coefs.RData")  # weekly_coefs with columns incl. week, estimate
# event_total: your feols event-study model

# --- Prepare global weights over ALL 52 weeks
weights_all <- weekly_coefs %>%
  mutate(global_weight = exp(estimate) / sum(exp(estimate))) %>%
  select(week, global_weight)

# Keep only weeks 10–25 and tag Jalali year (rough split)
wk_subset <- weights_all %>%
  filter(between(week, 10, 25)) %>%
  mutate(jyear = if_else(week <= 12, 1397L, 1398L))

# --- Extract event-study coefficients for periods 0..15 and map to weeks 10..25
coef_tbl <- broom::tidy(event_total) %>%
  filter(grepl("^period::[0-9]+:treated$", term)) %>%
  mutate(period = as.integer(gsub("period::([0-9]+):treated", "\\1", term))) %>%
  filter(period >= 0, period <= 15) %>%
  rename(ES_coefs = estimate) %>%     # treat ES_coefs as % loss in logs
  arrange(period) %>%
  mutate(week = 10 + period)          # period 0..15 -> weeks 10..25


# ===== Run all scenarios =====
compute_loss_for(33170175000, 35743662000,
                 "Province total consumption")
compute_loss_for(20565999000, 22007469000,
                 "City total consumption")
compute_loss_for(5936571000,  6453556000,
                 "Province commercial consumption")
compute_loss_for(1366476000,  1622563000,
                 "City commercial consumption")