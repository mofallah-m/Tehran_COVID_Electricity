# Project: COVID-19 Economic Activity via Electricity Consumption (Tehran)

This repo estimates commercial activity losses in Tehran during the first COVID-19 wave using high-frequency electricity consumption as a proxy for output. It builds weekly panels, runs DiD/event-study regressions, performs pre-trend/placebo checks, and converts estimated losses into KWh/GDP terms. Outputs (figures & LaTeX tables) are written to `./Figures` and `./Tables`.

## Repository structure

```
.
├── 00_master_run_all.R (the driver).R   # Orchestrates end-to-end run
├── 01_libraries_functions.R             # Packages + helper functions (plots, tables, printing)
├── 02_import_and_clean.R                # Load raw RData, build weekly panels, summary stats, baseline weights
├── 03_methods.R                         # Core DiD + event-study estimation & tables/figures
├── 04_tests.R                           # Pre-trend joint F test, placebo ES, alternative controls
├── 05_economics_cost.R                  # Converts ES losses to KWh (and supports GDP mapping)
├── Data/                                # Input & intermediate RData files (see below)
├── Figures/                             # Auto-generated PDFs (event-study plots, etc.)
└── Tables/                              # Auto-generated LaTeX tables
```

## Data inputs (expected in `./Data`)

- `Data94to98.RData` and `Data97to01.RData` — raw/combined consumption datasets.  
- Created during the pipeline:
  - `df_weekly.RData` — weekly, regression-ready panel (02 → used in 03).
  - `df_reg.RData` — aligned panel for DiD/ES and tests (03 → used in 04).
  - `weekly_coefs.RData` — smoothed and raw baseline seasonal weights from pre-COVID weeks (02 → used in 05).

## Outputs

- **Figures/**
  - `event_total.pdf`, `event_mid.pdf`, `event_peak.pdf`, `event_low.pdf` (event-study plots).
  - 2019 vs 2020 aligned-week comparison plots for total/mid/peak/off-peak.
- **Tables/**
  - `summary_stats.tex` (formatted with padded columns).
  - `eventstudy_coefs.tex` (baseline ES by relative week; ref = −1).
  - `eventstudy_coefs_altcontrol.tex` (ES with control = mean of 2018, 2019, 2021).

## How to run

1. **Set the working directory** inside `00_master_run_all.R (the driver).R`:
   ```r
   setwd("~/Desktop/TeIAS/Electricity")  # edit as needed
   ```
2. Ensure required input files exist in `./Data` (see list above).
3. Run the driver in a fresh R session:
   ```r
   rm(list=ls()); graphics.off()
   source("00_master_run_all.R (the driver).R")
   ```
4. Check `./Figures` and `./Tables` for outputs.

## Software & packages

- **R ≥ 4.1** recommended.
- CRAN packages used (loaded in `01_libraries_functions.R`):
  - `dplyr`, `tidyr`, `data.table`, `lubridate`, `jalcal`, `skimr`, `knitr`, `kableExtra`,
  - `ggplot2`, `patchwork`, `viridis`,  
  - `fixest`, `broom`, `modelsummary`, `stargazer`,  
  - `purrr`, `mgcv`, `car`.
- Install any missing packages, e.g.:
  ```r
  pkgs <- c("dplyr","tidyr","data.table","lubridate","jalcal","skimr","knitr","kableExtra",
            "ggplot2","patchwork","viridis","fixest","broom","modelsummary","stargazer",
            "purrr","mgcv","car")
  install.packages(setdiff(pkgs, rownames(installed.packages())))
  ```

## Pipeline overview

### 1) Libraries & helpers (`01_libraries_functions.R`)
- Loads all packages.
- Utility functions:
  - Number formatting (`fmt_num`).
  - Plot savers for event-study figures (single aesthetic, consistent layout).
  - LaTeX table helpers (modelsummary output tweaks).
  - Cost-conversion helpers (printing KWh/GDP loss summaries).

### 2) Import & cleaning (`02_import_and_clean.R`)
- Creates `Figures/`, `Tables/`, `Data/` if missing.
- Loads raw inputs (`Data94to98.RData`, `Data97to01.RData`) and row-binds.
- Builds **pre-COVID weekly fixed-effects baseline**:
  - Reconstructs weekly series; smooths seasonal pattern with `mgcv::gam(...)`.
  - Fits a pre-COVID weekly FE model, extracts **weekly coefficients**, saves as `Data/weekly_coefs.RData` and smooths them for weights.
- Produces **summary statistics** (written to `Tables/summary_stats.tex`, with widened `tabcolsep` and `tabular*` width control).
- Generates **2019 vs 2020 aligned-week plots** for total/mid/peak/off-peak loads.
- Saves weekly regression-ready data as `df_weekly.RData`.

### 3) Methods (DiD & Event Study) (`03_methods.R`)
- Loads `df_weekly.RData`.
- Computes **pre-period (2019)** moments by ID on a **log1p** scale:
  - `pre_mean`, `pre_var` (variance set to 0 if <2 obs), and splits IDs into terciles (e.g., by variance) for heterogeneity buckets.
- Builds `df_reg.RData` with aligned treated vs counterfactual series (and factors like `month`, `year`, etc.).
- Estimates:
  - Baseline **DiD** and **event-study** models using `fixest::feols` for:
    - **Total**, **Mid-Peak**, **Peak**, **Off-Peak** consumption.
  - Saves **event-study plots** (one PDF per load type).
  - Writes **event-study LaTeX table** (`Tables/eventstudy_coefs.tex`; ref period = −1).

### 4) Robustness & validation tests (`04_tests.R`)
- **Pre-trend joint F-test**: Tests that ES **leads** (periods −2 to −7) are jointly zero; uses `car::linearHypothesis`.
- **Placebo event study**: Restricts to 2019-12-01…2020-06-21, renames series (`total_treated`, `total_control`, etc.), and re-runs ES.
- **Alternative control**: Uses the mean of 2018, 2019, and 2021 as the counterfactual baseline; outputs `Tables/eventstudy_coefs_altcontrol.tex`.

### 5) Economic cost conversion (`05_economics_cost.R`)
- Loads `weekly_coefs.RData` (smoothed **seasonal weights**) and uses `event_total` ES results.
- Constructs **global weights** across 52 weeks, then restricts to **weeks 10–25** (the first-wave window).
- Maps **ES periods 0..15** to calendar weeks `10..25`, merges with weights, and **converts ES coefficients into KWh losses** by scaling with baseline 1397/1398 consumption levels.
- The helper `compute_loss_for()` prints total KWh losses for each scenario. The script demonstrates four examples:
  - Province total vs. city total consumption (2018/19 → 2019/20),
  - Province vs. city **commercial** consumption.

> If you additionally want **GDP** losses, provide the GDP-per-KWh benchmark (e.g., Tehran province 1397 GDP / 1397 commercial KWh) and multiply the estimated KWh loss by that ratio. You can embed that directly in `compute_loss_for()` or print alongside KWh.

## Key modeling choices (at a glance)

- **Outcome**: log or log1p of weekly electricity consumption (total/mid/peak/off-peak).
- **Design**: DiD with event-study around the COVID shock window (weeks 10–25).
- **Fixed effects**: `power_grid_location` and year/month FEs where applicable.
- **Clustering**: by customer/site (`id`).
- **Seasonality control**: pre-COVID weekly effects smoothed by `mgcv::gam` to build **weights**.
- **Validation**: joint F-test for pre-trends; placebo ES; alternative baseline control (mean of 2018, 2019, 2021).

## Common pitfalls & fixes

- **Mismatched filename**: Fix `05_economic_cost.R` vs `05_economics_cost.R`.
- **Missing packages**: Run the bulk installer above.
- **Working directory**: Update `setwd(...)` in the driver to your project root.
- **Missing `Data/*.RData`**: Ensure raw RData files exist; otherwise adapt `02_import_and_clean.R` to your raw sources.

## Reproducibility notes

- All outputs are programmatically generated; delete `Figures/` and `Tables/` to rebuild cleanly.
- Randomness is not used in the core pipeline; no seed required.
- The weekly weight construction uses smoothing (`mgcv::gam`) — keep package versions noted if you need byte-for-byte reproducibility.
