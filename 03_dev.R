# ============================================================
# dev.R — Main and hospital-specific model development
# ============================================================
# Run 01_explore.R and 02_breakpoint.R first, from the project root.
# Uses their cleaned `data` and saved breakpoint specification.
# All model selection uses development data only. Zero counts are retained.
# The workflow follows dev.R, dev_nd2.R and dev_bnd.R:
# benchmarks -> DR grid -> weather screening -> 26 combinations -> freeze.
# Corrected selection may differ from previously reported model rankings.

library(fpp3)
library(purrr)
library(splines)
library(prophet)
if (requireNamespace("ggtime", quietly = TRUE)) library(ggtime)

# ------------------------------------------------------------
# 1. Inputs and shared origin-week predictors
# ------------------------------------------------------------

dev_output_dir <- file.path("results", "development")
dir.create(dev_output_dir, recursive = TRUE, showWarnings = FALSE)
breakpoint_spec <- readRDS(file.path("results", "breakpoints", "breakpoint_spec.rds"))

# Count distinct holiday dates, rather than holiday names.
holiday_weekly <- generated_holidays |>
  filter(country == "VN", year >= 2015, year <= 2024) |>
  mutate(ds = as.Date(ds),
         week_start = floor_date(ds, "week", week_start = 1)) |>
  distinct(ds, .keep_all = TRUE) |>
  count(week_start, name = "n_holiday_days")

weekly_data <- data |>
  as_tibble() |>
  select(week_start, week_id, phase, y_lri, y_hosp1, y_hosp2,
         temp_mean, rh_mean, rain_sum) |>
  arrange(week_start) |>
  left_join(holiday_weekly, by = "week_start") |>
  mutate(
    n_holiday_days = replace_na(n_holiday_days, 0L),
    T0 = temp_mean, T1 = lag(temp_mean),
    T01 = (T0 + T1) / 2, T02 = (T0 + T1 + lag(temp_mean, 2)) / 3,
    H0 = rh_mean, H1 = lag(rh_mean),
    H01 = (H0 + H1) / 2, H02 = (H0 + H1 + lag(rh_mean, 2)) / 3,
    R0 = rain_sum, R1 = lag(rain_sum),
    R01 = R0 + R1, R02 = R0 + R1 + lag(rain_sum, 2)
  )

weather_domains <- list(
  temperature = c("T0", "T1", "T01", "T02"),
  humidity = c("H0", "H1", "H01", "H02"),
  rainfall = c("R0", "R1", "R01", "R02")
)
pdq_grid <- tibble(p = c(0, 1, 2, 0, 1, 2), d = 0,
                   q = c(0, 0, 0, 1, 1, 1))
dr_grid <- crossing(df = 1:4, K = 1:6, holiday = c(FALSE, TRUE), pdq_grid)

# ------------------------------------------------------------
# 2. Small helpers used by all three analyses
# ------------------------------------------------------------

fit_dr <- function(sample, lambda, terms, K, p, d, q) {
  rhs <- c(terms, sprintf("fourier(period = 'year', K = %s)", K),
           sprintf("pdq(%s,%s,%s)", p, d, q), "PDQ(0,0,0)", "1")
  formula <- as.formula(paste("box_cox(y_target + 1, lambda) ~",
                             paste(rhs, collapse = " + ")))
  sample |> model(DR = ARIMA(formula))
}

rank_models <- function(candidates) {
  ranked <- candidates |>
    mutate(info = map(fit, glance)) |>
    unnest(info) |>
    arrange(AICc)
  # A failed fit must not silently disappear from the prespecified search.
  if (any(!is.finite(ranked$AICc))) stop("A candidate failed; inspect model warnings.")
  ranked |> mutate(rank = row_number(), delta_AICc = AICc - min(AICc))
}

freeze_spline <- function(x, df, prefix) {
  basis <- ns(x, df = df)
  list(knots = attr(basis, "knots"),
       boundary = attr(basis, "Boundary.knots"),
       columns = paste0(prefix, seq_len(ncol(basis))))
}

add_spline <- function(tbl, variable, spec) {
  basis <- ns(tbl[[variable]], knots = spec$knots,
              Boundary.knots = spec$boundary)
  tbl[spec$columns] <- as.data.frame(basis)
  tbl
}

ljung_box_result <- function(fit, dof) {
  fit |> augment() |>
    features(.innov, ljung_box, lag = 24, dof = dof) |>
    as_tibble() |> select(lb_stat, lb_pvalue)
}

# ------------------------------------------------------------
# 3. Develop one outcome; repeat independently for each hospital
# ------------------------------------------------------------

develop_series <- function(series_id, input_spec) {
  lambda <- input_spec$lambda
  breaks <- sort(input_spec$regime_start_dates)

  # A breakpoint is the first week of its new regime.
  origin_data <- weekly_data |>
    mutate(
      y = .data[[input_spec$outcome_column]],
      regime = factor(findInterval(as.numeric(week_start), as.numeric(breaks)) + 1L,
                      levels = seq_len(length(breaks) + 1L)),
      target_week = lead(week_start), y_target = lead(y),
      target_week_id = lead(week_id), target_regime = lead(regime),
      target_holiday = lead(n_holiday_days), target_phase = lead(phase)
    )

  # Benchmarks retain all development observations, as in the original scripts.
  benchmark_data <- origin_data |>
    filter(phase == "development") |> as_tsibble(index = week_start)
  benchmarks <- benchmark_data |> model(
    SNAIVE = SNAIVE(box_cox(y + 1, lambda) ~ lag(52)),
    ETS = ETS(box_cox(y + 1, lambda), ic = "aicc")
  )

  # Index DR models by the target week; weather remains at the origin or earlier.
  dev <- origin_data |>
    filter(target_phase == "development", !is.na(y_target)) |>
    as_tsibble(index = target_week)

  # 3A. Joint selection: trend df x Fourier K x holiday x ARMA order (288).
  baseline_rank <- dr_grid |>
    mutate(fit = pmap(list(df, K, holiday, p, d, q),
      function(df, K, holiday, p, d, q) {
        terms <- c(sprintf("ns(target_week_id, df = %s)", df), "target_regime",
                   if (holiday) "target_holiday")
        fit_dr(dev, lambda, terms, K, p, d, q)
      })) |>
    rank_models()
  best <- baseline_rank |> slice(1)

  # Retain the original top-three baseline residual checks.
  baseline_diagnostics <- map_dfr(1:3, function(i) {
    ljung_box_result(baseline_rank$fit[[i]],
                     baseline_rank$p[i] + baseline_rank$q[i]) |>
      mutate(rank = i, .before = 1)
  })

  # Freeze the selected time basis before reducing the sample for weather.
  trend_spec <- freeze_spline(dev$target_week_id, best$df, "trend_")
  dev <- add_spline(dev, "target_week_id", trend_spec)
  backbone <- c(trend_spec$columns, "target_regime",
                if (best$holiday) "target_holiday")
  fit_selected <- function(sample, weather_terms = character()) {
    fit_dr(sample, lambda, c(backbone, weather_terms),
           best$K, best$p, best$d, best$q)
  }

  # 3B. Eight candidates per domain, each evaluated on the same domain sample.
  # Store weather bases as well, so combinations and validation reuse them.
  weather_splines <- list()
  screening <- list()
  for (domain in names(weather_domains)) {
    variables <- weather_domains[[domain]]
    domain_data <- dev |> filter(if_all(all_of(variables), ~ !is.na(.x)))
    for (variable in variables) {
      spline_spec <- freeze_spline(domain_data[[variable]], 3, paste0(variable, "_s"))
      weather_splines[[variable]] <- spline_spec
      dev <- add_spline(dev, variable, spline_spec)
      domain_data <- add_spline(domain_data, variable, spline_spec)
    }
    screening[[domain]] <- crossing(variable = variables, form = c("linear", "ns3")) |>
      mutate(
        term = if_else(form == "linear", variable, paste0("ns(", variable, ", df = 3)")),
        fit = map2(variable, form, function(variable, form) {
          terms <- if (form == "linear") variable else weather_splines[[variable]]$columns
          fit_selected(domain_data, terms)
        })
      ) |> rank_models() |> mutate(domain = domain, .before = 1)
  }
  screening_rank <- bind_rows(screening)
  retained <- map_dfr(screening, ~ slice_head(.x, n = 2))

  # 3C. Two retained terms per domain: 3 x 3 x 3 - 1 = 26 combinations.
  # NA denotes omission of that domain; the empty model is excluded here.
  weather_grid <- crossing(
    temperature = c(NA_character_, retained$term[retained$domain == "temperature"]),
    humidity = c(NA_character_, retained$term[retained$domain == "humidity"]),
    rainfall = c(NA_character_, retained$term[retained$domain == "rainfall"])
  ) |> filter(if_any(everything(), ~ !is.na(.x)))
  common_data <- dev |>
    filter(if_all(all_of(unique(retained$variable)), ~ !is.na(.x)))

  # Convert readable formula labels to the stored spline columns.
  expand_terms <- function(terms) {
    terms <- terms[!is.na(terms)]
    unlist(lapply(terms, function(term) {
      row <- retained[match(term, retained$term), ]
      if (row$form == "linear") row$variable else weather_splines[[row$variable]]$columns
    }), use.names = FALSE)
  }
  weather_rank <- weather_grid |>
    mutate(fit = pmap(list(temperature, humidity, rainfall),
      function(temperature, humidity, rainfall) {
        fit_selected(common_data, expand_terms(c(temperature, humidity, rainfall)))
      })) |> rank_models()
  finalists <- weather_rank |> slice_head(n = 3)

  # 3D. Refit DR_BASE on exactly the same sample as W1-W3 for Delta AICc.
  final_fits <- c(list(DR_BASE = fit_selected(common_data)),
                 setNames(finalists$fit, paste0("W", 1:3)))
  final_comparison <- imap_dfr(final_fits, function(fit, model_name) {
    glance(fit) |> as_tibble() |> select(AIC, AICc, BIC, sigma2) |>
      bind_cols(ljung_box_result(fit, best$p + best$q)) |>
      mutate(model = model_name, n = nrow(common_data), .before = 1)
  }) |>
    mutate(delta_AICc_vs_DR = AICc - AICc[model == "DR_BASE"],
           residual_sd = sqrt(sigma2),
           sd_reduction_pct = 100 * (1 - residual_sd / residual_sd[model == "DR_BASE"]))

  # Save the selected definitions, not the previously reported winners.
  weather_labels <- setNames(pmap(finalists[c("temperature", "humidity", "rainfall")],
    function(temperature, humidity, rainfall) {
      terms <- c(temperature, humidity, rainfall)
      terms[!is.na(terms)]
    }), paste0("W", 1:3))
  specification <- list(
    series = series_id, outcome_column = input_spec$outcome_column,
    lambda = lambda, offset = 1, horizon = 1,
    development_end = max(benchmark_data$week_start),
    regime_start_dates = breaks, breakpoint_convention = "first_week_new_regime",
    seasonal_lag = 52,
    ets = fabletools::model_sum(benchmarks$ETS[[1]]),
    dr = best |> select(df, K, holiday, p, d, q),
    trend = trend_spec, weather_splines = weather_splines,
    backbone_terms = backbone,
    weather_labels = weather_labels,
    weather_terms = lapply(weather_labels, expand_terms),
    n_baseline = nrow(dev), n_weather = nrow(common_data),
    common_target_weeks = common_data$target_week,
    primary_comparator = "DR_BASE", primary_weather = "W1"
  )

  # Export compact tables; model objects stay available in the R session.
  series_dir <- file.path(dev_output_dir, series_id)
  dir.create(series_dir, showWarnings = FALSE)
  tables <- list(
    baseline_grid = baseline_rank,
    baseline_diagnostics = baseline_diagnostics,
    weather_screening = screening_rank,
    retained_weather = retained,
    weather_combinations = weather_rank,
    final_comparison = final_comparison,
    benchmark_summary = glance(benchmarks),
    # Preserve the benchmark convention used in the existing reporting code.
    # Report these settings explicitly; DR diagnostics use lag 24 and p + q.
    benchmark_diagnostics = benchmarks |> augment() |>
      features(.innov, ljung_box, lag = 52, dof = 0) |>
      mutate(lag = 52L, dof = 0L),
    benchmark_coefficients = tidy(benchmarks),
    final_coefficients = imap_dfr(final_fits, ~ tidy(.x) |> mutate(model = .y)),
    finalist_definitions = tibble(model = names(weather_labels),
                                 weather = map_chr(weather_labels, ~ paste(.x, collapse = " + ")))
  )
  for (table_name in names(tables)) {
    write.csv(tables[[table_name]] |> as_tibble() |> select(where(~ !is.list(.x))),
              file.path(series_dir, paste0(table_name, ".csv")), row.names = FALSE)
  }
  # Residual diagnostics are part of the manuscript, rather than input QC.
  for (model_name in names(final_fits)) {
    ggsave(file.path(series_dir, paste0(model_name, "_residuals.pdf")),
           gg_tsresiduals(final_fits[[model_name]], lag_max = 104), width = 9, height = 7)
  }
  for (model_name in c("SNAIVE", "ETS")) {
    ggsave(file.path(series_dir, paste0(model_name, "_residuals.pdf")),
           gg_tsresiduals(select(benchmarks, all_of(model_name)), lag_max = 104),
           width = 9, height = 7)
  }
  print(final_comparison)
  list(specification = specification, tables = tables,
       benchmarks = benchmarks, final_fits = final_fits)
}

# ------------------------------------------------------------
# 4. Run the three development analyses and freeze specifications
# ------------------------------------------------------------

development_results <- list()
for (series_id in c("MAIN", "ND2", "BND")) {
  message("Developing models: ", series_id)
  development_results[[series_id]] <- develop_series(series_id, breakpoint_spec[[series_id]])
}
development_spec <- lapply(development_results, function(x) x$specification)
saveRDS(development_spec, file.path(dev_output_dir, "development_spec.rds"))
write.csv(holiday_weekly, file.path(dev_output_dir, "holiday_calendar.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(dev_output_dir, "sessionInfo.txt"))

# Validation must read development_spec.rds and reuse its spline bases and terms.
# Refit coefficients at each origin, but do not repeat model selection.
# Update manuscript tables from these outputs if corrected results differ.
