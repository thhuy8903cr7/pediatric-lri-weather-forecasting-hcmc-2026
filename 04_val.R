# ============================================================
# val.R — Rolling validation for MAIN, ND2 and BND
# ============================================================
# Run 01_explore.R first to create cleaned `data` in this R session.
# Run 03_dev.R once beforehand to save results/development outputs.
# Working directory: project root, as for the preceding scripts.
# Validation: 2023-2024; expanding history; weekly refit; horizon = 1.
# All structures and spline bases are read from development_spec.rds.

library(fpp3)
library(purrr)
library(splines)
library(distillery)

# ------------------------------------------------------------
# 1. Inputs and shared origin-week weather
# ------------------------------------------------------------

val_output_dir <- file.path("results", "validation")
dir.create(val_output_dir, recursive = TRUE, showWarnings = FALSE)
dev_dir <- file.path("results", "development")
validation_spec <- readRDS(file.path(dev_dir, "development_spec.rds"))
holiday_weekly <- read.csv(file.path(dev_dir, "holiday_calendar.csv")) |>
  mutate(week_start = as.Date(week_start))
validation_years <- 2023:2024
bootstrap_B <- 2000L
bootstrap_block <- 11L
bootstrap_seed <- 20260827L

# Use the same holiday calendar and exposure definitions as development.
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

# ------------------------------------------------------------
# 2. Helpers: frozen splines and paired circular block bootstrap
# ------------------------------------------------------------

apply_spline <- function(tbl, variable, spline_spec) {
  basis <- ns(tbl[[variable]], knots = spline_spec$knots,
              Boundary.knots = spline_spec$boundary)
  tbl[spline_spec$columns] <- as.data.frame(basis)
  tbl
}

make_cbb_indices <- function(n, B, block_length) {
  replicate(B, {
    starts <- sample(seq_len(n), ceiling(n / block_length), replace = TRUE)
    indices <- unlist(lapply(starts, function(s) {
      ((s - 1L + 0:(block_length - 1L)) %% n) + 1L
    }))
    indices[seq_len(n)]
  })
}

bootstrap_delta <- function(predictions, model_name) {
  paired <- predictions |>
    filter(.model %in% c("DR_BASE", model_name)) |>
    select(target_week, .model, error) |>
    pivot_wider(names_from = .model, values_from = error) |>
    arrange(target_week)
  errors <- as.matrix(select(paired, DR_BASE, all_of(model_name)))
  delta_rmse <- function(data, ...) {
    sqrt(mean(data[, 2]^2)) - sqrt(mean(data[, 1]^2))
  }
  # Identical seed and paired week indices for each weather comparison.
  set.seed(bootstrap_seed)
  indices <- make_cbb_indices(nrow(errors), bootstrap_B, bootstrap_block)
  boot <- booter(x = errors, statistic = delta_rmse, shuffle = indices)
  ci <- quantile(boot$results, c(0.025, 0.975))
  tibble(.model = model_name, n = nrow(errors),
         B = bootstrap_B, block_length = bootstrap_block, seed = bootstrap_seed,
         Delta_RMSE = unname(boot$original.est),
         CI_low = unname(ci[1]), CI_high = unname(ci[2]))
}

# ------------------------------------------------------------
# 3. Validate one series using its own frozen development choices
# ------------------------------------------------------------

validate_series <- function(series_id, spec) {
  lambda <- spec$lambda
  offset <- spec$offset
  seasonal_lag <- spec$seasonal_lag
  breaks <- sort(spec$regime_start_dates)

  aligned <- weekly_data |>
    mutate(
      y = .data[[spec$outcome_column]],
      regime = factor(findInterval(as.numeric(week_start), as.numeric(breaks)) + 1L,
                      levels = seq_len(length(breaks) + 1L)),
      target_week = lead(week_start), y_target = lead(y),
      target_week_id = lead(week_id), target_regime = lead(regime),
      target_holiday = lead(n_holiday_days)
    ) |>
    filter(!is.na(target_week))

  # Evaluate the saved bases; never estimate new knots during validation.
  aligned <- apply_spline(aligned, "target_week_id", spec$trend)
  for (variable in names(spec$weather_splines)) {
    aligned <- apply_spline(aligned, variable, spec$weather_splines[[variable]])
  }
  targets <- aligned |>
    filter(year(target_week) %in% validation_years) |>
    transmute(origin_week = week_start, target_week, observed = y_target)
  first_target <- min(targets$target_week)
  last_target <- max(targets$target_week)

  # Preserve all selected components, including holiday inclusion and ARMA order.
  make_dr_formula <- function(extra = character()) {
    rhs <- c(spec$backbone_terms, extra,
             sprintf("fourier(period = 'year', K = %s)", spec$dr$K),
             sprintf("pdq(%s,%s,%s)", spec$dr$p, spec$dr$d, spec$dr$q),
             "PDQ(0,0,0)", "1")
    as.formula(paste("box_cox(y_target + offset, lambda) ~",
                     paste(rhs, collapse = " + ")))
  }
  form_base <- make_dr_formula()
  form_w1 <- make_dr_formula(spec$weather_terms$W1)
  form_w2 <- make_dr_formula(spec$weather_terms$W2)
  form_w3 <- make_dr_formula(spec$weather_terms$W3)

  # Convert the saved ETS label to an explicitly fixed formula.
  # For example ETS(A,N,N) becomes error('A') + trend('N') + season('N').
  ets_parts <- strsplit(gsub("^ETS\\(|\\)$| ", "", spec$ets), ",")[[1]]
  form_ets <- as.formula(sprintf(
    "box_cox(y_target + offset, lambda) ~ error('%s') + trend('%s') + season('%s')",
    ets_parts[1], ets_parts[2], ets_parts[3]
  ))
  required <- unique(c(spec$backbone_terms, unlist(spec$weather_terms)))

  # Preserve the original validation history: benchmarks start at the first
  # available direct target (week 2); DR starts on the common weather sample.
  bench_source <- aligned |>
    filter(target_week < last_target) |>
    select(target_week, y_target) |> as_tsibble(index = target_week)
  dr_source <- aligned |>
    filter(target_week < last_target,
           target_week %in% spec$common_target_weeks | target_week > spec$development_end,
           if_all(all_of(required), ~ !is.na(.x))) |>
    select(target_week, y_target, all_of(required)) |> as_tsibble(index = target_week)

  cv_bench <- bench_source |>
    stretch_tsibble(.init = sum(bench_source$target_week < first_target), .step = 1)
  cv_dr <- dr_source |>
    stretch_tsibble(.init = sum(dr_source$target_week < first_target), .step = 1)
  future_bench <- new_data(cv_bench, n = 1)
  lookup <- aligned |> select(target_week, all_of(required))
  future_dr <- new_data(cv_dr, n = 1) |>
    left_join(lookup, by = "target_week")

  # Essential alignment check: every fit ends one week before its target.
  history <- bind_rows(
    as_tibble(cv_bench) |> mutate(family = "benchmark"),
    as_tibble(cv_dr) |> mutate(family = "DR")
  ) |>
    group_by(family, .id) |>
    summarise(training_start = min(target_week), training_end = max(target_week),
              n_train = n(), .groups = "drop") |>
    mutate(target_week = training_end + 7)
  stopifnot(nrow(targets) == 104L,
            all(targets$target_week == targets$origin_week + 7),
            setequal(future_bench$target_week, targets$target_week),
            setequal(future_dr$target_week, targets$target_week),
            all(history$target_week %in% targets$target_week),
            !anyNA(select(as_tibble(future_dr), all_of(required))))

  # Weekly parameter estimation; no reselection of model structures.
  fit_bench <- cv_bench |> model(
    SNAIVE = SNAIVE(box_cox(y_target + offset, lambda) ~ lag(seasonal_lag)),
    ETS = ETS(form_ets)
  )
  fit_dr <- cv_dr |> model(
    DR_BASE = ARIMA(form_base), W1 = ARIMA(form_w1),
    W2 = ARIMA(form_w2), W3 = ARIMA(form_w3)
  )
  forecasts <- bind_rows(
    forecast(fit_bench, new_data = future_bench) |> as_tibble(),
    forecast(fit_dr, new_data = future_dr) |> as_tibble()
  )

  # Fable returns means and distributions on the original count scale.
  # Do not invert Box-Cox again, subtract the offset again, or round forecasts.
  predictions <- forecasts |>
    left_join(targets, by = "target_week") |>
    transmute(.id, .model, origin_week, target_week, observed,
              predicted = .mean, error = observed - .mean,
              median = quantile(y_target, 0.50),
              lo80 = quantile(y_target, 0.10), hi80 = quantile(y_target, 0.90),
              lo95 = quantile(y_target, 0.025), hi95 = quantile(y_target, 0.975)) |>
    arrange(target_week, .model)
  counts <- predictions |> count(.model, target_week)
  stopifnot(nrow(counts) == 6L * nrow(targets), all(counts$n == 1L),
            all(is.finite(as.matrix(select(predictions, observed, predicted,
                                          median, lo80, hi80, lo95, hi95)))))

  # Direct aggregation equals sqrt(mean(fold_RMSE^2)) for one-step folds.
  metrics <- predictions |> group_by(.model) |>
    summarise(n = n(), RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)),
              .groups = "drop") |>
    mutate(Delta_RMSE = RMSE - RMSE[.model == "DR_BASE"],
           RMSE_improvement_pct = -100 * Delta_RMSE / RMSE[.model == "DR_BASE"])

  # The original MAIN analysis includes interval scores and paired bootstrap.
  # Hospital sensitivity analyses retain point metrics, as in the original files.
  bootstrap <- NULL
  if (series_id == "MAIN") {
    predictions <- predictions |> mutate(
      IS80 = hi80 - lo80 + 10 * pmax(lo80 - observed, 0) + 10 * pmax(observed - hi80, 0),
      IS95 = hi95 - lo95 + 40 * pmax(lo95 - observed, 0) + 40 * pmax(observed - hi95, 0),
      WIS = (0.5 * abs(observed - median) + 0.10 * IS80 + 0.025 * IS95) / 2.5
    )
    interval_metrics <- predictions |> group_by(.model) |>
      summarise(WIS = mean(WIS),
                Coverage80 = mean(observed >= lo80 & observed <= hi80),
                Coverage95 = mean(observed >= lo95 & observed <= hi95), .groups = "drop")
    bootstrap <- map_dfr(c("W1", "W2", "W3"), ~ bootstrap_delta(predictions, .x))
    metrics <- metrics |>
      left_join(interval_metrics, by = ".model") |>
      left_join(select(bootstrap, .model, CI_low, CI_high), by = ".model") |>
      mutate(Relative_CI_low_pct = 100 * CI_low / RMSE[.model == "DR_BASE"],
             Relative_CI_high_pct = 100 * CI_high / RMSE[.model == "DR_BASE"])
  }

  series_dir <- file.path(val_output_dir, series_id)
  dir.create(series_dir, showWarnings = FALSE)
  write.csv(predictions, file.path(series_dir, "forecasts.csv"), row.names = FALSE)
  write.csv(metrics, file.path(series_dir, "metrics.csv"), row.names = FALSE)
  write.csv(history, file.path(series_dir, "training_windows.csv"), row.names = FALSE)
  if (!is.null(bootstrap)) {
    write.csv(bootstrap, file.path(series_dir, "bootstrap.csv"), row.names = FALSE)
  }
  # Export the exact formulas alongside scores to identify each W label.
  formula_table <- tibble(
    model = c("DR_BASE", "W1", "W2", "W3", "ETS"),
    formula = map_chr(list(form_base, form_w1, form_w2, form_w3, form_ets),
                      ~ paste(deparse(.x), collapse = " ")),
    lambda = lambda, offset = offset
  )
  write.csv(formula_table, file.path(series_dir, "model_formulas.csv"), row.names = FALSE)
  print(metrics)
  list(metrics = metrics, predictions = predictions, bootstrap = bootstrap)
}

# ------------------------------------------------------------
# 4. Run MAIN and the two hospital sensitivity analyses
# ------------------------------------------------------------

validation_results <- list()
for (series_id in c("MAIN", "ND2", "BND")) {
  message("Rolling validation: ", series_id)
  validation_results[[series_id]] <- validate_series(series_id, validation_spec[[series_id]])
}
validation_metrics <- imap_dfr(validation_results, ~ mutate(.x$metrics, series = .y, .before = 1))
write.csv(validation_metrics, file.path(val_output_dir, "validation_metrics.csv"), row.names = FALSE)
saveRDS(validation_spec, file.path(val_output_dir, "specification_used.rds"))
capture.output(sessionInfo(), file = file.path(val_output_dir, "sessionInfo.txt"))

# ------------------------------------------------------------
# 5. Figure 2: MAIN observed counts, forecast means and 95% intervals
# ------------------------------------------------------------

plot_data <- validation_results$MAIN$predictions |>
  filter(.model %in% c("ETS", "DR_BASE", "W1")) |>
  mutate(.model = factor(.model, levels = c("ETS", "DR_BASE", "W1")))
model_colours <- c(Observed = "black", ETS = "#A97800", DR_BASE = "#C65A4B", W1 = "#2C9C95")
figure_2 <- ggplot(plot_data, aes(target_week)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = .model), alpha = 0.22, show.legend = FALSE) +
  geom_line(aes(y = predicted, colour = .model), linewidth = 0.6) +
  geom_line(aes(y = observed, colour = "Observed"), linewidth = 0.55) +
  facet_wrap(~ .model, ncol = 1) +
  scale_colour_manual(name = NULL, values = model_colours) +
  scale_fill_manual(values = model_colours) +
  labs(x = NULL, y = "Weekly LRI admissions") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey96", colour = NA),
        panel.border = element_rect(colour = "grey55", fill = NA))
ggsave(file.path(val_output_dir, "Figure_2.png"), figure_2, width = 10, height = 8, dpi = 600)
ggsave(file.path(val_output_dir, "Figure_2.pdf"), figure_2, width = 10, height = 8)

# Compare these outputs with the manuscript before replacing reported values.
# Updated development choices can change validation scores and W1-W3 labels.
