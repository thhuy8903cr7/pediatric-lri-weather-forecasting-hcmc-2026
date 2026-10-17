# ============================================================
# 02_breakpoint.R
# Development-only structural breakpoint analysis for pediatric LRI
# ============================================================
# Run 01_explore.R first, in the same R session and from the project root.
# Required inputs: data and lambda_summary from that script.
# No data are reloaded and no lambda values or breakpoint dates are entered here.
# Zero-admission weeks are retained as recorded.

library(fpp3)
library(breakfast)
if (requireNamespace("ggtime", quietly = TRUE)) library(ggtime)

# ------------------------------------------------------------
# 1. Settings and inputs
# ------------------------------------------------------------

breakpoint_output_dir <- file.path("results", "breakpoints")
breakpoint_development_start <- as.Date("2015-01-05")
breakpoint_development_end <- as.Date("2022-12-26")
fourier_candidates <- 1:6
primary_seed <- 111L
stability_seeds <- 1:50

# Keep the original breakfast defaults for these two parameters.
# Other sol.wcm/model.gsa arguments also retain their package defaults.
wcm_intervals <- 100L
gsa_max_ar <- 10L

if (!exists("data", inherits = FALSE) || !exists("lambda_summary", inherits = FALSE)) {
  stop("Run 01_explore.R first. Objects 'data' and 'lambda_summary' are required.")
}
if (!is.data.frame(data) || !is.data.frame(lambda_summary)) {
  stop("The inputs from 01_explore.R must be data frames.")
}

required_columns <- c("week_start", "phase", "y_lri", "y_hosp1", "y_hosp2")
stopifnot(
  all(required_columns %in% names(data)),
  all(c("series", "lambda_estimated", "n_development") %in% names(lambda_summary))
)
stopifnot(
  inherits(data$week_start, "Date"),
  !anyNA(data[required_columns]),
  !anyDuplicated(data$week_start),
  nrow(lambda_summary) == 3,
  setequal(lambda_summary$series, c("MAIN", "ND2", "BND")),
  is.numeric(lambda_summary$lambda_estimated),
  all(is.finite(lambda_summary$lambda_estimated))
)

# Sort once before constructing any time series or mapping indices to dates.
breakpoint_data <- data |>
  as_tibble() |>
  arrange(week_start)

stopifnot(
  all(diff(breakpoint_data$week_start) == 7),
  all(lubridate::wday(breakpoint_data$week_start, week_start = 1) == 1),
  all((breakpoint_data$phase == "development") ==
        (breakpoint_data$week_start <= breakpoint_development_end))
)

development_data <- breakpoint_data |>
  filter(phase == "development")

stopifnot(
  min(development_data$week_start) == breakpoint_development_start,
  max(development_data$week_start) == breakpoint_development_end,
  all(lambda_summary$n_development == nrow(development_data)),
  all(vapply(development_data[c("y_lri", "y_hosp1", "y_hosp2")], is.numeric, logical(1)))
)

series_settings <- tibble(
  series = c("MAIN", "ND2", "BND"),
  outcome_column = c("y_lri", "y_hosp1", "y_hosp2")
) |>
  left_join(lambda_summary, by = "series")

dir.create(breakpoint_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Small helpers shared by all three series
# ------------------------------------------------------------

# breakfast returns the last observation of the previous segment.
# Report the next observation as the first week of the new regime.
make_breakpoint_table <- function(indices, weeks) {
  indices <- sort(as.integer(indices))
  stopifnot(!anyNA(indices), !anyDuplicated(indices), all(indices >= 1), all(indices < length(weeks)))

  tibble(
    cpt_index = indices,
    last_week_old_regime = weeks[indices],
    first_week_new_regime = weeks[indices + 1L]
  )
}

detect_breakpoints <- function(x, seed) {
  set.seed(seed)
  path <- sol.wcm(x, M = wcm_intervals)
  model.gsa(path, p.max = gsa_max_ar)
}

# Record an entire set of breakpoint dates, including runs with no breakpoint.
format_breakpoint_set <- function(indices, weeks) {
  dates <- weeks[sort(as.integer(indices)) + 1L]
  if (length(dates) == 0) return("None")
  paste(format(dates, "%Y-%m-%d"), collapse = "; ")
}

# ------------------------------------------------------------
# 3. Select nuisance seasonality and detect breakpoints
# ------------------------------------------------------------

analyse_breakpoints <- function(weekly_data, outcome_column, lambda) {
  dev <- weekly_data |>
    mutate(
      y = .data[[outcome_column]],
      y_bc = box_cox(y + 1, lambda),
      week = yearweek(week_start)
    ) |>
    as_tsibble(index = week)

  stopifnot(all(is.finite(dev$y)), all(dev$y >= 0), all(is.finite(dev$y_bc)))

  # Compare K = 1,...,6 on exactly the same development observations.
  # This nuisance Fourier order is separate from the forecasting model's K.
  fourier_fits <- lapply(fourier_candidates, function(k) {
    dev |>
      model(adjustment = TSLM(y_bc ~ trend() + fourier(period = "year", K = k)))
  })

  fourier_rank <- purrr::map2_dfr(
    fourier_candidates, fourier_fits,
    function(k, fit) glance(fit) |> mutate(K = k)
  ) |>
    select(K, AIC, AICc, BIC) |>
    arrange(AICc, K) |>
    mutate(rank = row_number(), delta_AICc = AICc - min(AICc)) |>
    select(rank, K, AIC, AICc, delta_AICc, BIC)

  stopifnot(all(is.finite(fourier_rank$AICc)))
  selected_K <- fourier_rank$K[1]
  selected_fit <- fourier_fits[[match(selected_K, fourier_candidates)]]

  # Join residuals by week so their positions remain aligned with source dates.
  adjustment_residuals <- residuals(selected_fit) |>
    as_tibble() |>
    select(week, .resid)
  stopifnot(!anyDuplicated(adjustment_residuals$week))

  diagnostic_data <- dev |>
    as_tibble() |>
    left_join(adjustment_residuals, by = "week") |>
    as_tsibble(index = week)
  stopifnot(nrow(diagnostic_data) == nrow(dev), all(is.finite(diagnostic_data$.resid)))

  x <- diagnostic_data$.resid
  primary_fit <- detect_breakpoints(x, primary_seed)
  primary_indices <- sort(as.integer(primary_fit$cpts))
  breakpoints <- make_breakpoint_table(primary_indices, dev$week_start)
  stopifnot(nrow(breakpoints) == primary_fit$no.of.cpt)

  # Use seed 111 for the primary analysis, as in the original code.
  # The 50 additional runs assess algorithmic stability; they do not select
  # whichever set of dates best matches the manuscript or validation results.
  seed_indices <- lapply(stability_seeds, function(seed) {
    fit <- detect_breakpoints(x, seed)
    indices <- sort(as.integer(fit$cpts))
    make_breakpoint_table(indices, dev$week_start)
    stopifnot(length(indices) == fit$no.of.cpt)
    indices
  })

  stability_runs <- tibble(
    seed = stability_seeds,
    n_cpt = lengths(seed_indices),
    breakpoint_set = vapply(
      seed_indices, format_breakpoint_set, character(1), weeks = dev$week_start
    ),
    matches_primary = vapply(
      seed_indices, function(indices) identical(indices, primary_indices), logical(1)
    )
  )

  # Keep every detected location, even when the number of breaks differs by seed.
  # Runs with zero breakpoints remain represented in stability_runs.
  stability_dates <- purrr::map2_dfr(stability_seeds, seed_indices, function(seed, indices) {
    make_breakpoint_table(indices, dev$week_start) |>
      mutate(seed = seed, .before = 1)
  })

  stability_counts <- stability_runs |>
    count(n_cpt, name = "n_seeds") |>
    mutate(proportion = n_seeds / length(stability_seeds))

  stability_sets <- stability_runs |>
    count(n_cpt, breakpoint_set, name = "n_seeds", sort = TRUE) |>
    mutate(proportion = n_seeds / length(stability_seeds))

  stability_locations <- stability_dates |>
    count(cpt_index, first_week_new_regime, name = "n_seeds", sort = TRUE) |>
    mutate(proportion = n_seeds / length(stability_seeds))

  # These are diagnostics, not a requirement that the residuals be white noise.
  # WCM.gSa allows serial dependence in the noise component.
  stopifnot(length(primary_fit$est) == length(x), all(is.finite(primary_fit$est)))
  diagnostic_data <- diagnostic_data |>
    mutate(
      estimated_mean_shift = primary_fit$est,
      residual_after_break = .resid - estimated_mean_shift
    )

  plots <- list(
    acf_before = ACF(diagnostic_data, .resid, lag_max = 104) |>
      autoplot() + labs(title = "Before removing detected mean shifts"),
    acf_after = ACF(diagnostic_data, residual_after_break, lag_max = 104) |>
      autoplot() + labs(title = "After removing detected mean shifts"),
    mean_shifts = ggplot(diagnostic_data, aes(week_start, .resid)) +
      geom_line(linewidth = 0.4, colour = "grey50") +
      geom_step(aes(y = estimated_mean_shift), linewidth = 0.6, colour = "#A65A55") +
      labs(x = NULL, y = "Adjusted Box-Cox residual", title = "Detected mean shifts") +
      theme_bw()
  )

  # Preserve exact estimates and dates for later development/validation scripts.
  # No fixed number of regimes is imposed; zero detected breaks means one regime.
  specification <- list(
    outcome_column = outcome_column,
    lambda = lambda,
    outcome_offset = 1,
    development_weeks = dev$week_start,
    nuisance_fourier_K = selected_K,
    primary_seed = primary_seed,
    stability_seeds = stability_seeds,
    wcm_M = wcm_intervals,
    gsa_p_max = gsa_max_ar,
    breakfast_version = as.character(packageVersion("breakfast")),
    breakpoint_indices = primary_indices,
    regime_start_dates = breakpoints$first_week_new_regime,
    breakpoint_convention = "First week of the new regime",
    n_regimes = length(primary_indices) + 1L
  )

  list(
    specification = specification,
    fourier_rank = fourier_rank,
    breakpoints = breakpoints,
    stability_runs = stability_runs,
    stability_dates = stability_dates,
    stability_counts = stability_counts,
    stability_sets = stability_sets,
    stability_locations = stability_locations,
    diagnostic_data = as_tibble(diagnostic_data),
    plots = plots
  )
}

# ------------------------------------------------------------
# 4. Apply the same procedure to MAIN, ND2 and BND
# ------------------------------------------------------------

breakpoint_results <- list()
for (i in seq_len(nrow(series_settings))) {
  series_id <- series_settings$series[i]
  message("Breakpoint analysis: ", series_id)
  breakpoint_results[[series_id]] <- analyse_breakpoints(
    weekly_data = development_data,
    outcome_column = series_settings$outcome_column[i],
    lambda = series_settings$lambda_estimated[i]
  )
}

breakpoint_spec <- lapply(breakpoint_results, function(result) result$specification)
breakpoint_summary <- purrr::imap_dfr(breakpoint_results, function(result, series_id) {
  spec <- result$specification
  tibble(
    series = series_id,
    lambda = spec$lambda,
    nuisance_fourier_K = spec$nuisance_fourier_K,
    n_development = length(spec$development_weeks),
    n_breakpoints = length(spec$breakpoint_indices),
    n_regimes = spec$n_regimes,
    primary_breakpoint_set = format_breakpoint_set(spec$breakpoint_indices, spec$development_weeks),
    n_stability_runs = nrow(result$stability_runs),
    n_matching_primary = sum(result$stability_runs$matches_primary)
  )
})

# ------------------------------------------------------------
# 5. Export results for review and subsequent model development
# ------------------------------------------------------------

table_names <- c(
  "fourier_rank", "breakpoints", "stability_runs", "stability_dates",
  "stability_counts", "stability_sets", "stability_locations", "diagnostic_data"
)
for (table_name in table_names) {
  combined_table <- purrr::imap_dfr(breakpoint_results, function(result, series_id) {
    result[[table_name]] |> mutate(series = series_id, .before = 1)
  })
  write.csv(
    combined_table,
    file.path(breakpoint_output_dir, paste0(table_name, ".csv")),
    row.names = FALSE
  )
}

write.csv(breakpoint_summary, file.path(breakpoint_output_dir, "breakpoint_summary.csv"), row.names = FALSE)
saveRDS(breakpoint_spec, file.path(breakpoint_output_dir, "breakpoint_spec.rds"))

for (series_id in names(breakpoint_results)) {
  for (plot_name in names(breakpoint_results[[series_id]]$plots)) {
    ggsave(
      file.path(breakpoint_output_dir, paste0(series_id, "_", plot_name, ".pdf")),
      plot = breakpoint_results[[series_id]]$plots[[plot_name]], width = 9, height = 5
    )
  }
}

print(breakpoint_summary)
message("Breakpoint outputs: ", normalizePath(breakpoint_output_dir, winslash = "/"))

# Next steps must read the generated specification instead of typing old dates.
# For a date b in regime_start_dates, the new regime starts at week_start >= b.
# Update Figure 1 and the manuscript only after checking the actual outputs.
# Stability across seeds describes algorithmic variability, not statistical
# confidence intervals for the breakpoint locations.
