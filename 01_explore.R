# ============================================================
# 01_explore.R
# Exploratory analysis of weekly pediatric LRI admissions
# ============================================================
# Run from the project root after checking the paths below.
# This script describes the data; it does not select forecast models.
# Original observations, including zero-admission weeks, are retained.
# Figure 1 displays the previously reported, fixed regime boundaries.

library(fpp3)

# Recent feasts versions provide time-series plotting through ggtime.
# Older fpp3 installations already expose these plotting functions.
if (requireNamespace("ggtime", quietly = TRUE)) library(ggtime)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

# Point this to the intended version of the input data.
# For the future GitHub repository, use "data/weekly_master.csv".
data_path <- "data/weekly_master.csv"
output_dir <- file.path("results", "explore")
development_end <- as.Date("2022-12-26")

# Display settings only: these reported dates are not re-estimated here.
# Update from the frozen specification if breakpoint decisions change later.
figure_breakpoints <- as.Date(c("2020-02-17", "2022-04-11"))

series_labels <- c(
  MAIN = "Combined LRI admissions",
  ND2 = "Children's Hospital 2",
  BND = "Hospital for Tropical Diseases"
)

weather_labels <- c(
  temp_mean = "Temperature (degrees C)",
  rh_mean = "Relative humidity (%)",
  rain_sum = "Rainfall (mm/week)"
)

# ------------------------------------------------------------
# 2. Read, order and check the weekly data
# ------------------------------------------------------------

if (!file.exists(data_path)) stop("Input file not found: ", data_path)
data <- read.csv(data_path, stringsAsFactors = FALSE)

required_columns <- c(
  "week_start", "iso_year", "week_id", "phase",
  "y_lri", "y_hosp1", "y_hosp2", "temp_mean", "rh_mean", "rain_sum"
)
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

# Accept the original month/day/year dates or unambiguous ISO dates.
data <- data |>
  mutate(
    week_start = as.Date(lubridate::parse_date_time(
      week_start, orders = c("mdy", "ymd"), quiet = TRUE
    ))
  ) |>
  arrange(week_start)

numeric_columns <- c(
  "iso_year", "week_id", "y_lri", "y_hosp1", "y_hosp2",
  "temp_mean", "rh_mean", "rain_sum"
)

stopifnot(
  nrow(data) > 104,
  !anyNA(data[required_columns]),
  all(vapply(data[numeric_columns], is.numeric, logical(1))),
  all(vapply(data[numeric_columns], function(x) all(is.finite(x)), logical(1))),
  !anyDuplicated(data$week_start),
  all(diff(data$week_start) == 7),
  all(lubridate::wday(data$week_start, week_start = 1) == 1),
  all(data$iso_year == lubridate::isoyear(data$week_start)),
  all(data$week_id == seq_len(nrow(data))),
  all(data$y_lri == data$y_hosp1 + data$y_hosp2),
  all(as.matrix(data[c("y_lri", "y_hosp1", "y_hosp2")]) >= 0),
  all(data$phase %in% c("development", "validation", "test")),
  all((data$phase == "development") == (data$week_start <= development_end)),
  min(data$week_start) == as.Date("2015-01-05"),
  max(data$week_start) == as.Date("2024-12-23")
)

# Preserve source phase labels for compatibility with the existing scripts.
# For reporting, pool source validation/test into the 2023-2024 period.
data <- data |>
  mutate(
    phase = factor(phase, levels = c("development", "validation", "test")),
    analysis_period = if_else(
      phase == "development", "Development", "Validation"
    )
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Descriptive summaries for all three admission series
# ------------------------------------------------------------

# Long format avoids repeating the same summaries for each hospital.
outcome_data <- data |>
  select(week_start, iso_year, phase, analysis_period, y_lri, y_hosp1, y_hosp2) |>
  pivot_longer(c(y_lri, y_hosp1, y_hosp2), names_to = "outcome", values_to = "value") |>
  mutate(series = recode(outcome, y_lri = "MAIN", y_hosp1 = "ND2", y_hosp2 = "BND"))

# Input may already be grouped by series, year or analysis period.
summarise_values <- function(df) {
  df |>
    summarise(
      n_weeks = n(),
      total = sum(value),
      mean = mean(value),
      sd = sd(value),
      median = median(value),
      q1 = quantile(value, 0.25),
      q3 = quantile(value, 0.75),
      min = min(value),
      max = max(value),
      n_zero = sum(value == 0),
      .groups = "drop"
    )
}

outcome_by_year <- outcome_data |>
  group_by(series, iso_year) |>
  summarise_values()

outcome_by_period <- outcome_data |>
  group_by(series, analysis_period) |>
  summarise_values()

weather_data <- data |>
  select(week_start, analysis_period, all_of(names(weather_labels))) |>
  pivot_longer(all_of(names(weather_labels)), names_to = "variable", values_to = "value")

# Totals are omitted because they are not useful for temperature/humidity.
weather_overall <- weather_data |>
  group_by(variable) |>
  summarise_values() |>
  select(-total)

weather_by_period <- weather_data |>
  group_by(variable, analysis_period) |>
  summarise_values() |>
  select(-total)

# List consecutive zero runs without assigning a cause or recoding them as NA.
find_zero_runs <- function(df) {
  df |>
    arrange(week_start) |>
    mutate(
      is_zero = value == 0,
      run_id = cumsum(is_zero != lag(is_zero, default = FALSE))
    ) |>
    filter(is_zero) |>
    group_by(run_id) |>
    summarise(
      start_week = min(week_start),
      end_week = max(week_start),
      n_weeks = n(),
      .groups = "drop"
    ) |>
    select(-run_id)
}

zero_runs <- outcome_data |>
  group_by(series) |>
  group_modify(~ find_zero_runs(.x)) |>
  ungroup()

# ------------------------------------------------------------
# 4. Development-only Guerrero estimates
# ------------------------------------------------------------

# Use the full development history, as in the original exploration script.
# No validation observations contribute to these estimates.
lambda_summary <- outcome_data |>
  filter(phase == "development") |>
  group_by(series) |>
  summarise(
    n_development = n(),
    lambda_estimated = guerrero(value + 1, .period = 52),
    .groups = "drop"
  )

print(lambda_summary)

# Retain the full Guerrero estimates for the subsequent analysis steps.
# Do not substitute manually entered or rounded values in downstream scripts.
# Keep the legacy object names needed by the existing breakpoint script.
data_y1 <- data |> mutate(y_hosp = y_hosp1)
data_y2 <- data |> mutate(y_hosp = y_hosp2)
lambda_nd2 <- lambda_summary$lambda_estimated[lambda_summary$series == "ND2"]
lambda_bnd <- lambda_summary$lambda_estimated[lambda_summary$series == "BND"]
lambda_guerrero <- lambda_summary |>
  filter(series == "MAIN") |>
  transmute(lambda = lambda_estimated)
lambda_main <- lambda_guerrero$lambda

# ------------------------------------------------------------
# 5. Table 1 and explicit sample counts
# ------------------------------------------------------------

# Lagged exposures define the descriptive sample only.
# Table 1 still summarizes raw weekly admissions and weather measurements.
# Do not overwrite data: later stages need their own explicit sample rules.
table_data <- data |>
  mutate(
    T01 = (temp_mean + lag(temp_mean, 1)) / 2,
    H02 = (rh_mean + lag(rh_mean, 1) + lag(rh_mean, 2)) / 3,
    R02 = rain_sum + lag(rain_sum, 1) + lag(rain_sum, 2)
  ) |>
  filter(!is.na(T01), !is.na(H02), !is.na(R02))

format_iqr <- function(x, digits = 0) {
  sprintf(
    "%.*f [%.*f-%.*f]",
    digits, median(x), digits, quantile(x, 0.25), digits, quantile(x, 0.75)
  )
}

table1 <- bind_rows(
  table_data |> mutate(analysis_period = "Overall"),
  table_data
) |>
  group_by(analysis_period) |>
  summarise(
    Weeks = as.character(n()),
    `LRI admissions, median [IQR]` = format_iqr(y_lri),
    `Temperature (degrees C), mean +/- SD` = sprintf("%.1f +/- %.1f", mean(temp_mean), sd(temp_mean)),
    `Humidity (%), mean +/- SD` = sprintf("%.1f +/- %.1f", mean(rh_mean), sd(rh_mean)),
    `Rainfall (mm/week), median [IQR]` = format_iqr(rain_sum, digits = 1),
    .groups = "drop"
  ) |>
  pivot_longer(-analysis_period, names_to = "Characteristic", values_to = "Value") |>
  pivot_wider(names_from = analysis_period, values_from = Value) |>
  select(Characteristic, Overall, Development, Validation)

sample_counts <- bind_rows(
  data |> count(analysis_period, name = "n_weeks") |> mutate(sample = "Source series"),
  table_data |> count(analysis_period, name = "n_weeks") |> mutate(sample = "Table 1")
) |>
  select(sample, analysis_period, n_weeks)

print(table1)
print(sample_counts)

# ------------------------------------------------------------
# 6. Time-series plots and development-only diagnostics
# ------------------------------------------------------------

# Keep the original plotting scope: the combined descriptive plots show the
# full series; hospital seasonal/monthly plots show development data only.
# All ACF, PACF, CCF and correlations use development data only.
explore_series <- function(outcome_column, label, development_plots = FALSE) {
  ts_data <- data |>
    mutate(y = .data[[outcome_column]], week = yearweek(week_start)) |>
    as_tsibble(index = week)
  dev <- ts_data |> filter(phase == "development")
  plot_data <- if (development_plots) dev else ts_data
  scope <- if (development_plots) "development" else "2015-2024"

  plots <- list(
    time = autoplot(plot_data, y) +
      labs(title = paste(label, scope, sep = ": "), x = NULL, y = "Admissions"),
    seasonal = gg_season(plot_data, y, period = "1 year", labels = "right", labels_repel = TRUE) +
      labs(title = paste(label, "seasonal pattern"), x = NULL, y = "Admissions") +
      theme(legend.position = "none"),
    monthly = plot_data |>
      mutate(month = lubridate::month(week_start, label = TRUE, abbr = TRUE)) |>
      ggplot(aes(month, y)) + geom_boxplot() +
      labs(title = paste(label, "monthly distribution"), x = NULL, y = "Admissions"),
    acf = ACF(dev, y, lag_max = 104) |> autoplot() +
      labs(title = paste(label, "development ACF")),
    pacf = PACF(dev, y, lag_max = 104) |> autoplot() +
      labs(title = paste(label, "development PACF")),
    ccf_temperature = CCF(dev, y, temp_mean, lag_max = 12) |> autoplot() +
      labs(title = paste(label, "development CCF: temperature")),
    ccf_humidity = CCF(dev, y, rh_mean, lag_max = 12) |> autoplot() +
      labs(title = paste(label, "development CCF: humidity")),
    ccf_rainfall = CCF(dev, y, rain_sum, lag_max = 12) |> autoplot() +
      labs(title = paste(label, "development CCF: rainfall"))
  )

  correlation <- dev |>
    as_tibble() |>
    select(y, temp_mean, rh_mean, rain_sum) |>
    cor()

  list(plots = plots, correlation = correlation)
}

exploration <- list(
  MAIN = explore_series("y_lri", series_labels[["MAIN"]]),
  ND2 = explore_series("y_hosp1", series_labels[["ND2"]], development_plots = TRUE),
  BND = explore_series("y_hosp2", series_labels[["BND"]], development_plots = TRUE)
)

# Full-period hospital time plots allow direct inspection of service disruptions.
plot_hospitals <- ggplot(
  filter(outcome_data, series != "MAIN"), aes(week_start, value)
) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Weekly admissions") +
  theme_bw()

# Figure 1 follows the original manuscript layout: four panels, three regime
# backgrounds, fixed breakpoint markers and the validation boundary.
figure_labels <- c(
  y_lri = "LRI admissions",
  temp_mean = "Temperature (\u00B0C)",
  rh_mean = "Relative humidity (%)",
  rain_sum = "Rainfall (mm/week)"
)
figure_colours <- c(
  y_lri = "#A65A55", temp_mean = "#C28A5B",
  rh_mean = "#6F918A", rain_sum = "#6C879C"
)

overview_data <- table_data |>
  select(week_start, y_lri, all_of(names(weather_labels))) |>
  pivot_longer(-week_start, names_to = "variable", values_to = "value") |>
  mutate(variable = factor(variable, levels = names(figure_labels)))

validation_start <- min(table_data$week_start[table_data$analysis_period == "Validation"])
bp1 <- figure_breakpoints[1]
bp2 <- figure_breakpoints[2]
stopifnot(bp1 < bp2, bp2 < validation_start)

# Show the year at ISO week 1 and W27 at the mid-year tick, as in the old code.
axis_data <- table_data |>
  mutate(
    iso_week = lubridate::isoweek(week_start),
    year_label = as.character(lubridate::isoyear(week_start)),
    label = if_else(iso_week == 1, year_label, "W27")
  ) |>
  filter(iso_week %in% c(1, 27))

boundary_lines <- tibble(
  date = c(bp1, bp2, validation_start),
  type = factor(
    c("Breakpoint 1", "Breakpoint 2", "Validation boundary"),
    levels = c("Breakpoint 1", "Breakpoint 2", "Validation boundary")
  )
)

plot_overview <- ggplot(overview_data, aes(week_start, value)) +
  annotate(
    "rect", xmin = min(table_data$week_start), xmax = bp1,
    ymin = -Inf, ymax = Inf, fill = "#F5ECEA", alpha = 0.30
  ) +
  annotate(
    "rect", xmin = bp1, xmax = bp2,
    ymin = -Inf, ymax = Inf, fill = "#EEF2E8", alpha = 0.30
  ) +
  annotate(
    "rect", xmin = bp2, xmax = max(table_data$week_start),
    ymin = -Inf, ymax = Inf, fill = "#EAF0F5", alpha = 0.30
  ) +
  geom_line(aes(colour = variable), linewidth = 0.5, show.legend = FALSE) +
  geom_vline(
    data = boundary_lines,
    aes(xintercept = date, linetype = type),
    colour = "grey20", linewidth = 0.55, show.legend = TRUE
  ) +
  facet_wrap(
    ~ variable, ncol = 1, scales = "free_y",
    labeller = labeller(variable = figure_labels)
  ) +
  scale_colour_manual(values = figure_colours, guide = "none") +
  scale_linetype_manual(values = c(
    "Breakpoint 1" = "dashed",
    "Breakpoint 2" = "dotdash",
    "Validation boundary" = "solid"
  )) +
  scale_x_date(breaks = axis_data$week_start, labels = axis_data$label) +
  labs(
    x = NULL, y = NULL, linetype = NULL,
    caption = "Figure 1. Time series of weekly LRI admissions and meteorological variables"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 8),
    plot.caption = element_text(face = "bold", hjust = 0.5, margin = margin(t = 10))
  )

plot_weather <- weather_data |>
  mutate(variable = factor(variable, levels = names(weather_labels), labels = unname(weather_labels))) |>
  ggplot(aes(week_start, value)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_bw()

# ------------------------------------------------------------
# 7. Export tables and plots
# ------------------------------------------------------------

tables <- list(
  outcome_by_year = outcome_by_year,
  outcome_by_period = outcome_by_period,
  weather_overall = weather_overall,
  weather_by_period = weather_by_period,
  zero_admission_runs = zero_runs,
  lambda_summary = lambda_summary,
  sample_counts = sample_counts,
  Table1 = table1
)

for (name in names(tables)) {
  write.csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
}
write.table(table1, file.path(output_dir, "Table1.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

for (series in names(exploration)) {
  result <- exploration[[series]]
  write.csv(result$correlation, file.path(output_dir, paste0(series, "_correlation.csv")))

  for (plot_name in names(result$plots)) {
    ggsave(
      file.path(output_dir, paste0(series, "_", plot_name, ".pdf")),
      plot = result$plots[[plot_name]], width = 9, height = 5
    )
  }
}

ggsave(file.path(output_dir, "hospital_time_series.pdf"), plot_hospitals, width = 10, height = 6)
ggsave(file.path(output_dir, "weather_time_series.pdf"), plot_weather, width = 10, height = 7)
ggsave(file.path(output_dir, "Figure_1.pdf"), plot_overview, width = 10, height = 8)
ggsave(file.path(output_dir, "Figure_1.png"), plot_overview, width = 10, height = 8, dpi = 300)

print(plot_overview)
message("Exploration outputs: ", normalizePath(output_dir, winslash = "/"))

# Outputs remain available as data frames/lists for interactive inspection.
# Correlations and CCF plots are descriptive, not causal estimates or validation
# results. Do not use validation-period patterns to reselect model specifications.
