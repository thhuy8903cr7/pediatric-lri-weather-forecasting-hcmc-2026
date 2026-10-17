# Weekly pediatric LRI forecasting

Code for the main analysis and hospital-specific sensitivity analyses.

**Data:** The dataset is not publicly available and is shared only upon request.

## Scripts

- **01_explore.R:** Prepares weekly data, summarizes admissions and weather, estimates Box-Cox parameters, and produces Table 1 and Figure 1.
- **02_breakpoint.R:** Detects development-period structural breakpoints for each outcome, evaluates stability across 50 seeds, and saves regime definitions.
- **03_dev.R:** Selects benchmark and dynamic regression models, screens weather predictors, ranks 26 weather combinations, checks residuals, and saves frozen specifications for MAIN, ND2, and BND.
- **04_val.R:** Performs weekly expanding-window validation in 2023-2024, calculates forecast accuracy, adds interval scores and bootstrap confidence intervals for MAIN, and exports validation tables and Figure 2.

## Run the code

1. Install the required R packages: `fpp3`, `purrr`, `ggtime`, `breakfast`, `prophet`, and `distillery`.
2. Set the working directory to the repository root and update `data_path` in `01_explore.R` to the location of the supplied dataset.
3. Run the four scripts in order, in the **same R session**:

```r
source("Main_analysis/01_explore.R")
source("Main_analysis/02_breakpoint.R")
source("Main_analysis/03_dev.R")
source("Main_analysis/04_val.R")
```

Tables, figures, and model specifications are saved in `results/`.
