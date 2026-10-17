# Weekly pediatric LRI forecasting

Code for the main analysis and hospital-specific sensitivity analyses.

**Data:** The dataset is not publicly available and is shared only upon request.

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
