# Apply High-Pass Filtering to fMRI Data

This function applies high-pass filtering to fMRI data to remove
low-frequency noise and drift. It supports two filtering methods: cosine
basis functions and polynomial regressors.

## Usage

``` r
high_pass_filter(X, high_pass_model = "cosine", frame_times = NULL, ...)
```

## Arguments

- X:

  A data frame or matrix containing the data to be filtered. If it
  contains columns 'subjects' and 'run', the function will apply
  filtering separately for each subject-run combination.

- high_pass_model:

  A character string specifying the high-pass filtering method. Options
  are 'cosine' (default) or 'poly' for polynomial regressors.

- frame_times:

  A numeric vector of time points for each frame. If NULL, the function
  will attempt to extract this from a 'time' column in X.

- ...:

  Additional arguments passed to the function.

## Value

A data frame or matrix with the same structure as X, but with
high-frequency components removed from the data columns.

## Examples

``` r
# Create a simple example data frame with drift
set.seed(123)
n_frames <- 100
time <- seq(0, 99)

# Create a signal with low-frequency drift
drift <- 0.1 * time
signal <- sin(2 * pi * 0.1 * time) + drift
noise <- rnorm(n_frames, 0, 0.5)
data <- signal + noise

# Create a data frame
df <- data.frame(
  time = time,
  signal = data
)

# Apply high-pass filtering using cosine basis functions
filtered_df <- high_pass_filter(df, high_pass_model = "cosine")
```
