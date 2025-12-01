# Split fMRI Timeseries Data by ROI Columns

This function splits a timeseries data frame containing multiple ROI
columns into a list of data frames, where each data frame contains the
common columns (subjects, run, time) and one ROI column.

## Usage

``` r
split_timeseries(timeseries, columns = NULL)
```

## Arguments

- timeseries:

  A data frame containing fMRI timeseries data with required columns
  'subjects', 'run', and 'time', plus one or more ROI columns.

- columns:

  A character vector specifying which columns to split by. If NULL
  (default), all columns except 'subjects', 'run', and 'time' will be
  used.

## Value

A named list of data frames, where each data frame contains the common
columns (subjects, run, time) and one ROI column. The names of the list
elements correspond to the ROI column names.

## Examples

``` r
# Create a simple example timeseries with multiple ROIs
set.seed(123)
n_frames <- 100

# Create a data frame with multiple ROIs
timeseries <- data.frame(
  subjects = rep(1, n_frames),
  run = rep(1, n_frames),
  time = seq(0, n_frames-1),
  ROI1 = rnorm(n_frames),
  ROI2 = rnorm(n_frames),
  ROI3 = rnorm(n_frames)
)

# Split the timeseries by all ROI columns
split_data <- split_timeseries(timeseries)
```
