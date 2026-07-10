# Extract pointers of custom C++ trend kernels from trend list or emc object

Extract the pointers so they can be re-added to an emc object after
loading it from disk.

## Usage

``` r
get_custom_kernel_pointers(input_data)
```

## Arguments

- input_data:

  Either an emc object or a trend list.

## Value

A list-of-lists of custom pointers, one sub-list per model. Trends
without a custom kernel return NULL. Returns NULL if no custom pointers
are found.
