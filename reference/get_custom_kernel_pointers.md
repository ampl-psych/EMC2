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

A list of custom pointers. The list is of the same size as total number
of trends; trends without a custom kernel return NULL.
