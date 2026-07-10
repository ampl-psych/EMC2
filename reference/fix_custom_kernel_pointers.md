# Reset pointers of custom C++ trend kernels to an emc object

When an emc object is loaded from disk, or returned by forked processes,
the pointers to custom kernels need to be re-created. This is a
convenience function to do this.

## Usage

``` r
fix_custom_kernel_pointers(emc, pointer_source, model_number = 1)
```

## Arguments

- emc:

  A target emc object with missing pointers

- pointer_source:

  Either a trend list with correct pointers, or an emc object with
  correct pointers

- model_number:

  If `pointer_source` is a raw trend list and `emc` is a joint model,
  specifies which model the trend belongs to. Defaults to 1.

## Value

An emc object with the custom pointers re-instated.
