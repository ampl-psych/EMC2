# Reset pointers of custom C++ trend kernels to an emc object

When an emc object is loaded from disk, or returned by forked processes,
the pointers to custom kernels need to be re-created. This is a
convenience function to do this.

## Usage

``` r
fix_custom_kernel_pointers(emc, pointer_source)
```

## Arguments

- emc:

  A target emc object with missing pointers

- pointer_source:

  Either a trend object with correct pointers or another emc object with
  correct pointers

## Value

An emc object with the custom pointers re-instated.
