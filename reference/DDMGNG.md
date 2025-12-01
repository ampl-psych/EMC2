# The GNG (go/nogo) Diffusion Decision Model

In the GNG paradigm one of the two possible choices results in a
response being withheld (a non-response), which is indicated in the data
by an NA for the rt, with the corresponding level of the R (response)
factor still being specified. For example, suppose the go response is
coded as "yes" and nogo is coded as "no", then for a non-response (R,rt)
= ("no",NA) and for a response e.g., (R,rt) = ("yes",1.36). The GNG
paradigm must also have a response

## Usage

``` r
DDMGNG()
```

## Value

A model list with all the necessary functions to sample

## Details

The model used is described in the following paper, with the addition of
modeling the TIMEOUT (which is considered but not used in this paper).

Gomez, P., Ratcliff, R., & Perea, M. (2007). A Model of the Go/No-Go
Task. Journal of Experimental Psychology: General, 136(3), 389–413.
https://doi.org/10.1037/0096-3445.136.3.389

The likelihood of non-responses requires and evaluation of the DDM cdf,
specifically 1 - p(hitting the yes boundary before TIMEOUT).

To use these models three functions must be supplied in the design's
function argument with the names TIMEOUT, Rnogo and Rgo. For example,
assuming a 2.5 second timeout, and R factor with levels c("no","yes")
and "no" mapping to a non-response.

TIMEOUT=function(d)rep(2.5,nrow(d))
Rnogo=function(d)factor(rep("no",nrow(d)),levels=c("no","yes"))
Rgo=function(d)factor(rep("yes",nrow(d)),levels=c("no","yes")))

See the help for DDM for further details. At present this model is not
fully implemented in C, so is a little slower to use than the DDM, but
not greatly.

## window (i.e., a length of time, TIMEOUT period, after which withholding is

assumed).

## Examples

``` r
dGNG <- design(Rlevels = c("left","right"),
               factors=list(subjects=1,S=c("left","right")),
               functions=list(
               TIMEOUT=function(d)rep(2.5,nrow(d)),
               # no go response level
               Rnogo=function(d)factor(rep("left",nrow(d)),levels=c("left","right")),
               # go response level
               Rgo=function(d)factor(rep("right",nrow(d)),levels=c("left","right"))),
               formula=list(v~S,a~1, Z~1, t0~1),
               model=DDMGNG)
#> Parameter(s) SZ, s, st0, sv not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v"        "v_Sright" "a"        "Z"        "t0"      
#> 
#>  Design Matrices: 
#> $v
#>      S v v_Sright
#>   left 1        0
#>  right 1        1
#> 
#> $a
#>  a
#>  1
#> 
#> $Z
#>  Z
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 
#> $SZ
#>  SZ
#>   1
#> 
#> $s
#>  s
#>  1
#> 
#> $st0
#>  st0
#>    1
#> 
#> $sv
#>  sv
#>   1
#> 

p_vector <- sampled_pars(dGNG)
```
