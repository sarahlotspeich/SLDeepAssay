`SLDeepAssay`: A package for maximum likelihood estimation from serial
dilution assays with partial deep sequencing
================
Sarah C. Lotspeich and Brian D. Richardson

## Installation

Installation of the `SLDeepAssay` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done in the following way.

``` r
# Install the package
devtools::install_github(repo = "sarahlotspeich/SLDeepAssay", 
                         ref = "main")
```

``` r
# Then load it
library(SLDeepAssay)
```

The `SLDeepAssay` package contains parallel functions designed to
analyze (or generate) assay data from single or multiple dilution
levels. A common naming convention used throughout is to use the
suffixes `_sd` and `_md` to denote functions for single- and
multiple-dilution assay data, respectively.

## Built-In Functions to Simulate Assay Data

We begin by demonstrating how to simulate assay data from a single
dilution level using the `simulate_assay_sd()` function. There are 4
required arguments,

- `M`: the total number of wells (a scalar),
- `n`: the number of existing distinct viral lineages (DVL) (a
  scalar),  
- `lambda`: the DVL-specific rate parameters (a vector of length `n`
  with all elements $>0$)
- `q`: the proportion of p24-positive wells to be deep sequenced (a
  scalar between $0$ and $1$),

and 2 optional arguments,

- `remove_undetected`: a logical indicator for whether DVL that are not
  detected in any wells should be deleted (`TRUE`) or not (`FALSE`), and
- `seed = NULL`: the desired random seed (a scalar).

The function returns a list containing two versions of the simulated
data, corresponding to the quantitative viral outgrowth assay (QVOA)
(`$any_DVL`) and Ultra Deep Sequencing Assay of the Outgrowth Virus
(UDSA) (`$DVL_specific`) measures discussed in the manuscript.

For demonstration, we simulate a single-dilution assay of `M = 12` total
wells (each with `dilution = 1` million cells in it), with `n = 6`
underlying DVL, and assume that `q = 0.5` (i.e., 50%) of p24-positive
wells will undergo deep sequencing. To specify `lambda`, recall the
connection for a single dilution level between the DVL-specific rates
and the overall IUPM = `sum(lambda)/dilution`. Therefore, if we want to
simulate data where the overall IUPM = 0.5, we set
`lambda = rep(0.5 / 6, 6)` (or any other way to distribute a total mass
of 0.5 across the 6 DVL).

We now demonstrate a few ways to simulate, depending on our choices for
the optional arguments.

``` r
assay = simulate_assay_sd(M = 12, 
                          lambda = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = TRUE)
assay
```

    ## $any_DVL
    ##  [1] 1 1 1 0 0 0 0 0 0 1 1 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    0    0    1    0    0    0    0    0    0    NA    NA    NA
    ## [2,]    0    1    0    0    0    0    0    0    0    NA    NA    NA
    ## [3,]    1    0    0    0    0    0    0    0    0    NA    NA    NA

A few things to notice here…

- `assay$any_DVL` is a vector containing indicators for each well, where
  `1` indicates that the well tested positive for at least one DVL in
  the QVOA.
- `assay$DVL_specific` only has 2 rows when we set `n = 6`. Where did
  they go?? Since `remove_undetected = TRUE` here, the other 4 rows were
  deleted because they must not have been detected in either of the 2
  p24-positive and deep sequenced wells. (We will confirm this in a
  moment.)
- `assay$DVL_specific` has columns of `NA` values at the end. This isn’t
  an error! The wells represented by columns 10–12 are the $M_P - m$
  p24-positive but not deep sequenced wells. They are correctly captured
  as being infected with *something* in `assay$any_DVL`, though.

### The `remove_undetected` option

Now what if I leave the `seed` the same, but instead set
`remove_undetected = FALSE`?

``` r
assay = simulate_assay_sd(M = 12,
                          lambda = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = FALSE)
assay
```

    ## $any_DVL
    ##  [1] 1 1 0 0 0 0 0 0 0 0 1 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    0    1    0    0    0    0    0    0    0     0    NA    NA
    ## [2,]    0    0    0    0    0    0    0    0    0     0    NA    NA
    ## [3,]    1    0    0    0    0    0    0    0    0     0    NA    NA
    ## [4,]    0    0    0    0    0    0    0    0    0     0    NA    NA
    ## [5,]    0    0    0    0    0    0    0    0    0     0    NA    NA
    ## [6,]    0    1    0    0    0    0    0    0    0     0    NA    NA

Notice that `assay$any_DVL` is the same as before, which makes sense
since it’s the same random seed. However, `assay$DVL_specific` now has
all `n = 6` rows, even though rows 1–2 and 4–5 contain all 0s.

### The `seed` option

When should I bother to specify the `seed` myself? Well, for
reproducibility it’s a great idea! That way, every time you run your
code, you know that you’ll be getting the same simulated assay. Suppose
I run my code once, change my mind about something, and run it again. If
I don’t set the seed, I’ll be dealing with completely different
datasets!

``` r
assay = simulate_assay_sd(M = 12,
                          lambda = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = TRUE)
assay
```

    ## $any_DVL
    ##  [1] 1 1 0 0 0 0 0 0 0 0 0 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    1    0    0    0    0    0    0    0    0     0     0    NA
    ## [2,]    0    1    0    0    0    0    0    0    0     0     0    NA

``` r
assay = simulate_assay_sd(M = 12,
                          lambda = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = TRUE)
assay
```

    ## $any_DVL
    ##  [1] 1 1 0 0 0 0 0 0 0 0 0 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    1    1    0    0    0    0    0    0    0     0     0    NA

It tends to be safer to set the `seed`, but if you’re using this
function inside of a loop or something as part of a large-scale
simulation study it might make more sense to do so outside of the
function.

### Multiple dilution levels

Now if we instead want to simulate assay data from multiple dilution
levels, we would use the `simulate_assay_md()` function. There are still
4 required arguments, but now they are

- `M`: the total number of wells (a scalar),
- `tau`:
- `q`: the proportion of p24-positive wells to be deep sequenced (a
  scalar between $0$ and $1$), and
- `dilutions`:

The 2 optional arguments are the same as with the `simulate_assay_sd()`
function.

## Built-In Functions for Estimation

Now that we have our assay data, we turn to primary purpose of the
`SLDeepAssay` package: to efficiently estimate the rate of infection
from dilution assays with deep sequencing data using maximum likelihood
estimation (MLE). We can obtain point estimates, standard errors, and
confidence intervals all with a single call to the
`fit_SLDeepAssay_sd()` function. There are actually two ways to supply
the assay data for estimation.

### Supplying full assay data

The first argument `assay` expects a matrix like the one
`simulate_assay_sd()` returned in `assay$DVL_specific` above. If this is
supplied, then the only other *required* argument is `dilution`, which
is the number of cells per well (in millions). Then, you can skip
parameters `M`, `n`, `MP`, `m`, and `Y` (they will be discussed in the
next section), and all that is left are optional arguments:

- `corrected`: a logical indicator for whether the bias-corrected MLE
  should be returned (`TRUE`) or not (`FALSE`),
- `maxit`: the maximum number of iterations allowed to find the MLE,
- `lb`: the lower bound for the parameters, and
- `ub`: the upper bound for the parameters.

The last three optional parameters are all passed to the `optim()`
function, which is used to find the MLE inside of `simulate_assay_sd()`.
Recall that data were simulated for wells with 1 million cells per well,
so `dilution = 1` and the MLE can be interpreted as the infectious units
per million (IUPM).

We now estimate the IUPM based on the last `assay` generated above,
leaving all optional arguments at their defaults.

``` r
res = fit_SLDeepAssay_sd(assay = assay$DVL_specific, 
                         dilution = 1)
res
```

    ## $mle
    ## [1] 0.287682
    ## 
    ## $se
    ## [1] 0.1667235
    ## 
    ## $ci
    ## [1] 0.09238609 0.89581604
    ## 
    ## $mle.bc
    ## [1] 0.24143
    ## 
    ## $ci.bc
    ## [1] 0.0623706 0.9345501

Finally, we interpret the output as follows.

- `res$mle` is the uncorrected MLE.
- `res$se` is the estimated standard error (assumed to be the same for
  both the bias corrected and uncorrected MLE).
- `res$ci` is the 95% confidence interval (CI) for the uncorrected MLE.
- `res$mle.bc` is the bias-corrected MLE.
- `res$ci.bc` is the 95% CI for the bias-corrected MLE.

### Supplying summarized assay data
