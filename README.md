`SLDeepAssay`: A package for maximum likelihood estimation from serial
dilution assays with partial deep sequencing
================
Sarah C. Lotspeich and Brian D. Richardson

## Installation

Installation of `SLDeepAssay` from GitHub requires the
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
analyze (or generate) assay data from the quantitative viral outgrowth assay
(QVOA) and the ultra deep sequencing assay of the outgrowth virus
(UDSA). The methods implemented here are described in the described in
the paper, “Quantifying the HIV reservoir with dilution assays and deep
viral sequencing,” published in Biometrics (2024) at <https://academic.oup.com/biometrics/article/80/1/ujad018/7609164?login=false>.

A common naming convention used throughout is to use the suffixes `_sd`
and `_md` to denote functions for single- and multiple-dilution assay
data, respectively.

## Data Simulation

We begin by demonstrating how to simulate assay data from a single
dilution level using the `simulate_assay_sd()` function.

There are 3 required arguments,

- `M`: the total number of wells (a scalar),
- `tau`: the mean counts of infected cells per well with each DVL (a
  vector of with all elements $>0$)
- `q`: the proportion of p24-positive wells to be deep sequenced (a
  scalar between $0$ and $1$),

and 8 optional arguments,

- `u`: the total number of wells (a scalar),
- `sens_QVOA`, `spec_QVOA`: the sensitivity and specificity of the
  QVOA,respectively (scalars between 0 and 1, inclusive). Defaults are
  1, which correspond to perfect assays.
- `sens_UDSA`, `spec_UDSA`: the sensitivity and specificity of the
  UDSA,respectively (scalars between 0 and 1, inclusive). Defaults are
  1, which correspond to perfect assays.
- `sequence_all`: Logical, if `sequence_all = FALSE` (the default), then
  only p24-positive wells are considered for UDSA. If
  `sequence_all = TRUE` and `q = 1`, then all wells will undergo UDSA.
- `k`: overdispersion parameter (a positive scalar). Default is
  `k = Inf`, which corresponds to no overdispersion.
- `remove_undetected`: a logical indicator for whether DVL that are not
  detected in any wells should be deleted (`TRUE`) or not (`FALSE`).

The function returns a list containing two versions of the simulated
data, corresponding to the quantitative viral outgrowth assay (QVOA)
(`$any_DVL`) and Ultra Deep Sequencing Assay of the Outgrowth Virus
(UDSA) (`$DVL_specific`) measures.

For demonstration, we simulate a single-dilution assay of `M = 12` total
wells (each with `u = 1` million cells in it), with `n = 6` underlying
DVLs, and assume that `q = 0.5` (i.e., 50%) of p24-positive wells will
undergo deep sequencing. We also assume that both assays have perfect
sensitivity and specificity and that cell counts are not overdispersed.
Lastly, if we want to simulate data where the overall IUPM = 0.5, we set
`tau = rep(0.5 / 6, 6)` (or any other way to distribute a total mass of
0.5 across the 6 DVL).

``` r
set.seed(1)
assay_sd = simulate_assay_sd(M = 12, 
                          tau = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = TRUE)
assay_sd
```

    ## $any_DVL
    ##  [1] 1 1 0 0 0 0 0 0 0 0 0 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    0    1    0    0    0    0    0    0    0     0     0    NA
    ## [2,]    1    0    0    0    0    0    0    0    0     0     0    NA

A few things to notice here…

- `assay$any_DVL` is a vector containing indicators for each well, where
  `1` indicates that the well tested positive for at least one DVL in
  the QVOA.
- `assay$DVL_specific` only has 2 rows when we set `n = 6`. Where did
  they go? Since `remove_undetected = TRUE` here, the other 4 rows were
  deleted because they must not have been detected in either of the 2
  p24-positive and deep sequenced wells. (We will confirm this in a
  moment.)
- `assay$DVL_specific` has columns of `NA` values at the end. This isn’t
  an error! The wells represented by columns 10–12 are the $M_P - m$
  p24-positive but not deep sequenced wells. They are correctly captured
  as being infected with some DVL in `assay$any_DVL`, though.

### The `remove_undetected` option

Now what if I leave the `seed` the same, but instead set
`remove_undetected = FALSE`?

``` r
set.seed(1)
assay_sd2 = simulate_assay_sd(M = 12,
                          tau = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = FALSE)
assay_sd2
```

    ## $any_DVL
    ##  [1] 1 1 0 0 0 0 0 0 0 0 0 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    0    1    0    0    0    0    0    0    0     0     0    NA
    ## [2,]    1    0    0    0    0    0    0    0    0     0     0    NA
    ## [3,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [4,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [5,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [6,]    0    0    0    0    0    0    0    0    0     0     0    NA

Notice that `assay$any_DVL` is the same as before, which makes sense
since we used the same random number seed. However, `assay$DVL_specific`
now has all `n = 6` rows, even though rows 1–2 and 4–5 contain all 0s.

### Multiple dilution levels

Now if we instead want to simulate assay data from multiple dilution
levels, we can use the `simulate_assay_md()` function. There are still 4
required arguments, but now they are

- `M`: the total number of wells at each dilution level (a vector of
  length D)
- `tau`: DVL-specific IUPMs (a vector of length ). (Note: All elements
  in must be \> 0.)
- `q`: proportions of p24-positive wells that underwent UDSA at each
  dilution level (a vector of length D), and
- `u`: a vector of dilution levels (in millions of cells per well).

The 2 optional arguments `k` and `remove_undetected` are the same as
with the `simulate_assay_sd()` function. The function
`simulate_assay_md_imperfect()` function has additional arguments
`sens_QVOA`, `spec_QVOA`, `sens_UDSA`, `spec_UDSA` to generate data at
multiple dilution levels from imperfect assays.

For demonstration, we simulate a multiple-dilution assay with 3 dilution
levels, `M = c(8, 10, 12)` total wells in the 3 dilution levels
(`u = c(0.25, 0.5, 1)` million cells), with `n = 6` underlying DVLs, and
assume that `q = c(0, 0.5, 1)` (i.e., 0% of positive wells in the first
dilution level undergo deep sequencing, 50% of positive wellls at the
second, etc.). We again assume that both assays have perfect sensitivity
and specificity, that cell counts are not overdispersed, and that the
true IUPM is 0.5, split evenly among the 6 DVLs.

``` r
set.seed(1)
assay_md = simulate_assay_md(M = c(8, 10, 12), 
                          tau = rep(0.5 / 6, 6), 
                          q = c(0, 0.5, 1), 
                          u = c(0.25, 0.5, 1),
                          remove_undetected = TRUE)
assay_md
```

    ##      u  M n MN MP m   q Y1 Y2 Y3 Y4
    ## 1 0.25  8 4  7  1 0 0.0  0  0  0  0
    ## 2 0.50 10 4  8  2 1 0.5  0  0  0  1
    ## 3 1.00 12 4  9  3 3 1.0  2  1  2  0

Note that, in the multiple dilution setting, a data frame with summary
data of the simulated results is returned. This summary contains one row
per dilution level and the following columns: `M` (total number of
wells), `n` (number of DVLs, `MN` (number of p24-negative wells), `m`
(number of deep sequenced wells), `Y1`,…, `Yn` (counts of wells positive
for DVL i, (i = 1,…,n), and `u` (dilution levels, in millions of cells
per well).

## Estimation

Now that we have our assay data, we turn to the primary purpose of the
`SLDeepAssay` package: to efficiently estimate the number of infectious
cells per million using maximum likelihood estimation (MLE). We can
obtain point estimates, standard errors, and confidence intervals all
with a single call to the `fit_SLDeepAssay_sd()` function. There are two
ways to supply the assay data for estimation.

### Supplying full assay data

The first argument `assay` expects a matrix like the one
`simulate_assay_sd()` returned in `assay$DVL_specific` above. If this is
supplied, then the only other *required* argument is `u`, which is the
dilution level in millions of cells per well. Then, you can skip
parameters `M`, `MP`, `m`, and `Y` (they will be discussed in the next
section), and all that is left are optional arguments:

- `corrected`: a logical indicator for whether the bias-corrected MLE
  should be returned (`TRUE`) or not (`FALSE`),
- `maxit`: the maximum number of iterations allowed to find the MLE,

The `maxit` parameter is passed to the `optim()` function, which is used
to find the MLE inside of `simulate_assay_sd()`. Recall that data were
simulated for wells with 1 million cells per well, so `u = 1` and the
MLE can be interpreted as the infectious units per million (IUPM).

We now estimate the IUPM based on the single dilution data `assay_sd`
generated above, leaving all optional arguments at their defaults.

``` r
res_sd = fit_SLDeepAssay_sd(assay= assay_sd$DVL_specific, 
                            u = 1)
res_sd
```

    ## $mle
    ## [1] 0.273577
    ## 
    ## $se
    ## [1] 0.1581683
    ## 
    ## $ci
    ## [1] 0.08809639 0.84957352
    ## 
    ## $mle_bc
    ## [1] 0.2410677
    ## 
    ## $ci_bc
    ## [1] 0.06662735 0.87221920

Finally, we interpret the output as follows.

- `res$mle` is the uncorrected MLE.
- `res$se` is the estimated standard error (assumed to be the same for
  both the bias corrected and uncorrected MLE).
- `res$ci` is the 95% confidence interval (CI) for the uncorrected MLE.
- `res$mle_bc` is the bias-corrected MLE.
- `res$ci_bc` is the 95% CI for the bias-corrected MLE.

### Supplying summarized assay data

Instead of supplying the raw assay data in the argument `assay`, we can
supply summary data using the following parameters: - `M`: the total
number of wells (a scalar), - `MP`: the total number of p24-positive
wells (a scalar), - `m`: the total number of wells deep sequenced (a
scalar between 0 and `MP`), - `Y`: he numbers of wells (without missing
data) that were infected with each DVL (a vector of length `n`)

As expected, the two methods to supply the same information to
`fit_SLDeepAssay_sd()` give the same results.

``` r
M <- ncol(assay_sd$DVL_specific)
MP <- sum(colSums(assay_sd$DVL_specific) > 0, na.rm = TRUE)
m <- MP - sum(is.na(colSums(assay_sd$DVL_specific)))
Y <- rowSums(assay_sd$DVL_specific, na.rm = TRUE)

res_sd = fit_SLDeepAssay_sd(assay = NULL, 
                            M = M,
                            MP = MP,
                            m = m,
                            Y = Y,
                            u = 1)
res_sd
```

    ## $mle
    ## [1] 0.273577
    ## 
    ## $se
    ## [1] 0.1585768
    ## 
    ## $ci
    ## [1] 0.08783896 0.85206342
    ## 
    ## $mle_bc
    ## [1] 0.2278709
    ## 
    ## $ci_bc
    ## [1] 0.05825487 0.89134397

### Multiple Dilutions

We can similarly estimate the IUPM using the simulated multiple dilution
level data set `assay_md` using the function `fit_SLDeepAssay_md()`.

As in the single dilution case, there are two options to pass assay data
to `fit_SLDeepAssay_md()`. The first is to use summary data of the form
output by `simulate_assay_md()`, making this option easy to use with
simulated data.

The second is to supply the following two parameters:

- `assay`: a list of data frames with assay data from each dilution
  level.
- `u`: vector of dilution levels, in millions of cells per well.

This second option may be easier to use in practice, when raw data are
collected from real assays at multiple dilution levels. Below is an
example of using the summary data method. Note that the output
`simulate_assay_md()` is of the same format as that from
`simulate_assay_sd()`.

``` r
res_md <- fit_SLDeepAssay_md(assay_summary = assay_md)
res_md
```

    ## $mle
    ## [1] 0.4495328
    ## 
    ## $se
    ## [1] 0.1594212
    ## 
    ## $ci
    ## [1] 0.2243332 0.9008019
    ## 
    ## $mle_bc
    ## [1] 0.4352598
    ## 
    ## $ci_bc
    ## [1] 0.2123156 0.8923091

## Testing for Overdisperion

As described in the accompanying paper, we may be concerned with
infectious cell counts being overdispersed, i.e., following a negative
binomial distribution as opposed to a Poisson. When data are available
from multiple dilution levels, we can test the null hypothesis of no
overdispersion via a likelihood ratio test with the function
`lrt_SLDeepAssay_md()`.

The arguments for `lrt_SLDeepAssay_md()` are the same as those for
`fit_SLDeepAssay_md()`, along with three additional optional arguments
needed to maximize the negative binomial likelihood. The defaults for
these arguments tend to work well in most cases.

- `lb` Lower-bound on the IUPM (passed to ). Default is `lb = 1E-6`.
- `ub` Upper-bound on the IUPM (passed to ). Default is `ub = Inf`.
- `k0` initial value for the dispersion parameter k in optimization
  procedure. A value of k = `Inf` corresponds to no overdispersion.
  Default is `k = 1`.

``` r
lrt_res <- lrt_SLDeepAssay_md(assay_summary = assay_md)

lrt_res
```

    ## $mle
    ## [1] 0.4495329
    ## 
    ## $mle_bc
    ## [1] 0.4352599
    ## 
    ## $mle_negbin
    ## [1] 0.4495329
    ## 
    ## $mle_gamma
    ## [1] 0
    ## 
    ## $lrt_stat
    ## [1] 0

## Imperfect Assays

Here we simulate data assuming imperfect assay, by providing
sensitivities and specificities of the assays of 90%.

``` r
assay_imp <- simulate_assay_sd(M = 24,
                               tau = rep(0.5 / 6, 6),
                               q = 1,
                               sens_QVOA = 0.9,
                               spec_QVOA = 0.9,
                               sens_UDSA = 0.9,
                               spec_UDSA = 0.9)
```

Then we can analyze these data (i) naively assuming that the assays have
perfect sensitivity and (ii) assuming they have 90% sensitivity and
specificity.

The naive fit uses the original `fit_SLDeepAssay_sd()` function. Note
that the estimated IUPM and corresponding confidence interval are far
from the true IUPM of 0.5.

``` r
# naive fit (assuming perfect assays)
res_imp1 <- fit_SLDeepAssay_sd(assay = assay_imp$DVL_specific)

res_imp1
```

    ## $mle
    ## [1] 1.971558
    ## 
    ## $se
    ## [1] 0.3979641
    ## 
    ## $ci
    ## [1] 1.327371 2.928375
    ## 
    ## $mle_bc
    ## [1] 1.856905
    ## 
    ## $ci_bc
    ## [1] 1.220011 2.826283

The `fit_SLDeepAssay_sd_imperfect()` function estimates the IUPM under
given values for the sensitivities and specificities of the assays. Note
that this function, unlike `fit_SLDeepAssay_sd()` requires QVOA and UDSA
data, which are stored in `assay_imp$any_DVL` and
`assay_imp$DVL_specific`, respectively.

With the sensitivites and specificities correctly adjusted for, the
estimated IUPM is much closer to 0.5, and the 95% confidence interval
contains the true value.

``` r
# fit accounting for imperfect assays
res_imp2 <- fit_SLDeepAssay_sd_imperfect(assay_QVOA = assay_imp$any_DVL,
                                         assay_UDSA = assay_imp$DVL_specific,
                                         sens_QVOA = 0.9,
                                         spec_QVOA = 0.9,
                                         sens_UDSA = 0.9,
                                         spec_UDSA = 0.9)

res_imp2
```

    ## $mle
    ## [1] 0.6238332
    ## 
    ## $se
    ## [1] 0.2251968
    ## 
    ## $ci
    ## [1] 0.307464 1.265735
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
