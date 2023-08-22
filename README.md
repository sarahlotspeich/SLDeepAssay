# Quantifying the HIV reservoir with dilution assays and deep viral sequencing

This repository contains `R` code and simulation data to reproduce results from the [manuscript](https://arxiv.org/abs/2209.04716) by Sarah C. Lotspeich, Brian D. Richardson, Pedro L. Baldoni, Kimberly P. Enders, and Michael G. Hudgens (2023+). 

For the `imputeCensRd` package, which implements the conditional mean imputation approaches from the paper, can be found at the end of this README. 

## Tables 

**Table 1.** Simulation results with a single dilution level, assuming a constant rate of infected cells for all distinct viral lineages. The true IUPM in all settings was $T = 1$.

  - [Script (Run Simulations)](sims/SIMS-SINGLE-DILUTION.R)
  - [Script (Make Table)](tables_revision/Table1_revision.R)
  - [Data (Simulation Results)](sim_data/sd_sim_data.csv)

**Table 2.** Simulation results with multiple dilution levels and a constant rate of infected cells for all distinct viral lineages. In all settings, the true IUPM was $T = 1$ and the proportions of positive wells that were deep sequenced at the three dilution levels were $\pmb{q} = (0, 0.5, 1)$.

  - [Script (Run Simulations)](sims/SIMS-MULTIPLE-DILUTIONS.R)
  - [Script (Make Table)](tables_revision/Table2_revision.R)
  - [Data (Simulation Results)](sim_data/md_sim_data.csv)

**Table S1.** Simulation results with a single dilution level and a non-constant rate of infected cells for all distinct viral lineages. The true overall IUPM in all settings was $T = 1$.

  - [Script (Run Simulations)](sims/SIMS-SINGLE-DILUTION.R)
  - [Script (Make Table)](tables_revision/TableS1_revision.R)
  - [Data (Simulation Results)](sim_data/sd_sim_data.csv)

**Table S2.** Simulation results with a single dilution level and a constant rate of infected cells for all distinct viral lineages. The true overall IUPM in all settings was $T = 0.5$.

  - [Script (Run Simulations)](sims/SIMS-SINGLE-DILUTION.R)
  - [Script (Make Table)](tables_revision/TableS2_revision.R)
  - [Data (Simulation Results)](sim_data/sd_sim_data.csv)

**Table S3.** Simulation results with multiple dilution levels and a non-constant rate of infected cells for all distinct viral lineages. In all settings, the true IUPM was $T = 1$ and the proportions of p24-positive wells that were deep-sequenced at the three dilution levels were $\pmb{q} = (0, 0.5, 1)$.

  - [Script (Run Simulations)](sims/SIMS-MULTIPLE-DILUTIONS.R)
  - [Script (Make Table)](tables_revision/TableS3_revision.R)
  - [Data (Simulation Results)](sim_data/md_sim_data.csv)

**Table S4.** Summary of assay data collected from 17 individuals receiving care at the University of North Carolina HIV Cure Center

  - [Script (Make Table)](tables_revision/TableS4_revision.R)
  - [Data (Summarized)](real_data_application_revision/tableS4_data.csv)

**Table S5.** Estimated infectious units per million (IUPM) of HIV (with 95\% confidence intervals) for 17 individuals receiving care at the University of North Carolina HIV Cure Center

  - [Script (Run Analysis)](real_data_application_revision/DATA-ANALYSIS_REVISION.R)
  - [Script (Make Table)](tables_revision/TableS5_revision.R)
  - [Data (Analysis Results)](sim_data/tableS5_data.csv)

## Figures 

**Figure 1.** Illustration of the data collection scheme from the QVOA and UDSA at a single dilution level. (No code or data to share.)

**Figure 2.** Estimated infectious units per million (IUPM) with 95\% confidence intervals for 17 people living with HIV in the University of North Carolina HIV Cure Center Study. The IUPM and confidence interval were log transformed for comparisons of precision.

  - [Script (Run Simulations)](Sim-Scripts/FigureS3-Extrapolation-Methods-Weibull.R)
  - [Script (Make Figure)](Figure-Scripts/FigureS3-Extrapolation-Methods-Weibull.R)
  - [Data (Simulation Results)](Figure-Data/data_FigureS3.csv)  

# `SLDeepAssay`: A package for maximum likelihood estimation from serial dilution assays with partial deep sequencing

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
dilution level using the `simulate_assay_sd()` function.

There are 3 required arguments,

- `M`: the total number of wells (a scalar),
- `tau`: the mean counts of infected cells per well with each DVL (a
  vector of with all elements $>0$)
- `q`: the proportion of p24-positive wells to be deep sequenced (a
  scalar between $0$ and $1$),

and 1 optional argument,

- `remove_undetected`: a logical indicator for whether DVL that are not
  detected in any wells should be deleted (`TRUE`) or not (`FALSE`).

The function returns a list containing two versions of the simulated
data, corresponding to the quantitative viral outgrowth assay (QVOA)
(`$any_DVL`) and Ultra Deep Sequencing Assay of the Outgrowth Virus
(UDSA) (`$DVL_specific`) measures discussed in the manuscript.

For demonstration, we simulate a single-dilution assay of `M = 12` total
wells (each with `u = 1` million cells in it), with `n = 6` underlying
DVL, and assume that `q = 0.5` (i.e., 50%) of p24-positive wells will
undergo deep sequencing. Lastly, if we want to simulate data where the
overall IUPM = 0.5, we set `tau = rep(0.5 / 6, 6)` (or any other way to
distribute a total mass of 0.5 across the 6 DVL).

We now demonstrate a few ways to simulate, depending on our choices for
the optional arguments.

``` r
assay = simulate_assay_sd(M = 12, 
                          tau = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = TRUE)
assay
```

    ## $any_DVL
    ##  [1] 1 1 1 0 0 0 0 0 0 1 1 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    0    1    0    0    0    0    0    0    0    NA    NA    NA
    ## [2,]    0    0    1    0    0    0    0    0    0    NA    NA    NA
    ## [3,]    1    0    0    0    0    0    0    0    0    NA    NA    NA
    ## [4,]    0    0    1    0    0    0    0    0    0    NA    NA    NA

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
                          tau = rep(0.5 / 6, 6), 
                          q = 0.5, 
                          remove_undetected = FALSE)
assay
```

    ## $any_DVL
    ##  [1] 1 0 0 0 0 0 0 0 0 0 0 1
    ## 
    ## $DVL_specific
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## [1,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [2,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [3,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [4,]    1    0    0    0    0    0    0    0    0     0     0    NA
    ## [5,]    0    0    0    0    0    0    0    0    0     0     0    NA
    ## [6,]    0    0    0    0    0    0    0    0    0     0     0    NA

Notice that `assay$any_DVL` is the same as before, which makes sense
since it’s the same random seed. However, `assay$DVL_specific` now has
all `n = 6` rows, even though rows 1–2 and 4–5 contain all 0s.

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
supplied, then the only other *required* argument is `u`, which is the
dilution level in millions of cells per well. Then, you can skip
parameters `M`, `MP`, `m`, and `Y` (they will be discussed in the next
section), and all that is left are optional arguments:

- `corrected`: a logical indicator for whether the bias-corrected MLE
  should be returned (`TRUE`) or not (`FALSE`),
- `maxit`: the maximum number of iterations allowed to find the MLE,
- `lb`: the lower bound for the parameters, and
- `ub`: the upper bound for the parameters.

The last three optional parameters are all passed to the `optim()`
function, which is used to find the MLE inside of `simulate_assay_sd()`.
Recall that data were simulated for wells with 1 million cells per well,
so `u = 1` and the MLE can be interpreted as the infectious units per
million (IUPM).

We now estimate the IUPM based on the last `assay` generated above,
leaving all optional arguments at their defaults.

``` r
res = fit_SLDeepAssay_sd(assay = assay$DVL_specific, 
                         u = 1)
res
```

    ## $mle
    ## [1] 0.182324
    ## 
    ## $se
    ## [1] 0.1307294
    ## 
    ## $ci
    ## [1] 0.04472169 0.74330946
    ## 
    ## $mle_bc
    ## [1] 0.5609867
    ## 
    ## $ci_bc
    ## [1] 0.3552982 0.8857518

Finally, we interpret the output as follows.

- `res$mle` is the uncorrected MLE.
- `res$se` is the estimated standard error (assumed to be the same for
  both the bias corrected and uncorrected MLE).
- `res$ci` is the 95% confidence interval (CI) for the uncorrected MLE.
- `res$mle_bc` is the bias-corrected MLE.
- `res$ci_bc` is the 95% CI for the bias-corrected MLE.

### Supplying summarized assay data
