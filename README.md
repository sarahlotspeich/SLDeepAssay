# Quantifying the HIV reservoir with dilution assays and deep viral sequencing

This repository contains `R` code and simulation data to reproduce results from the [manuscript](https://arxiv.org/abs/2209.04716) by Sarah C. Lotspeich, Brian D. Richardson, Pedro L. Baldoni, Kimberly P. Enders, and Michael G. Hudgens (2023+). 

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

  - [Script (Run Analysis)](real_data_application_revision/DATA-ANALYSIS-REVISION.R)
  - [Script (Make Table)](tables_revision/TableS5_revision.R)
  - [Data (Analysis Results)](real_data_application_revision/tableS5_data.csv)

**Table S6.** Empirical power of the likelihood ratio test for overdispersion

  - [Script (Run Simulations)](sims/SIMS-OVERDISPERSION.R)
  - [Script (Make Table)](tables_revision/TableS6_revision.R)
  - [Data (Analysis Results)](sim_data/overdispersion_sim_data.csv)

**Table S7.** Overdispersion LRT results for HIV Application

  - [Script (Run Test)](real_data_application_revision/DATA-ANALYSIS-REVISION.R)
  - [Script (Make Table)](tables_revision/TableS7_revision.R)
  - [Data (Test Results)](real_data_application_revision/tableS7_data.csv)

**Table S8.** Sensitivity Analysis for Subject C13

  - [Script (Run Sensitivity Analysis)](real_data_application_revision/DATA-ANALYSIS-REVISION.R)
  - [Script (Make Table)](tables_revision/TableS8_revision.R)
  - [Data (Test Results)](real_data_application_revision/tableS8_data.csv)

## Figures 

**Figure 1.** Illustration of the data collection scheme from the QVOA and UDSA at a single dilution level. (No code or data to share.)

**Figure 2.** Estimated infectious units per million (IUPM) with 95\% confidence intervals for 17 people living with HIV in the University of North Carolina HIV Cure Center Study. The IUPM and confidence interval were log transformed for comparisons of precision.

  - [Script (Run Analysis)](real_data_application_revision/DATA-ANALYSIS-REVISION.R)
  - [Script (Make Figure)](figures_revision/Figure2.R)
  - [Data (Analysis Results)](real_data_application_revision/figure2_data.csv)

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
