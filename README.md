# SimCD
![GitHub Logo](/Miscel/Fig1_A4_cropped.png)
## Overview
**Sim**ultaneous **C**lustering and condition-specific **D**ifferential (**SimCD**) expression analysis is a unified Bayesian method based on a hierarchical gamma-negative binomial (hGNB) model, to simultaneously
perform clustering and differential expression analysis for single-cell RNA-seq data. SimCD is capable of including both gene- and cell-level biological explanatory variables
to better model scRNA-seq data and it also obviates the need for any sophisticated pre-processing steps.

## Quick Start
The main function scnbr_V4 takes as input:

- **Counts**: A matrix of scRNA-seq count data, rows are corresponding to genes and columns are corresponding to samples.
- **X**: A design matrix for cell level covariates such as condition, age and sex (i.e. $x_{j}^{(1)}$).
- **Z**: A design matrix for gene level covariates such as gene length or GC content (i.e. $x_{g}^{(2)}$).
- **K**: Latent space dimension.
- **Burnin**: Number of burn-in iterations in MCMC.
- **Collections**: Number of collected posterior samples after burn-in iterations.
- **PGTruncation**: The truncation level used for genrating random numbers from Polya-Gamma distribution
- **randtry**: To be used for set.seed() function.

## Example
