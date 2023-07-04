# SimCD
![GitHub Logo](/Miscel/Fig1_A4_cropped.png)
## Overview
**Sim**ultaneous **C**lustering and condition-specific **D**ifferential (**SimCD**) expression analysis is a unified Bayesian method based on a hierarchical gamma-negative binomial (hGNB) model, to simultaneously
perform clustering and differential expression analysis for single-cell RNA-seq data. SimCD is capable of including both gene- and cell-level biological explanatory variables
to better model scRNA-seq data and it also obviates the need for any sophisticated pre-processing steps.

## Installation
After loading the SimCD R and C codes available in the "Src" directory of this github repo, the only remaining installation step is to create shared object files (.so files) from the c codes so that they can be later loaded into R using "dyn.load" function. The shared object files can be created from the C source codes (.c files) using the below script:

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
Given a count gene expression matrix from single cell RNA-sequencing data as well as  cell and gene level biological covariates, SimCD simultaneously cluster cells and detect condition-specific differentially expressed genes for scRNA-seq count data. Here, in below R codes we have applied SimCD to a pre-processed CORTEX scRNA-seq data. The pre-processed CORTEX dataset gene expression data, cell-level covariates and gene-level covariates are in `count_cortex`, `xcov` and `zcov` parameters in the below R code resepectively. We have run the `scnbr_v4` function which is the main SimCD function to infer the latent parameters with latent space dimension (`K`) equal to 10. Here, in this example we have collected 1000 MCMC samples (`Collections`) after 1000 burn-in iterations (`Burnin`). 

``` r
library(doParallel)
library(foreach)
dyn.load('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_sum.so')
dyn.load('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_vector.so')
dyn.load('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_sum_matrix.so')
dyn.load('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_matrix.so')
dyn.load('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_MultR.so')
source('/scratch/user/naminiyakan/hgnbcode/LFC/dirrnd.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/KLsym.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/logcosh.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/logOnePlusExp.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/PolyaGamRnd_Gam.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/scnbr_v5.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_sum.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_sum_matrix.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_vector.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_matrix.R')
source('/scratch/user/naminiyakan/hgnbcode/LFC/CRT_MultR.R')
load("/scratch/user/naminiyakan/hgnbcode/LFC/Cortex_data.RData")
reshgnb<- scnbr_v4(counts=count_cortex, X=xcov, Z=zcov, K=10, ncores=24, Burnin = 1000L, Collections = 1000L, PGTruncation = 10L, randtry = 2020)
save.image("/scratch/user/naminiyakan/hgnbcode/LFC/res_cortex_K10.RData")

```
