# SimCD
![GitHub Logo](/Miscel/Fig1_A4_cropped.png)
## Overview
**Sim**ultaneous **C**lustering and condition-specific **D**ifferential (**SimCD**) expression analysis is a unified Bayesian method based on a hierarchical gamma-negative binomial (hGNB) model, to simultaneously
perform clustering and differential expression analysis for single-cell RNA-seq data. SimCD is capable of including both gene- and cell-level biological explanatory variables
to better model scRNA-seq data and it also obviates the need for any sophisticated pre-processing steps.

## Installation
After loading the SimCD R and C codes available in the "Src" directory of this github repo, the only remaining installation step is to create shared object files (.so files) from the c codes so that they can be later loaded into R using "dyn.load" function. The shared object files can be created from the C source codes (.c files) using the below script:
```
R CMD SHLIB [options] [-o dllname] files
```

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
Given a count gene expression matrix from single cell RNA-sequencing data as well as  cell and gene level biological covariates, SimCD simultaneously cluster cells and detect condition-specific differentially expressed genes for scRNA-seq count data. Here, in below R codes we have applied SimCD to a pre-processed Src scRNA-seq data. The pre-processed CORTEX dataset gene expression data, cell-level covariates and gene-level covariates are in `count_cortex`, `xcov` and `zcov` parameters in the below R code resepectively. We have run the `scnbr_v4` function which is the main SimCD function to infer the latent parameters with latent space dimension (`K`) equal to 10. Here, in this example we have collected 1000 MCMC samples (`Collections`) after 1000 burn-in iterations (`Burnin`). 

``` r
library(doParallel)
library(foreach)
# Loading the shared object files created from .c files based on the steps described in the SimCD installation procedure:
dyn.load('/scratch/user/naminiyakan/SimCD/Src/CRT_sum.so')
dyn.load('/scratch/user/naminiyakan/SimCD/Src/CRT_vector.so')
dyn.load('/scratch/user/naminiyakan/SimCD/Src/CRT_sum_matrix.so')
dyn.load('/scratch/user/naminiyakan/SimCD/Src/CRT_matrix.so')
dyn.load('/scratch/user/naminiyakan/SimCD/Src/CRT_MultR.so')

source('/scratch/user/naminiyakan/SimCD/Src/dirrnd.R')
source('/scratch/user/naminiyakan/SimCD/Src/KLsym.R')
source('/scratch/user/naminiyakan/SimCD/Src/logcosh.R')
source('/scratch/user/naminiyakan/SimCD/Src/logOnePlusExp.R')
source('/scratch/user/naminiyakan/SimCD/Src/PolyaGamRnd_Gam.R')
source('/scratch/user/naminiyakan/SimCD/Src/scnbr_v5.R')
source('/scratch/user/naminiyakan/SimCD/Src/CRT_sum.R')
source('/scratch/user/naminiyakan/SimCD/Src/CRT_sum_matrix.R')
source('/scratch/user/naminiyakan/SimCD/Src/CRT_vector.R')
source('/scratch/user/naminiyakan/SimCD/Src/CRT_matrix.R')
source('/scratch/user/naminiyakan/SimCD/Src/CRT_MultR.R')

# Loading the pre-processed CORTEX dataset
load("/scratch/user/naminiyakan/SimCD/CORTEX/Cortex_data.RData")

# Running SimCD by passing the pre-processed CORTEX dataset gene expression data, cell-level covariates and gene-dependent covariates
# counts: scRNA-seq gene expression count data
# X: cell-level biological covariates such as treatment condition, age and sex
# Z: gene-level biological covariates such as gene length or GC content
# K: Latent space dimension
# Burnin: Number of burn-in iterations in MCMC sampling
# Collections:  Number of collected posterior samples after burn-in iterations
# PGTruncation: The truncation level used for genrating random numbers from Polya-Gamma distribution
# randtry: To be used for set.seed() function.
resSimCD<- scnbr_v4(counts=count_cortex, X=xcov, Z=zcov, K=10, ncores=24, Burnin = 1000L, Collections = 1000L, PGTruncation = 10L, randtry = 2020)

save.image("/scratch/user/naminiyakan/SimCD/CORTEX/res_cortex_K10.RData")

```

After that the Gibbs sampling inference procedure is done and SimCD's latent parameters are learned, we can visualize the cell embedding space by plotting the t-SNE visualization of latent space model parameter $\theta_j$ which can be later used for cell clustering. To do this, one can use the below R code:

``` r
library(Rtsne)
library(ggplot2)
# The SimCD-inferred embedding of cells can be accessed in reshgnb$Theta
tsne<-Rtsne(t(resSimCD$Theta),perplexity = 40,theta=0,pca=FALSE)

ggplot(as.data.frame(tsne$Y), aes((V1),(V2), color=factor(cell_type),show.legend =FALSE)) +
  labs(colour = "Cell type") +
  guides(color=FALSE) +
  geom_point(size=1.5) +
  xlab(paste0("Dim1")) +
  ylab(paste0("Dim2")) + 
  scale_colour_tableau() +
  theme_bw() + 
  coord_fixed()

```
![GitHub Logo](/Miscel/tsne_simcd_cortex_v2.png)

To do condition-specific differential expression (DE) analysis using SimCD, one can use the derived symmetric Kullback–Leibler (KL) divergence values for each gene based on inferred model parameters $\beta_{g}^{(1)}$ and rank them to extract the full list of DE genes across specific conditions. These values can be accessed from `resSimCD$kl` parameter.
