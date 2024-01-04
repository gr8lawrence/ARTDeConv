# ARTdeConv: Adaptive Regularized Tri-factor non-negative matrix factorization for deConvolution

![Schematic Representation of The ARTdeConv Workflow](./images/ARTdeConv_schema.png)

## Overview

ARTdeConv enables the deconvolution of bulk tissue gene expression using only the gene signatures of a subset of cell types. This capability is particularly valuable for scenarios where some cell types are not easily preserved for single-cell sequencing, resulting in unavailable gene signatures. It requires additional parameters such as means and ranges of cell-type proportions to facilitate accurate deconvolution results. This document offers a quick guide to using the package and outlines the prerequisites for the algorithm to function properly.

If you have any questions, feel free to create an Issue on the GitHub page or contact tianyi96@live.unc.edu.

## Installation

To install the package, use the following code:
```R
require(devtools)
install_github("https://github.com/gr8lawrence/ARTDeConv", dependencies = TRUE)
```

## Quick Start

At the top, a schematic of the algorithm is available. To begin using the ARTdeConv package, an example dataset is provided, comprising bulk expression data from 8 human Peripheral Blood Mononuclear Cell (PBMC) samples from 2 subjects (HD30, HD31) across 4 time points (Day 0, 1, 3, 7). The dataset also contains a matrix of gene signatures for 4 major PBMC cell types (T cell, B cell, monocyte, dendritic cell), along with pre-calculated means and ranges for these cell types. These data are sourced or derived from the published studies of Hoek et al. (2015) and Kleiveland et al. (2015).

After installing the package, the attached data can be viewed (as a list object)
```R
library(ARTdeConv)
deconv_ls
```

To start the deconvolution, we require three elements from the schema: bulk expression (`deconv_ls$bulk_mat`), gene signature expression (`deconv_ls$bulk_mat`), and means and ranges of cell types (`deconv_ls$M` and `deconv_ls$R`). Following the analysis in Section 3.2 of the ARTdeConv paper, we introduce an additional "cell type", named "others", which encompasses all unmeasured cell types.

Using the notation of the, we have the following dimension parameters:

 * $m$: 73 (the number of gene signatures);
 * $n$: 8 (the number of samples);
 * $K$: 5 (the total number of cell types);
 * $K_0$: 4 (the number of cell types whose gene signature expression are known);
 
ARTdeConv requires that the rows of the bulk matrix correspond to those of the gene signature matrix in terms of gene features (a condition already met in the attached processed data, as shown below): 

```R
## extract the bulk and gene signature expression
Y = deconv_ls$bulk_mat
Theta = deconv_ls$bulk_mat
nrow(Y) == nrow(Theta)
all.equal(rownames(Y), rownames(Theta)) # verify the matching rows
```

> Note: both `Y` and `Theta` have to be **matrix objects** in R. ARTdeConv currently does not support other formats of the bulk and signature matrices such as `ExpressionSet`.

It also requires the input signature matrix to have $K$ columns for all $K$ cell types, with the first $K_0$ columns representing the $K_0$ cell types with known gene signature expressions and the remaining columns padded with zeros. To fulfill this criterion, we preprocess the data:

```R
## pad Theta with one column of 0 for the "others" cell type
Theta_0 = cbind(Theta, 0)
colnames(Theta_0) = c(colnames(Theta), "others") 
```

> Note: $K_0 < K$ is a hard requirement for ARTdeConv. $K$ can be determined through expert knowledge on the tissue or as $K_0 + 1$ if the cell types with missing reference are lumped together, such as in the example here. If $K_0 = K$ is the case for the application, we recommend using a reference-based deconvolution approach instead.

The initial step involves cross-validation (CV), which helps identify the tuning parameters for the ultimate deconvolution. To expedite this process, we implement parallelization, which necessitates core registration. Alternatively, setting `parallel = FALSE` in the function (also applicable to `artdeconv()` in the subsequent code chunk) allows the process to run sequentially. 

The grid for tuning parameters and the number of folds also adhere to the example in Section 3.2 of the paper. The code for cross-validation is as follows:

```R
library(foreach)
library(doParallel)

cl = parallel::makeCluster(4)
registerDoParallel(cl)
cv_params = cv_artdeconv(Y = Y, 
                         Theta_0 = Theta_0,
                         k0 = ncol(Theta), 
                         meds = deconv_ls$M, 
                         ranges = deconv_ls$R, 
                         alpha1_range = 2^seq(-4, 0, 0.2), 
                         alpha2_range = 1e-12, 
                         beta_range = 2^seq(0, 4, 0.2),
                         tol = 1e-4,
                         n_fold = 4) 

```
After cross-validation, we can extract the optimal tuning parameters and use them for a single ARTdeConv run:

```R
best_fit = artdeconv(Y = Y, 
                     Theta_0 = Theta_0,
                     k0 = ncol(Theta), 
                     meds = deconv_ls$M, 
                     ranges = deconv_ls$R, 
                     alpha1 = cv_params$alpha1,
                     alpha2 = cv_params$alpha2,
                     beta = cv_params$beta,
                     tol = 1e-4)
```

Finally, we can obtain the deconvoluted proportions in a $K \times n$ matrix, with each row corresponding to each column of `Theta_0` (in that order) and each column to that of `Y`.

```R
P_hat = best_fit$P_hat
rownames(P_hat) = colnames(Theta_0)
colnames(P_hat) = colnames(Y)
P_hat
```

## More Resources

More vignettes are under development.

## Citation

To cite ARTdeConv, please use:

Liu, Tianyi, Quefeng Li, Xiaojing Zheng, and Fei Zou. 2023. “Adaptive Regularized Tri-Factor Non-Negative Matrix Factorization for Cell Type Deconvolution.” *bioRxiv*. https://doi.org/10.1101/2023.12.07.570631.

To cite the works of Hoek et al. and Kleiveland et al., please use:

Hoek, Kristen L., Parimal Samir, Leigh M. Howard, Xinnan Niu, Nripesh Prasad, Allison Galassie, Qi Liu, et al. 2015. “A Cell-Based Systems Biology Assessment of Human Blood to Monitor Immune Responses after Influenza Vaccination.” *PloS One* 10 (2): e0118528.

Kleiveland, Charlotte R. 2015. *Peripheral Blood Mononuclear Cells*. Springer.
