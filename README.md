# ARTdeConv: Adaptive Regularized Tri-factor non-negative matrix factorization for deConvolution

![Schematic Representation of The ARTdeConv Workflow](./images/ARTdeConv_schema.png)

## Overview


## Installation

To install the package, you can use the following code:
```R
require(devtools)
install_github("https://github.com/gr8lawrence/ARTDeConv", dependencies = TRUE)
```

## Quick Start

A schematic representation of the algorithm can be found at the top. To start using the algorithm, the ARTdeConv package includes an example dataset, which includes bulk expression from 8 human Peripheral Blood Mononuclear Cell (PBMC) samples from 2 subjects (HD30, HD31) across 4 time points (Day 0, 1, 3, 7). It also includes a matrix of gene signatures of 4 major PBMC cell types (T cell, B cell, monocyte, dendritic cell). A set of pre-calculated means and ranges for these cell types are also included. These data are based on the published studies of Hoek et al. (2015) and Kleiveland et al. (2015). 

After installing the package, the data can be viewed by
```R
library(ARTdeConv)
deconv_ls
```

To begin the deconvolution, we need to present the three ingredients shown in the schema, which are bulk expression(`deconv_ls$bulk_mat`), gene signature expression (`deconv_ls$bulk_mat`), and means and ranges of cell types (`deconv_ls$M` and `deconv_ls$R`). Following the analysis in the Section 3.2 of the ARTdeConv paper, we stipulate that there is one more "cell type" that encompasses all other cell types whose gene signatures are not measured, and we name it "others".

Using the notation of the, we have the following dimension parameters:

 * $m$: 73 (the number of gene signatures);
 * $n$: 8 (the number of samples);
 * $K$: 5 (the total number of cell types);
 * $K_0$: 4 (the number of cell types whose gene signature expression are known);
 
ARTdeConv requires the rows of the bulk matrix match those of the gene signature matrix (which is already satisfied in the processed data). It also requires the input signature matrix contains $K$ columns for all $K$ cell types, with the $K_0$ cell types with gene signature expression occupying the first $K_0$ columns, and the rest of the columns padded with 0. To meet this requirement, we pre-process the data:

```R
## extract the bulk and gene signature expression
Y = deconv_ls$bulk_mat
Theta = deconv_ls$bulk_mat
nrow(Y) == nrow(Theta_0) # verify the matching rows

## pad Theta with one column of 0 for the "others" cell type
Theta_0 = cbind(Theta, 0)

```

The first step is the cross-validation (CV), which can be used to determine the tuning parameters for the final deconvolution. To speed up the process, we use parallelization. One can disuse this by specifying `parallel = FALSE` in the function itself.
```R
2))

```


## More Resources


## Citation

To cite ARTdeConv, please use:


To cite the works of Hoek et al. and Kleiveland et al., please use:
