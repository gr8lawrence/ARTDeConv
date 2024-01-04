#' Example Data For Deconvolution
#' 
#' This example contains processed PBMC bulk expression and gene signature matrices for four major cell types
#' (T cell, B cell, monocyte, dendritic cell) from Hoek et al. (2015), as well as proportion mean and ranges 
#' based on data from Kleiveland et al. (2015). If you end up using the data, please cite the original authors.
#' 
#' @docType data
#' 
#' @usage deconv_ls
#' 
#' @format an object of class \code{"list"} that contains 4 fields: \code{bulk_mat}, \code{sig_mat}, \code{M} and
#'  \code{R}, corresponding to the bulk expression matrix, the signature matrix, the proportion mean values, and the proportion range values.
#' 
#' @keywords example data
"deconv_ls"