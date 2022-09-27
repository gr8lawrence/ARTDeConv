#' Returning a Delta Matrix 
#' 
#' Get a \eqn{M \times K} binary matrix (Delta matrix) where entries that are 1 correspond to 
#' known expression values in the signature matrix and 0 to unknown values; 
#' **we always assume the missing cell types are in the last columns**
#' 
#' @param k the total number of cell types (CTs) in the deconvolution
#' @param k0 the total number of known CTs in the deconvolution (\eqn{\leq K})
#' @param m the total number of genes in the deconvolution
#' @param m0 the total number of genes of the known CTs in the deconvolution (\eqn{\leq m})
#' 
#' @return a \eqn{M \times K} binary matrix
get_Delta <- function(k, k0, m, m0) {
  Delta = rbind(
    matrix(rep(c(rep(1, k0), rep(0, k - k0)), m0), m0, k, byrow=TRUE),
    matrix(0, m - m0, k)
  )
  return(Delta)
}

#' Function for Calculating the Frobenius Inner Products
#' 
#' Return the Frobenius (trace) inner products of two matrices of
#' the same size; will return an error message if the size does
#' not match
#' 
#' @param A a matrix
#' @param B a matrix of the same size of `A`
#' 
#' @return the Frobenius inner products of `A` and `B`
inner_product_TL <- function(A, B) {
  if (prod(dim(A) == dim(B)) == 0) stop("Dimensions of two matrices must be the same")
  inner = sum(A * B)
  return(inner)
}
