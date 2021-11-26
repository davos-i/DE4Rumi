#' coldata for testing
#'
#' A dataset of coldata from sheep experiment
#'
#' @format A data frame with 239 rows and 11 variables:
#' \describe{
#'   \item{sample_names}{Names of samples in correct format}
#'   \item{Region_Diet}{Treatment information}
#'   ...
#' }
#' @source private
"coldata"


#' Counts matrix for testing
#'
#' A counts matrix from sheep experiment, generated via feature counts.
#'
#' @format A data frame with 27054 rows and 240 variables:
#' \describe{
#'   \item{gene_ensembl}{Ensembl IDs}
#'   \item{Samples}{Each column is an individual sample, names match sample_names of coldata}
#'   ...
#' }
#' @source private
"cts_all"
