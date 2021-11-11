#Helper functions for DESEQ2


#' Check formatting of supplied count matrix
#'
#' \code{check_count_matrix()} returns messages to the console to determine if count
#' matrix has been formatted correctly for use in downstream analysis
#'
#' The count matrix, normally as outputted from featureCounts, needs to be
#' formatted with:\itemize{
#' \item Column headings that match metadata
#' \item No rownames
#' \item First column titled gene_ensembl and contains all gene ensembl names
#' } This also provides a print out of the dimensions and predicted treatments
check_count_matrix <- function(data = cts){
  #check column names match metadata

  #check there is no rownames
  message(paste("Are rownames present?", tibble::has_rownames()))

}
