#Helper functions for DESEQ2


#' Check formatting of supplied count matrix
#'
#' \code{check_count_matrix()} returns messages to the console to determine if
#' count matrix has been formatted correctly for use in downstream analysis
#'
#' The count matrix, normally as outputted from featureCounts, needs to be
#' formatted with:\itemize{
#' \item Column headings that match metadata
#' \item No rownames
#' \item First column titled gene_ensembl and contains all gene ensembl names
#' } This also provides a print out of the dimensions and predicted treatments
#' @export
check_count_matrix <- function(data = cts){
  #check column names match metadata

  #check there is no rownames
  if(tibble::has_rownames(data) == TRUE){
    message("Rownames detected but should not be used.
            Use tibble::rownames_to_column('gene_ensemble')
            to correct dataframe")
  } else {
  message("Rownames not detected: PASS")
  }

}
