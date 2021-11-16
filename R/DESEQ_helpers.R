#Helper functions for DESEQ2


#' Check formatting of supplied count matrix
#'
#' \code{check_count_matrix()} returns messages to the console to determine if
#' count matrix has been formatted correctly for use in downstream analysis, and
#'  to see if column names match colData (metadata for each column).
#'
#' The count matrix, normally as outputted from featureCounts, needs to be
#' formatted with:\itemize{ \item Column headings that match metadata \item No
#' rownames \item First column titled gene_ensembl and contains all gene ensembl
#' names }
#' This function also prints the dimensions and predicted treatments to
#' console.
#'
#' @param count_data Dataframe of counts data to check.
#' @param colData Dataframe of column annotation data. each column should be
#'   metadata about the sample/animal. Particularly include all treatment info
#'   and animal ID in columns.
#' @param column_with_col_names non-string. Name of the column in \code{colData}
#'   that matches names of columns in \code{count_data}. Should not include
#'   \code{ gene_ensembl} column as this is automatically detected.
#'
#' @return Returns a dataframe of \code{count_data} with columns ordered to
#' match names in the \code{column_with _col_names} column of \code{colData}.
#'
#' @export

check_count_matrix <- function(count_data,
                               colData,
                               column_with_col_names = sample_names) {
  #check column names match metadata
  #select column names (not include gene_ensembl column)
  count_data_names <- count_data %>%
    dplyr::select(-.data$gene_ensembl) %>%
    colnames()

  colData_names <- colData %>%
    dplyr::select({{ column_with_col_names }}) %>%
    magrittr::extract2(ensym(column_with_col_names))


  #1 do names ALL match (not necessarily in order)
  names_logical <- all(count_data_names %in% colData_names)

  if(names_logical == TRUE){
    message(crayon::black$bgGreen$bold("PASS: All columns in count_data have matching data in colData"))
  } else if(names_logical == FALSE){
    #2 if not, which ones are missing
     message("FAIL: Not all columns in count_data have matching data in colData")
    message("Missing samples:")
    print(count_data_names[which(!(count_data_names %in% colData_names))])
    message("Function will now output count_data without these columns!!")
  }

  #re-order count matrix columns to match coldata
  #(and only select columns with data)
  count_data_out <- count_data %>%
    dplyr::select(.data$gene_ensembl, colData_names)

  #check there is no rownames
  if (tibble::has_rownames(count_data) == TRUE) {
    message(
      paste(
        "FAIL: Rownames detected but should not be used.",
        "Use tibble::rownames_to_column('gene_ensemble') to correct dataframe."
      )
    )
  } else {
    message(crayon::black$bgGreen$bold("PASS: Rownames not detected"))
  }

  #Check if first column name is "gene_ensembl"
  col_name_first_col <- colnames(count_data)[1]
  if (col_name_first_col == "gene_ensembl") {
    message(crayon::black$bgGreen$bold("PASS: First column name is 'gene_ensembl'"))
  } else {
    message(
      paste0(
        "FAIL: First column name is '",
        col_name_first_col,
        "' Should be 'gene_ensembl'"
      )
    )
  }

  #output dimensions
  message(crayon::green("\n Dimensions of raw counts dataset BEFORE sorting:"))
  message(crayon::green(paste(
    "genes: ", dim(count_data)[1],
    "    samples: ", dim(count_data)[2] - 1)
    ))

  message(crayon::blue(paste("\n Number of samples in colData:"),
                       length(colData_names)))

  message(crayon::green("\n Dimensions of raw counts dataset AFTER sorting:"))
  message(crayon::green(paste(
    "genes: ", dim(count_data_out)[1],
    "    samples: ", dim(count_data_out)[2] - 1)
  ))
return(count_data_out)
}



#'Collects human equivalent gene names and descriptions
#'
#'\code{annotate_gene_ensembl} creates a dataframe of human equivalent gene
#'names and description souced from gprofiler2::gconvert.
#'
#'The dataframe produced by this is used for downstream annotation of files.
#'It is also attached to the data when it is in the SummarizedExperiment class.
#'
#'@param data dataframe. Should be counts matrix with first
#'column name as 'gene_ensembl'
#'@param organism string. From gprofiler2: Organism names are constructed by
#'concatenating the first letter of the name and the family. For sheep data:
#'"oaries". For cattle data:"btaurus". Default is "oaries".
#'
#'@export
#'
annotate_gene_ensembl <- function(data, organism = "oaries") {
  message(paste("Using organism:", organism))
  converted0 <- gprofiler2::gconvert(
    query = data$gene_ensembl,
    organism = organism,
    target = "ENSG",
    # or ENSG | AFFY_BOVINE
    filter_na = T,
    mthreshold = 1
  )

  converted_genes <-
    converted0 %>%
    dplyr::select(gene_ensembl = input, gene_name = name, description)

  message(
    crayon::blue(
      "Important - this produces duplicate gene names from annotation - however gene_ensembl remain unique.
      This is why it's important to use gene_ensembl for all tasks and annotate output at end."
    )
  )
  print(paste("duplicate genes in annotation: ",
        nrow(converted_genes[duplicated(converted_genes$gene_name), ])))

  converted_genes
}



#'Make pairwise combinations
#'
#'\code{make_pairwise_combinations()} creates a list of all possible unique
#'combinations of treatment comparisons
#'
#'This takes the colData from a dds object and finds the column name defined by
#'\code{contrast_factor} to generate all unique treatment names as input to a call
#'to \code{combn()} to get all unique pairwise comparisons.
#'Returns a list of characters of length 2. \cr
#'e.g. \code{[1] "LIV_HCP-HP-UMEI" "LIV_LCP-LP-UMEI"}
#'
#'@param coldata either a dataframe or a call to colData() on a DESeqDataSet
#'@param contrast_factor string. Name of column containing treatments to contrast
#'@param top_level_column_name non-string.
#'Name of column containing top level variable. Defaults to: \code{Region}.
#'@param top_level_filter string. Name of top level (e.g. region) to filter data
#'with.
#'
#'@export
#'
make_pairwise_combinations <-
  function(coldata,
           contrast_factor#,
           #top_level_column_name = Region,
          # top_level_filter
          ){
    unique_contrast_levels <-
      coldata %>%
      as.data.frame() %>%
     # dplyr::filter({{ top_level_column_name }} == top_level_filter) %>%
      magrittr::extract2(contrast_factor) %>%
      as.character() %>% #removes factors if using colData(dds) as input
      unique()

    print(paste("Contrast levels:",
                paste(unique_contrast_levels, collapse = ", ")
    )
    )

    utils::combn(unique_contrast_levels, 2, simplify = FALSE)  #this makes the combos
  }



#'Calculate DESEQ results
#'
#'Internal function to generate a list of DESEQ2 result objects for multiple
#'pairwise comparisons
#'
#'This takes a dds object, which is pre-filtered to only contain one top-level
#'(e.g. Tissue Region) and uses the \code{make_pairwise_combinations} function
#'to generate all possible pairwise combinations, then use this to make all res
#'objects. It is designed to be used within larger function to automate DESeq2.
#'
#'@param dds_object The dds object, after DESeq2 has been run (e.g. dds_wald).
#'  Normally with only 1 top-level filter (e.g. Tissue Region)
#'@param contrast_factor Non-string. Name of the column containing the factor
#'  that will be contrasted. Should match one of the factors in DESeq2 design
#'  (e.g. Region_Diet)
#'@param combinations Either NA (default) or a list with each element 2 strings
#'(e.g. "LIV_HCP-HP-UMEI" "LIV_LCP-LP-UMEI")
#'@param use_IHW_filtering Logical. Inherits from larger function. If \code{TRUE}
#'it will use a call to \code{IHW::ihw()} when making results.
#'
#'@return Returns a named list with each element the results object output from
#'call to \code{DESeq2::results()}
#'



.calculate_DESEQ_results <-
  function(dds_object,
           contrast_factor,  #Name of column containing treatments to contrast
           combinations = NA, #If NA, uses make_pairwise_combinations()
           alpha = 0.05,
           use_IHW_filtering){

    cf <- rlang::enquo(contrast_factor)

    #make combinations if not already provided by user
    if(is.na(combinations)){
      combinations <-
        make_pairwise_combinations(SummarizedExperiment::colData(dds_object),
                                   contrast_factor = rlang::as_label(cf) )

    }


    #function to iterate each of the pairwise combinations over
    .f_make_results <-
      function(x){
        contrast_numerator <- x[1]
        contrast_denominator <- x[2]
        contrast_factor_string <- rlang::as_label(cf)
        a <- alpha

        if (use_IHW_filtering == TRUE) {
          f_IHW <- function(results_object, alpha) {IHW::ihw(results_object, alpha = alpha)}
        }  else {
          f_IHW <- function(results_object){results_object}
        }


        res <- DESeq2::results(dds_object,
                               contrast = c(contrast_factor_string,
                                            contrast_numerator,
                                            contrast_denominator),
                               alpha = a) %>%   f_IHW(alpha = a)

        message(crayon::green(paste("\n\n",res@elementMetadata$description[2])))

        DESeq2::summary(res, alpha = a)
        return(res)
      }

    #run with output to a list
    # res_out <-
    #   pbapply::pblapply(combinations, #iterate over all combinations
    #                     .f_make_results #using this internal function
    #   )
    res_out <-
      purrr::map(combinations, #iterate over all combinations
                        .f_make_results #using this internal function
      )

    #set names of the res_out list
    res_out <- setNames(res_out,
                        lapply(combinations,
                               function(x) {paste(x, collapse = " vs ")})
    )


    #res_out is a list of all res objects for each comparison, named.
    return(res_out)
  }


#'Title
#'
#'\code{words_to_look_like_code} short description
#'
#'long description
#'
#'@param data dataframe.
#'@return
#'
#'
#'@export
