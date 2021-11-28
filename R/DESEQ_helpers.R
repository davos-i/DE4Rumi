#Helper functions for DESEQ2


utils::globalVariables("where")


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
#'   that matches names of columns in \code{count_data}. This column hould not
#'    include and entry called "gene_ensembl", as this is automatically detected.
#'
#' @return Returns a dataframe of \code{count_data} with columns ordered to
#' match names in the \code{column_with_col_names} column of \code{colData}.
#'
#' @export

check_count_matrix <-
  function(count_data,
           colData,
           column_with_col_names) {


    column_as_string <- rlang::enquo(column_with_col_names) %>% rlang::as_label()
    #check column names match metadata
    #select column names (not include gene_ensembl column)
    count_data_names <- count_data %>%
      dplyr::select(-.data$gene_ensembl) %>%
      colnames()

    colData_names <- colData %>%
      dplyr::select({{ column_with_col_names }}) %>%
      magrittr::extract2(column_as_string)


    #1 do names ALL match (not necessarily in order)
    names_logical <- all(count_data_names %in% colData_names)
    #1b Do all colData names have a match in the counts data? If not, that's ok but good to list so user knows.
    colData_names_logical <- all(colData_names %in% count_data_names)
    colData_missing_names <- colData_names[which(!(colData_names %in% count_data_names))]

    if(names_logical == TRUE){
      message(crayon::black$bgGreen$bold("PASS: All columns in count_data have matching data in colData"))
      if(colData_names_logical == FALSE){
        message(crayon::red$bold("FAIL: colData has entries in column_with_col_names that do not match a column in count_data. Potentially missing data?"))
        message(crayon::red(paste("Entries in colData with no matching counts data:", paste(colData_missing_names, collapse = ", "))))
        stop("Checks terminated. See notes in Console. Run subset_colData() and try again.")
      }
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
      dplyr::select(.data$gene_ensembl, tidyselect::any_of(colData_names))

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



#' Subset colData automatically
#'
#' Run if [check_count_matrix] warns that there are entries in colData
#' that do not match columns in counts_matrix.
#'
#' Matches the entries in the \code{column_with_col_names} in \code{colData} to
#' the columns of \code{count_data}, and removes them from \code{colData}.
#'
#' @param count_data Dataframe of counts data to check.
#' @param colData Dataframe of column annotation data. each column should be
#'   metadata about the sample/animal. Particularly include all treatment info
#'   and animal ID in columns.
#' @param column_with_col_names non-string. Name of the column in \code{colData}
#'   that matches names of columns in \code{count_data}.
#'
#' @return Returns the filtered dataframe (to be used downstream).
#'
#' @export



subset_colData <-
function(count_data,
         colData,
         column_with_col_names) {

  column_as_string <- rlang::enquo(column_with_col_names) %>% rlang::as_label()
  #check column names match metadata
  #select column names (not include gene_ensembl column)
  count_data_names <- count_data %>%
    dplyr::select(-.data$gene_ensembl) %>%
    colnames()

  colData_names <- colData %>%
    dplyr::select({{ column_with_col_names }}) %>%
    magrittr::extract2(column_as_string)


  #1b Do all colData names have a match in the counts data? If not, that's ok but good to list so user knows.
  colData_names_logical <- all(colData_names %in% count_data_names)
  colData_missing_names <- colData_names[which(!(colData_names %in% count_data_names))]

  if(colData_names_logical == TRUE){
    message(crayon::green(("All entries in colData have matching count_data column. Run full check with check_counts_matrix() function. ")))
    return(colData)
  } else{
  message(crayon::red(paste("Entries in colData with no matching counts data:",
                            paste(colData_missing_names, collapse = ", "))))
    message(crayon::green(paste("Dimensions before filtering: ", paste(dim(colData), collapse = "   "))))
    colData_out <-
      colData %>%
      dplyr::filter({{ column_with_col_names }} %in% count_data_names )

    message(crayon::green(paste("Dimensions after filtering: ", paste(dim(colData_out), collapse = "   "))))
  }

  return(colData_out)
}




#'Collects human equivalent gene names and descriptions
#'
#'Creates a dataframe of human equivalent gene
#'names and description souced from gprofiler2::gconvert.
#'
#'The dataframe produced by this is used for downstream annotation of files. It
#'is also attached to the data when it is in the SummarizedExperiment class.
#'
#'@param data dataframe. Should be counts matrix with first column name as
#'  'gene_ensembl'
#'@param organism string. From gprofiler2: Organism names are constructed by
#'  concatenating the first letter of the name and the family. For sheep data:
#'  "oaries". For cattle data:"btaurus". Default is "oaries".
#'@param base_URL string. URL for the version of g:Profiler to use. This
#'  defaults to e100 as it matches the galaxy workflow, but also e103 broke a
#'  lot of names and is yet unresolved. Provided in the format:
#'  "http://biit.cs.ut.ee/gprofiler_archive3/e100_eg47_p14". \cr\cr Note that it is
#'  "http://" NOT "https://". \cr\cr Providing NA will default to g:Profiler2 default
#'  of "http://biit.cs.ut.ee/gprofiler". \cr\cr All available archives at:
#'  https://biit.cs.ut.ee/gprofiler/page/archives
#'
#'@export
#'
annotate_gene_ensembl <- function(data,
                                  organism = "oaries",
                                  base_URL = "http://biit.cs.ut.ee/gprofiler_archive3/e100_eg47_p14") {

  if(is.na(base_URL)){
    gprofiler2::set_base_url("http://biit.cs.ut.ee/gprofiler")
    message(paste("Using default g:Profiler URL: ", gprofiler2::get_base_url()))
  } else {
    gprofiler2::set_base_url(base_URL)
    message(paste("g:Profiler Version URL: ", gprofiler2::get_base_url()))
  }


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
    dplyr::select(gene_ensembl = .data$input,
                  gene_name = .data$name,
                  .data$description)

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
#'Creates a list of all possible unique
#'combinations of treatment comparisons
#'
#'This takes the colData from a dds object and finds the column name defined by
#'\code{contrast_factor} to generate all unique treatment names as input to a call
#'to [utils::combn] to get all unique pairwise comparisons.
#'Returns a list of characters of length 2. \cr
#'e.g. \code{[1] "LIV_HCP-HP-UMEI" "LIV_LCP-LP-UMEI"}
#'
#'@param coldata either a dataframe or a call to [SummarizedExperiment::colData]
#' on a DESeqDataSet
#'@param contrast_factor string. Name of column containing treatments to contrast
#'
#'
#'@export
#'
make_pairwise_combinations <-
  function(coldata,
           contrast_factor
          ){
    unique_contrast_levels <-
      coldata %>%
      as.data.frame() %>%
     # dplyr::filter({{ top_level_column_name }} == top_level_groups) %>%
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
#'(e.g. Tissue Region) and uses the [make_pairwise_combinations] function
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
#'@param alpha numeric. What is the p-value threshold to be used for
#'   determining significance. Used in call to [DESeq2::results] and
#'   others.
#'@param use_IHW_filtering Logical. Inherits from larger function. If \code{TRUE}
#'it will use a call to [IHW::ihw()] when making results.
#'
#'@return Returns a named list with each element the results object output from
#'call to [DESeq2::results]
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
          res <- DESeq2::results(dds_object,
                                 contrast = c(contrast_factor_string,
                                              contrast_numerator,
                                              contrast_denominator),
                                 alpha = a) %>%
            IHW::ihw(alpha = a)
          } else {
          res <- DESeq2::results(dds_object,
                                 contrast = c(contrast_factor_string,
                                              contrast_numerator,
                                              contrast_denominator),
                                 alpha = a)
        }

        message(crayon::green(paste("\n\n",res@elementMetadata$description[2])))

        DESeq2::summary(res, alpha = a)
        return(res)
      }


    res_out <-
      purrr::map(combinations, #iterate over all combinations
                        .f_make_results #using this internal function
      )

    #set names of the res_out list
    res_out <- stats::setNames(res_out,
                        lapply(combinations,
                               function(x) {paste(x, collapse = " vs ")})
    )


    #res_out is a list of all res objects for each comparison, named.
    return(res_out)
  }


#'Make SE object
#'
#'Takes gene expression count data, gene annotations and column data and returns
#'a SummarizedExperiment object to use with
#'\code{\link{auto_generate_DE_results}}
#'
#'\code{counts_data} and \code{colData} inputs to this function should be
#'pre-checked using \code{\link{check_count_matrix}}
#'
#'@param counts_data Dataframe of gene expression counts data.
#'@param gene_annotations Dataframe with columns "gene_ensembl", "gene_name" and
#'  "description", returns from a call to \code{\link{annotate_gene_ensembl}}
#'@param colData Dataframe of column annotation data. each column should be
#'  metadata about the sample/animal.
#'
#'@return Returns a SummarizedExperiment object of experimental data.
#'
#'
#'@export
#'
make_summarized_experiment_object <-
  function(counts_data,
           gene_annotations,
           colData){

    SummarizedExperiment::SummarizedExperiment(
      assays = counts_data %>%
        tibble::column_to_rownames(var = "gene_ensembl") %>%
        as.matrix(),
      rowData = gene_annotations,
      colData = colData)
  }





















