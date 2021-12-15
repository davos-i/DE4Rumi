#' Calculate RIF scores
#'
#' Calculates 2 types of Regulatory Impact Factors, based on Reverter work.
#'
#' RIF scores are ... Genes do not need to be DE to be able to have a regulatory
#' effect...(hence, all_genes_as_TF)
#'
#'
#'
#' @param DE_output object (class = list) returned from a call to
#'   [auto_generate_DE_results]
#' @param TFs character vector of Transcription Factors to use for RIF
#'   calculations. Optional and only used if \code{all_genes_as_TF = FALSE}
#' @param norm_exp_data Whole data normalisation output - use VST norm.
#' @param gene_annotations Output from call to [annotate_gene_ensembl]
#' @param all_genes_as_TF Logical. Should all genes (with row mean expression >
#'   2) be used as possible TFs in RIF?
#' @param colData coldata normally accessed via a call to
#'   [SummarizedExperiment::colData]
#' @param results_contrast_factor Name of Factor from
#'   \code{DESeq2_formula_design} that is used to generate pairwise
#'   comparisons.e.g. Region_Diet
#' @param samples_colname non-string. Name of \code{colData} column with
#'   sample names (used for filtering).
#'
#' @return Returns a list of dataframes, with each pairwise comparison having
#'   own set of RIF scores: RIF1, RIF2 and mean_RIF (mean of RIF1 and RIF2, for
#'   each gene)
#'
#' @export

calculate_RIF <- function(DE_output,
                           TFs,
                           norm_exp_data,
                           gene_annotations,
                           all_genes_as_TF = FALSE,
                           colData,
                           results_contrast_factor,
                           samples_colname){


  DE_data_list <- DE_output %>% purrr::map("DE_by_PIF_df")

  #Helpers
  col_of_sam <- rlang::enquo(samples_colname)
  samples_colname_string <- col_of_sam %>% rlang::as_label()
  cf <- rlang::enquo(results_contrast_factor)

  col_annot <- colData %>%
    as.data.frame() %>%
    dplyr::select(.data$sample_names, {{ cf }})


  # Define top level names
  top_level_names_DE_data_list <-
    names(DE_data_list) %>%
    stringr::str_extract(pattern = "^[^_]*")


  f_by_top_level <-
    function(.x, .y, cf, col_of_sam){

      list_of_DE_tables <- .x
      top_level_name <- .y

      ################################################################# #
      # Collect top_level_names from current iteration
      #
      # top_level_name <- names(list_of_DE_tables) %>%
      #   stringr::str_extract(pattern = "^[^_]*") %>%  #match 0 or more characters from beginning of string up to the first _
      #   unique()

      message(crayon::black$bgCyan$bold(paste("\n\n ******************* Start of - ",
                                              top_level_name,
                                              "******************* \n")))

      pairwise_names <- names(list_of_DE_tables) #%>%
      #stringr::str_split(pattern = "\\.", n = 2, simplify = TRUE) %>%
      #magrittr::extract(,2) %>% as.list()


      f_by_pairwise <-
        function(.x, .y, cf, col_of_sam, top_level_name){

          DE_table <- .x
          pairwise_comparison <- .y

          current_condition1vscondition2 <- stringr::str_split(pairwise_comparison, pattern = " vs ", simplify = T)

          ################################################################ #
          #Filter the whole coldata table to Conditions in the current current_condition1vscondition2
          #The %in% searches the whole dataframe (current_condition1vscondition2)
          #
          coldata_sub <-
            col_annot %>%
            dplyr::filter(stringr::str_detect(.data$sample_names,  top_level_name)) %>%
            dplyr::filter(!!cf  %in% current_condition1vscondition2) #%>% arrange()

          ################################################################ #
          #Filter normalised expression table for current iteration
          #
          columns_to_select <- coldata_sub[[samples_colname_string]]

          norm_exp_sub <-
            norm_exp_data %>%
            dplyr::select(.data$gene_ensembl, tidyselect::all_of(columns_to_select)) %>%
            dplyr::distinct(.data$gene_ensembl, .keep_all = T) %>%
            tibble::column_to_rownames("gene_ensembl") %>%
            as.matrix()

          ################################################################# #
          # Select DE genes as target for RIF
          # The DE_table parsed down to this function is already filtered for sig
          # Only need to select names

          Target <- DE_table$gene_ensembl

          ################################################################# #
          # All_genes_as_TF is used for basically all genes, excluding very
          # low mean counts.
          #
          if(all_genes_as_TF == TRUE) {
            TFs0 <- names(which(rowMeans(norm_exp_sub) > 2))
          } else if(all_genes_as_TF == FALSE) {
            if(is.na(TFs)){
              stop("No TFs provided and all_genes_as_TF == FALSE")
            } else{
              TFs0 <- TFs #This is normally a downloaded list of known TFs
            }
          }
          ################################################################# #


          if(length(Target) > 1){
            #target_list_output[[paste(region0,diet_compare0,Sys.time(),sep = "_")]] <- list(Target_genes = target_genes, TF = TFs0)

            ################################################################# #
            # Verify TF and Target genes
            #
            # Verifying which TFs and Target genes are in the subsetted normalized
            # data (as these lists come from databases, not all genes will be present in this dataset)
            #

            TFs0 <- rownames(norm_exp_sub)[rownames(norm_exp_sub) %in% TFs0]
            Target <- rownames(norm_exp_sub)[rownames(norm_exp_sub) %in% Target]

            # print details to console
            message(crayon::yellow("\n", paste0(top_level_name, " - ", pairwise_comparison),
                                   "\n", "Number of target genes: ", length(Target),
                                   "\n", "Number of TFs (or genes as TFS): ", length(TFs0)))

            ################################################################# #
            ## Order rows of normalized count data
            #
            RIF_input <- rbind(norm_exp_sub[Target,], norm_exp_sub[TFs0,])

            ################################################################# #
            # Remove rows of SD = 0
            # 1.Split dataset into 2 groups (condition 1 and condition 2)
            # get sample names for a and b
            # use sample names to split RIF_input into condition a and condition b

            condition1_names <-
              coldata_sub %>%
              dplyr::filter(!!cf  %in% current_condition1vscondition2[1]) %>%
              dplyr::pull(!!col_of_sam)

            condition2_names <-
              coldata_sub %>%
              dplyr::filter(!!cf  %in% current_condition1vscondition2[2]) %>%
              dplyr::pull(!!col_of_sam)


            list_separate_treatments <- list("a" = RIF_input[,which(colnames(norm_exp_sub) %in% condition1_names)],
                                             "b" = RIF_input[,which(colnames(norm_exp_sub) %in% condition2_names)])

            message(crayon::cyan(" names of condition 1: ",
                                 paste(colnames(list_separate_treatments$a), collapse = ", "),
                                 "\n", "names of condition 2: ", paste(colnames(list_separate_treatments$b),collapse = ", ")))


            ################################################################## #
            # 2. Identify genes with a SD of 0 within a condition
            #
            row_std <-  lapply(list_separate_treatments,
                               function(x){which(apply(x,1, stats::sd) == 0)})

            exclude <- unique(c(names(row_std$a), names(row_std$b)))

            #message(crayon::yellow("excluded genes (within condition SD == 0): ", paste(exclude, collapse = ", ")))
            message(crayon::yellow("Number of excluded genes (within condition SD == 0): ", length(exclude)))

            RIF_input <- RIF_input[!rownames(RIF_input) %in% exclude,]

            ################################################################# #
            # ############################ Run RIF  ######################### #
            #
            # Performing RIF analysis - outputs z-normalised RIF automatically

            n1 <- length(which(colnames(norm_exp_sub) %in% condition1_names))
            n2 <- length(which(colnames(norm_exp_sub) %in% condition2_names))


            if (!requireNamespace("CeTF", quietly = TRUE)) {
              stop("Package \"CeTF\" needed for RIF function to work. Please install it using BiocManager::install(\"CeTF\").",
                   call. = FALSE)
            }
            message(crayon::green("Running RIF... "))
            RIF_out <- CeTF::RIF(input = RIF_input,
                                 nta = length(Target),
                                 ntf = length(TFs0),
                                 nSamples1 = n1,
                                 nSamples2 = n2
            )

            #Calculate mean RIF scores, absolute of mean, and arrange
            RIF_out <- RIF_out %>%
              dplyr::rowwise() %>%
              dplyr::mutate(mean_RIF = mean(c(.data$RIF1,.data$RIF2), na.rm = TRUE)) %>%
              dplyr::ungroup() %>%
              dplyr::mutate(abs_mean_RIF = abs(.data$mean_RIF)) %>%
              dplyr::arrange(dplyr::desc(.data$abs_mean_RIF))

            #Annotate results
            RIF_out <- RIF_out %>%
              dplyr::left_join(gene_annotations, by = c("TF" = "gene_ensembl"))

            message(crayon::green("Running RIF... COMPLETE"))

          } else{
            message(crayon::magenta("\n",
                                    stringr::str_c(top_level_name, pairwise_comparison , sep = " - "),
                                    "\n", "Did not run RIF -- Number of target genes: ", length(Target)))
            #return nothing
            RIF_out <- NULL}
          return(RIF_out)
        }

      out <- purrr::map2(.x = list_of_DE_tables,
                         .y = pairwise_names,
                         f_by_pairwise,
                         cf = cf,
                         col_of_sam = col_of_sam,
                         top_level_name)

      ## remove NULL entries here
      out <- out %>% purrr::discard(is.null)

      message(crayon::black$bgCyan$bold(paste("\n\n ******************* End of - ",
                                              top_level_name,
                                              "******************* \n\n")))
      return(out)
    }
  out2 <- purrr::map2(.x = DE_data_list,
                      .y = top_level_names_DE_data_list,
                      .f = f_by_top_level,
                      cf = cf,
                      col_of_sam = col_of_sam)

  ###################################### #
  # Rename output
  # top_names <- names(DE_data_list) %>%
  #stringr::str_extract(pattern = "^[^_]*")
  new_names <- paste0(top_level_names_DE_data_list, " - RIF")
  out2 <- setNames(out2, new_names)

  return(out2)
}
