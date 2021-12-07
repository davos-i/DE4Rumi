# PCIT related functions

# 1. Top DE genes
# 2. Top PIF genes
# 3. Top RIF genes
# 4. Tissue Specific genes

# Input:
# 1. DE_out_object (need Raw PIF and DE_by_PIF ?)
# 2. RIF_output (can be run with either only TFs or all genes as TF)
# 3. Normalised counts (whole dataset)
# 4. Regions to include in dataset e.g. "LIV"

# Output (export to file):
# 1. PCIT result for Cytoscape
# 2. Venn plot
# 2. Attributes table for Cytoscape


#' Prepare PIF data for PCIT
#'
#' Internal function to prepare PIF data for PCIT
#'
#' @param DE_full_out Full output of [auto_generate_DE_results]
#' @param search_top_level character vector of top_level_groups (see
#'   [auto_generate_DE_results]) to search for, returning only tables from these
#'   groups. e.g. \code{c("ARC","LIV")}. If not provided it will return all DE
#'   tables.
#' @param PIF_sig Number. Threshold used for filtering the z-scored PIF values, so a
#'   \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas a
#'   value of 1.96 filters for a nominal p value of 0.05.
#' @param PIFn Number of values to return for both the max PIF and min PIF. E.g.
#' if \code{PIFn = 20}, then the 20 highest and 20 lowest are returned (if that many are present after filtering by PIF_sig)
#'
#' @return Returns a dataframe of all selected PIF genes,
#' with summary columns showing why they were selected.
#'
#'
.prep_PCIT_PIF <-
  function(DE_full_out,
           search_top_level = NULL,
           PIF_sig,
           PIFn){
    #################################################################### #
    # 1. Get data and filter by top level (group/region)
    raw_PIF_data <- GET_PIF_df(auto_DE_output = DE_full_out,
                               search_top_level )


    #################################################################### #
    # 2. Filter by PIF_sig and select top/bottom PIFn

    # Get names
    top_level_names_list <-
      names(raw_PIF_data) %>%
      stringr::str_extract(pattern = "^[^_]*") %>% as.list()

    pairwise_names <- names(raw_PIF_data) %>%
      stringr::str_split(pattern = " - ") %>%
      purrr::map(function(x){ x[[2]] })

    #iterate over all 3 (like map2 but for 3)
    list_in <-
      list(data_in = raw_PIF_data,
           top_name = top_level_names_list,
           pair_name = pairwise_names)

    filtered_PIF_data <-
      purrr::pmap(list_in,
                  function(data_in, top_name, pair_name){
                    #get column name of zPIF column
                    zPIF_name <- names(data_in)[stringr::str_detect(names(data_in), "zPIF")]

                    #filter, using absolute zPIF
                    x_filtered <-
                      data_in %>%
                      dplyr::filter(abs(.data[[zPIF_name]]) > PIF_sig) %>%
                      dplyr::select(.data$gene_ensembl,
                                    .data$gene_name,
                                    tidyselect::starts_with("zPIF"))

                    max_PIF <-
                      x_filtered %>% dplyr::slice_max(order_by = .data[[zPIF_name]],
                                                      n = PIFn)

                    min_PIF <-
                      x_filtered %>% dplyr::slice_min(order_by = .data[[zPIF_name]],
                                                      n = PIFn)

                    PIF_filt_out <-
                      dplyr::bind_rows(max_PIF, min_PIF) %>%
                      #rename zPIF.... column with "zPIF"
                      dplyr::rename_with(function(x){
                        stringr::str_trunc(x,4,"right", "")},
                        tidyselect::starts_with("zPIF")
                      ) %>%
                      dplyr::mutate(top_level = top_name,
                                    pairwise_comparison = pair_name)
                    PIF_filt_out
                  })

    #################################################################### #
    # 4. Create summary data
    # Collapse together into 1 dataframe to run summary,
    # Need gene_ensembl, gene_name, Top_levels, top_levels_n, Pairwise_comparisons
    #reduce
    PIF_df <- filtered_PIF_data %>%
      purrr::reduce(dplyr::full_join,
                    by = c("gene_ensembl",
                           "gene_name",
                           "zPIF",
                           "top_level",
                           "pairwise_comparison"))


    #Summarise
    PIF_filtered_summarised <-
      PIF_df %>%
      dplyr::group_by(.data$gene_ensembl, .data$gene_name) %>%
      dplyr::summarize(top_groups_included = paste(sort(unique(.data$top_level)),collapse=", "),
                       number_top_groups = dplyr::n_distinct(.data$top_level),
                       comparisons_included = paste(sort(unique(paste(.data$top_level,
                                                                      .data$pairwise_comparison,
                                                                      sep = "__"))),
                                                    collapse = ", ")) %>%
      dplyr::arrange(dplyr::desc(.data$number_top_groups))

    PIF_filtered_summarised
  }


#' Prepare PIF data for PCIT
#'
#' Internal function to prepare PIF data for PCIT
#' @param DE_full_out Full output of [auto_generate_DE_results]
#' @param search_top_level character vector of top_level_groups (see
#'   [auto_generate_DE_results]) to search for, returning only tables from these
#'   groups. e.g. \code{c("ARC","LIV")}. If not provided it will return all DE
#'   tables.
#' @param DEn Number of values to return, ranked on adjusted p value. E.g.
#' if \code{DEn = 20}, then the 20 genes with the lowest p value are returned
#' (there may be less than this due to prior filtering for significance).
#'
#' @return Returns a dataframe of all selected DE genes,
#' with summary columns showing why they were selected.


.prep_PCIT_DE <-
  function(DE_full_out,
           search_top_level = NULL,
           DEn){
    #################################################################### #
    # 1. Get data and filter by top level (group/region)
    raw_DE_data <- GET_DE_by_PIF_df(auto_DE_output = DE_full_out,
                               search_top_level )


    #################################################################### #
    # 2. Select min DEn, ranked on padj

    # Get names
    top_level_names_list <-
      names(raw_DE_data) %>%
      stringr::str_extract(pattern = "^[^_]*") %>% as.list()

    pairwise_names <- names(raw_DE_data) %>%
      stringr::str_split(pattern = "\\.", n = 2) %>%
      purrr::map(function(x){ x[[2]] })

    #iterate over all 3 (like map2 but for 3)
    list_in <-
      list(data_in = raw_DE_data,
           top_name = top_level_names_list,
           pair_name = pairwise_names)

    filtered_DE_data <-
      purrr::pmap(list_in,
                  function(data_in, top_name, pair_name){

                    min_DE <-
                      data_in %>% dplyr::slice_min(order_by = .data$padj,
                                                      n = DEn)

                    DE_filt_out <-
                      min_DE %>%
                      dplyr::mutate(top_level = top_name,
                                    pairwise_comparison = pair_name) %>%
                      dplyr::select(.data$gene_ensembl,
                                    .data$gene_name,
                                    .data$top_level,
                                    .data$pairwise_comparison)
                    DE_filt_out
                  })

        #################################################################### #
    # 4. Create summary data
    # Collapse together into 1 dataframe to run summary,
    # Need gene_ensembl, gene_name, Top_levels, top_levels_n, Pairwise_comparisons
    #reduce
    DE_df <- filtered_DE_data %>%
      purrr::reduce(dplyr::full_join,
                    by = c("gene_ensembl",
                           "gene_name",
                           "top_level",
                           "pairwise_comparison"))


    #Summarise
    DE_filtered_summarised <-
      DE_df %>%
      dplyr::group_by(.data$gene_ensembl, .data$gene_name) %>%
      dplyr::summarize(top_groups_included = paste(sort(unique(.data$top_level)),collapse=", "),
                       number_top_groups = dplyr::n_distinct(.data$top_level),
                       comparisons_included = paste(sort(unique(paste(.data$top_level,
                                                                      .data$pairwise_comparison,
                                                                      sep = "__"))),
                                                    collapse = ", ")) %>%
      dplyr::arrange(dplyr::desc(.data$number_top_groups))

    DE_filtered_summarised
  }


#' Prepare RIF data for PCIT
#'
#' Internal function to prepare RIF data for PCIT
#' @param RIF_output_data Full output of [calculate_RIF]
#' @param search_top_level character vector of top_level_groups (see
#'   [auto_generate_DE_results]) to search for, returning only tables from these
#'   groups. e.g. \code{c("ARC","LIV")}. If not provided it will return all DE
#'   tables.
#' @param RIF_sig Number. Threshold used for filtering the z-scored PIF values, so a
#'   \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas a
#'   value of 1.96 filters for a nominal p value of 0.05.
#' @param RIFn Number of values to return for both the max PIF and min PIF. E.g.
#' if \code{PIFn = 20}, then the 20 highest and 20 lowest are returned (if that many are present after filtering by PIF_sig)
#'
#' @return Returns a dataframe of all selected PIF genes,
#' with summary columns showing why they were selected.

.prep_PCIT_RIF <-
  function(RIF_output_data,
           search_top_level = NULL,
           RIF_sig,
           RIFn){
    #################################################################### #
    # 1. Get data and filter by top level (group/region)
    raw_RIF_data <- .search_list_helper(RIF_output_data,
                        search_term = search_top_level,
                        pattern_extract_top_level = "^[^ - ]*" ) %>%
      unlist(recursive = FALSE)


    #################################################################### #
    # 2. Filter by RIF_sig and select top/bottom RIFn

    # Get names
    top_level_names_list <-
      names(raw_RIF_data) %>%
      stringr::str_extract(pattern = "^[^ - ]*") %>% as.list()

    pairwise_names <- names(raw_RIF_data) %>%
      stringr::str_split(pattern = " - ") %>%
      purrr::map(function(x){ x[[2]] })

    #iterate over all 3 (like map2 but for 3)
    list_in <-
      list(data_in = raw_RIF_data,
           top_name = top_level_names_list,
           pair_name = pairwise_names)

    filtered_RIF_data <-
      purrr::pmap(list_in,
                  function(data_in, top_name, pair_name){
                    #get column name of zRIF column
                    #zRIF_name <- names(data_in)[stringr::str_detect(names(data_in), "mean_RIF")]
                    zRIF_name <- "mean_RIF"

                    #filter, using absolute zRIF
                    x_filtered <-
                      data_in %>%
                      dplyr::filter(abs(.data[[zRIF_name]]) > RIF_sig) %>%
                      dplyr::select("gene_ensembl" = .data$TF,
                                    .data$gene_name,
                                    .data[[zRIF_name]])

                    max_RIF <-
                      x_filtered %>% dplyr::slice_max(order_by = .data[[zRIF_name]],
                                                      n = RIFn)

                    min_RIF <-
                      x_filtered %>% dplyr::slice_min(order_by = .data[[zRIF_name]],
                                                      n = RIFn)

                    RIF_filt_out <-
                      dplyr::bind_rows(max_RIF, min_RIF) %>%
                     dplyr::mutate(top_level = top_name,
                                    pairwise_comparison = pair_name)
                    RIF_filt_out
                  })

    #################################################################### #
    # 4. Create summary data
    # Collapse together into 1 dataframe to run summary,
    # Need gene_ensembl, gene_name, Top_levels, top_levels_n, Pairwise_comparisons
    #reduce
    RIF_df <- filtered_RIF_data %>%
      purrr::reduce(dplyr::full_join,
                    by = c("gene_ensembl",
                           "gene_name",
                           "mean_RIF",
                           "top_level",
                           "pairwise_comparison"))


    #Summarise
    RIF_filtered_summarised <-
      RIF_df %>%
      dplyr::group_by(.data$gene_ensembl, .data$gene_name) %>%
      dplyr::summarize(top_groups_included = paste(sort(unique(.data$top_level)),collapse=", "),
                       number_top_groups = dplyr::n_distinct(.data$top_level),
                       comparisons_included = paste(sort(unique(paste(.data$top_level,
                                                                      .data$pairwise_comparison,
                                                                      sep = "__"))),
                                                    collapse = ", ")) %>%
      dplyr::arrange(dplyr::desc(.data$number_top_groups))

    RIF_filtered_summarised
  }





### make a full 'prepare' function, with export of data for use in cytoscape
#internal list that each 'selection' is appended to. E.g. user enters PIF, DE, RIF and/or TS (haven't done TS yet)

#' Prepare PCIT
#'
#' Function to prepare various data sources for selection for a Partical
#' Corellation Information Theory (PCIT) analysis. It is not necessarily useful
#' to use all DE, significant PIF or significant RIF values if there are many
#' hundreds or thousands as input to PCIT. This function will filter for
#' significance (PIF and RIF) and return a number of genes, based on user input.
#'
#'
#'
#'
#' @param data_type_selection A character vector of one or any of: "DE", "PIF",
#'   "RIF".
#' @param norm_exp_data Whole data normalisation output - use VST norm.
#' @param DE_full_out Full output of [auto_generate_DE_results].
#' @param RIF_output_data Full output of [calculate_RIF], if "RIF" included in
#'   \code{data_type_selection}.
#' @param gene_annotations Output from call to [annotate_gene_ensembl]
#' @param search_top_level character vector of top_level_groups (see
#'   [auto_generate_DE_results]) to search for, returning only tables from these
#'   groups. e.g. \code{c("ARC","LIV")}. If not provided it will return all top
#'   level groups in output.
#' @param filter_norm_data_columns logical. If \code{TRUE}, the
#'   \code{norm_exp_data} will be filtered to only include samples that are
#'   included in the top_level_groups selected by \code{search_top_level}. If
#'   \code{TRUE} and {search_top_level} is not provided or NULL, it will return
#'   only columns for top_level_groups that are within the DE_full_out.
#'   Sometimes a user may want to select genes on criteria from only a certain
#'   group, but still calculate PCIT across all normalised data.
#' @param DEn Number of values to return, ranked on adjusted p value. E.g. if
#'   \code{DEn = 20}, then the 20 genes with the lowest p value are returned
#'   (there may be less than this due to prior filtering for significance).
#' @param PIF_sig Number. Threshold used for filtering the z-scored PIF values,
#'   so a \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas
#'   a value of 1.96 filters for a nominal p value of 0.05.
#' @param PIFn Number of values to return for both the max PIF and min PIF. E.g.
#'   if \code{PIFn = 20}, then the 20 highest and 20 lowest are returned (if
#'   that many are present after filtering by PIF_sig)
#' @param RIF_sig Number. Threshold used for filtering the z-scored PIF values,
#'   so a \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas
#'   a value of 1.96 filters for a nominal p value of 0.05.
#' @param RIFn Number of values to return for both the max PIF and min PIF. E.g.
#'   if \code{PIFn = 20}, then the 20 highest and 20 lowest are returned (if
#'   that many are present after filtering by PIF_sig)
#' @param export_tables logical. \code{TRUE} indicates the tables of selected
#'   genes and filtered normalised counts tables are exported to
#'   \code{export_dir}. This is exported in formats useful for input to
#'   Cytoscape software, or for running PCIT using FORTRAN code (or windows
#'   executable) instead of in R.
#' @param export_venn_plot logical. \code{TRUE} indicates a venn diagram
#'   summarising the catagories of selected genes is exported to
#'   \code{export_dir}. This can only be exported, not viewed interactively
#'   within RStudio.
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created. Defaults to "./outputs/".
#'
#' @return Returns a list of dataframes of filtered output for all selected data
#'   types, as well as a normalised counts table filtered to include only
#'   selected genes, for use in PCIT.
#'
#' @export

PCIT_prepare_data <-
  function(data_type_selection = c("DE", "PIF", "RIF"),
           norm_exp_data,
           DE_full_out,
           RIF_output_data,
           gene_annotations,
           search_top_level = NULL,
           filter_norm_data_columns = TRUE,
           DEn,
           PIF_sig,
           PIFn,
           RIF_sig,
           RIFn,
           export_tables = FALSE,
           export_venn_plot = FALSE,
           export_dir = "./outputs/"
  ){
    if("DE" %in% data_type_selection){
      DE <-
        .prep_PCIT_DE(DE_full_out,
                      search_top_level,
                      DEn)
    } else {
      DE <- NULL
    }

    if("PIF" %in% data_type_selection){
      PIF <-
        .prep_PCIT_PIF(DE_full_out,
                     search_top_level,
                     PIF_sig,
                     PIFn)
    } else {
      PIF <- NULL
    }

    if ("RIF" %in% data_type_selection){
      RIF <-
        .prep_PCIT_RIF(RIF_output_data,
                     search_top_level,
                     RIF_sig,
                     RIFn)
    } else {
     RIF <- NULL
    }

    out_list <-
      list(RIF = RIF, PIF = PIF, DE = DE) %>%
      purrr::discard(is.null)

    selected_genes_detailed <-
      out_list %>%
      dplyr::bind_rows(.id = "selection") %>%
      dplyr::group_by(.data$gene_ensembl, .data$gene_name) %>%
      dplyr::summarise(selection = paste(sort(unique(.data$selection)), collapse = ", "),
                        top_groups_included = paste(sort(unique(.data$top_groups_included)), collapse = ", "),
                        comparisons_included = paste(sort(unique(.data$comparisons_included)), collapse = ", "),
                        selection_n = dplyr::n()
      ) %>%
      dplyr::left_join(gene_annotations, by = c("gene_ensembl", "gene_name")) %>%
      dplyr::arrange(dplyr::desc(.data$selection_n))

    ################################################################### #
    # Filter Normalised Counts
    #
    # If search level filter = NULL, then match normalised data to only include
    # the possible top levels from DE_out object.

    if(filter_norm_data_columns == TRUE){
      if(is.null(search_top_level)){
        top_level_names_in_DE <-
          names(DE_full_out) %>%
          stringr::str_extract(pattern = "^[^_]*")

        norm_selected_cols <-
          norm_exp_data %>%
          dplyr::select(where(is.character), tidyselect::contains(top_level_names_in_DE))

      } else{
        norm_selected_cols <-
          norm_exp_data %>%
          dplyr::select(where(is.character), tidyselect::contains(search_top_level))
      }
    } else{
      norm_selected_cols <- norm_exp_data
    }

    #Filter based on gene ensembls
    norm_selected <-
      norm_selected_cols %>%
      dplyr::filter(.data$gene_ensembl %in% selected_genes_detailed$gene_ensembl)


    message(crayon::green(paste("Number of selected genes within dataset: ", length(norm_selected$gene_ensembl))))
    message(crayon::green(paste("Number of samples within dataset: ", norm_selected %>% dplyr::select(where(is.numeric)) %>% ncol)))


    ################################################################### #
    # Check export dir, if needed.
    if(export_venn_plot == TRUE || export_tables == TRUE){
      if(!dir.exists(export_dir)){
        dir.create(export_dir, recursive = TRUE)
        message(crayon::red(paste("Directory created:", export_dir)))
      } else{message(crayon::green(paste(export_dir,"Directory exists")))}

    }


    #################################################################### #
    # Venn diagram

    if(export_venn_plot == TRUE){

      #check if more than 2 conditions
      if(length(data_type_selection) < 2){
        message(crayon::red("Cannot make Venn Diagram when length(data_type_selection) is less than 2."))
      } else {

        #check if required package is installed
        if (!requireNamespace("VennDiagram", quietly = TRUE)) {
          stop("Package \"VennDiagram\" needed for plotting Venn diagram. Please install it using install.packages(\"VennDiagram\").",
               call. = FALSE)
        }
        if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
          stop("Package \"RColorBrewer\" needed for plotting Venn diagram. Please install it using install.packages(\"RColorBrewer\").",
               call. = FALSE)
        }


        # Check if export dir exists

        all_selected_list <-
          out_list %>%
          purrr::map("gene_ensembl")


        myCol <- RColorBrewer::brewer.pal(length(all_selected_list), "Set1")

        VennDiagram::venn.diagram(
          x = all_selected_list,
          width = 3100,
          margin = .1,
          category.names = names(all_selected_list),
          filename = paste0(export_dir, "_Venn_genes_for_PCIT_",format(Sys.time(), "%Y%m%d_%H%M"),".png"),
          output = TRUE,
          fill = myCol,
          disable.logging = FALSE
        )
      }
    }
    ########################################################################
    # Export tables, in format for use in cytoscape/ fortran input
    #
    #export annotation file for cytoscape

    if(export_tables == TRUE){

    #export selected genes table
      data.table::fwrite(selected_genes_detailed,
           file = paste0(export_dir,
                         "PCIT_selected_genes_summary_",
                         format(Sys.time(), "%Y%m%d_%H%M"),".txt"),
           sep = "\t")

    #export normalised counts, filtered for PCIT and without annotation
      # 1. Using gene_ensembl
      # 2. Using gene_name


      #export for PCIT (without annotation) - ENSEMBL

      data.table::fwrite(norm_selected %>% dplyr::select(.data$gene_ensembl, where(is.double)),
             file = paste0(export_dir,
                           "PCIT-selected_norm_expression_GENE_ENSEMBL_",
                           format(Sys.time(), "%Y%m%d_%H%M"),".txt"),
             sep = " ",
             col.names = F)



      #export for PCIT (without annotation) - GENE_NAME


      data.table::fwrite(norm_selected %>% dplyr::select(.data$gene_name, where(is.double)),
             file = paste0(export_dir,
                           "PCIT-selected_norm_expression_GENE_NAME_",
                           format(Sys.time(), "%Y%m%d_%H%M"),".txt"),
             sep = " ",
             col.names = F)
      }

    out_list <- list("PCIT_summary_table" = selected_genes_detailed,
                     "PCIT_normalised_counts" = norm_selected)
   return(out_list)
  }






#' Run PCIT calculation
#'
#' Wrapper for [CeTF::PCIT] which is based on Reverter & Chan 2008. This can take
#' a long time, especially for larger datasets.
#'
#' @param norm_counts_for_PCIT Normalised counts for selected genes. Normally
#' this is produced by a call to [PCIT_prepare_data].
#' @param key_colname string of the name of column to be used as key, typically
#' \code{"gene_ensembl"} or \code{"gene_name"}.
#'
#'
#' @return Returns a list of:
#'
PCIT_calculate <-
  function(norm_counts_for_PCIT,
           key_colname = "gene_name")
  {
    if (!requireNamespace("CeTF", quietly = TRUE)) {
      stop("Package \"CeTF\" needed for RIF function to work. Please install it using BiocManager::install(\"CeTF\").",
           call. = FALSE)
    }


    if(inherits(norm_counts_for_PCIT, "data.frame") == FALSE){
      stop("norm_counts_for_PCIT is not of class 'data.frame'")
    }


    norm_selected <-
      norm_counts_for_PCIT %>%
      dplyr::select(.data[[key_colname]],
                    where(is.double))

    duplicated_rownames <- which(duplicated(norm_selected[[key_colname]]))
    if(length(duplicated_rownames) > 0){
      message(crayon::yellow(paste0("Removing rows with duplicated rownames: ",
                                    paste(norm_selected[[key_colname]][duplicated_rownames], collapse = ", "))))
      norm_selected <-
        norm_selected %>%
        dplyr::distinct(.data[[key_colname]],
                        .keep_all = TRUE)
    }

    pcit_in_matrix <-
      norm_selected %>%
        tibble::column_to_rownames(key_colname)

    pcit_out <-
    CeTF::PCIT(pcit_in_matrix,
               tolType = "mean")

    return(pcit_out)
  }








#' Filter PCIT output for Cytoscape
#'
#' Function that takes the output from PCIT analysis and filters out insignificant
#' correlations, and further filters based on a hard threshold to allow plotting
#' only the most meaningful correlations (both positive and negative) in Cytoscape.
#'
#' @param PCIT_out_data Full output of [PCIT_calculate]. A list of 3 elements: tab, adj_raw, adj_sig.
#' @param cor_threshold Threshold for filtering PCIT correlations (both positive and negative)
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created. Defaults to "./outputs/".
#'
#' @return Returns a ggplot/ggarrange object of density plots based on input
#' PCIT data. Also exports the required data for use in Cytoscape to
#'  \code{export_dir}.
#'

PCIT_filter_output <-
  function(PCIT_out_data,
           cor_threshold,
           export_dir
  ){

    #check CeTF is installed
    if (!requireNamespace("CeTF", quietly = TRUE)) {
      stop("Package \"CeTF\" needed for RIF function to work. Please install it using BiocManager::install(\"CeTF\").",
           call. = FALSE)
    }

    #select only significant differences
    full_output <- PCIT_out_data$tab %>% dplyr::filter(.data$corr2 != 0)

    #filter out only > or < 0.9
    filtered_output <- full_output %>%
      dplyr::filter(.data$corr1 >= cor_threshold | .data$corr1 <= -cor_threshold) %>%
      dplyr::select(-.data$corr2)

    #check number of nodes
    message(crayon::red("number of nodes:"))
    message(crayon::red(length(unique(c(filtered_output$gene1,filtered_output$gene2)))))


    filtered_output <- filtered_output %>% dplyr::rename(PCIT_cor = .data$corr1)

    ################################################################### #
    # Check export dir
    if(!dir.exists(export_dir)){
        dir.create(export_dir, recursive = TRUE)
        message(crayon::red(paste("Directory created:", export_dir)))
      } else{message(crayon::green(paste(export_dir,"Directory exists")))}

    ################################################################### #
    # output
    data.table::fwrite(filtered_output,
           paste0(export_dir,
                  "PCIT_FILTERED_",
                  cor_threshold, "_",
                  format(Sys.time(), "%Y%m%d_%H%M"),
                  ".txt"),
           sep = "\t")

    message(crayon::green(paste0("Exported data to : ",
                                 export_dir,
                                 "PCIT_FILTERED_",
                                 cor_threshold, "_",
                                 format(Sys.time(), "%Y%m%d_%H%M"),
                                 ".txt")))

    message(crayon::green("Creating summary density plots..."))

    ################################################################### #
    # Plot and return results
    CeTF::densityPlot(PCIT_out_data$adj_raw,
                      PCIT_out_data$adj_sig,
                      threshold = cor_threshold)
  }



### make a PCIT function that exports PCIT table for cytoscape

### make a function that combines the 2 with sensible defaults.
