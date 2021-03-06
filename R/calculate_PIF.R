#' Calculate PIF
#'
#' Calculates Phenotypic Impact Factor (PIF) for each pairwise comparison.
#'
#' The PIF value is the mean expression of Condition 1 and Condition 2 of the
#' pairwise comparison multiplied by the Log2 Fold Change. This favours genes
#' that are both highly abundant and highly different between the two
#' conditions. Designed to be used within the call to
#' \code{auto_generate_DE_results}, but can be used independently.
#'
#' @param DE_res_list List of DE results objects from DESEQ2. List names must be
#'   in the format Condition 1 vs Condition 2 e.g. "LIV_HCP.HP.UMEI vs
#'   LIV_LCP.LP.UMEI".
#' @param top_level_group_individual string to filter top level (e.g. Tissue Region).
#'   Inherits from \code{auto_generate_DE_results()}
#' @param log2_norm_data log2 Normalised data from current call to DESEQ2 (NOT normalised data from WHOLE dataset)
#' @param colData coldata normally accessed via a call to \code{SummarizedExperiment::colData()}
#' @param results_contrast_factor non-string. See \code{?auto_generate_DE_results}
#' @param export_tables logical. Should an excel file be exported of results? Inherits from \code{auto_generate_DE_results()}
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created. Inherits from \code{auto_generate_DE_results()}.
#'
#'
#' @return returns a list of dataframes with calculated mean expression, log2 FC and PIF columns.
#'
#' @export

calculate_PIF <- function(DE_res_list,
                                   top_level_group_individual,
                                   log2_norm_data,
                                   colData,
                                   results_contrast_factor,
                                   export_tables = FALSE,
                                   export_dir = "./outputs/"
                                   ){

  message(crayon::blue("PIF calculations - Preparing normalised data..."))

  cf <- rlang::enquo(results_contrast_factor)

  ########################################################################## #
  #Produce Mean Expression for each Treatment (results_contrast_factor)

  #Prepare normalised data for this run (top_level_group_individual)
  ## Select only the columns with expression data, and pivot longer
  col_names_current <- colnames(dplyr::select(log2_norm_data,
                                              tidyselect::starts_with(top_level_group_individual)))

  norm2 <- tidyr::pivot_longer(log2_norm_data,
                               cols = tidyselect::all_of(col_names_current),
                               names_to = "sample_names",
                               values_to = "expression")


  #Select results_contrast_factor and sample_names
  col_annot <- colData %>%
    as.data.frame() %>%
    dplyr::select(.data$sample_names, {{ cf }})

  # Add treatment annotation
  norm_annot <- dplyr::left_join(x = norm2,
                          y = col_annot,
                          by = "sample_names")


  #group by gene then results_contrast_factor and summarise with mean normalised counts
  norm_grouped <- norm_annot %>%
    dplyr::group_by(.data$gene_name, .data$gene_ensembl, {{ cf }})

  norm_summary <- dplyr::summarise(norm_grouped, mean = mean(expression, na.rm = T))

  #pivot table wider, so that for each row (gene) there is a a mean treatment expression
  norm_spread <- tidyr::pivot_wider(norm_summary,
                                    names_from = {{ cf }},
                                    values_from = mean)

  message(crayon::blue("PIF calculations - Preparing normalised data...COMPLETE"))
    ########################################################################## #
  # Iterate over each results object
  # 1. Get fold change for each pairwise comparison
  # 2. Get mean expression between the mean of Condition 1 and Condition 2
  # 3. Mean * FC = PIF

  #Define names of results list, as these names contain:
  #"Condition 1 vs Condition 2"
  DE_res_names <- names(DE_res_list)


  #function to iterate over, used with purrr::map2() so that object and name is
  # provided.

  .calculate_PIF_per_res_obj <- function(DE_res_object, DE_res_names){

    ################################ #
    #Prepare DE data
    #make names for filtering DE data
    DE_name_short <- stringr::str_remove_all(DE_res_names, pattern = " ")
    new_FC_col_name <- paste0("log2FC_", DE_name_short)
    new_P_col_name <- paste0("padj_", DE_name_short)

    #DE_res_object
    foldchange <- DE_res_object %>%
      BiocGenerics::as.data.frame() %>%
      tibble::rownames_to_column("gene_ensembl") %>%
      dplyr::select(.data$gene_ensembl,
                    !!new_FC_col_name := .data$log2FoldChange,
                    !!new_P_col_name := .data$padj)

    ################################ #
    # Join to normalised data
    norm_spread0 <- dplyr::left_join(norm_spread,
                              foldchange,
                              by = "gene_ensembl")

    ################################ #
    #Define new column names
    col1 <- stringr::str_split(DE_res_names, pattern = " vs ", simplify = T)[, 1] %>% rlang::sym()
    col2 <- stringr::str_split(DE_res_names, pattern = " vs ", simplify = T)[, 2] %>% rlang::sym()

    log2column <- paste0("log2FC_", col1, "vs", col2) %>% rlang::sym()
    padjcolumn <- paste0("padj_", col1, "vs", col2) %>% rlang::sym()

    mean_colname <- paste0("mean_", col1, "vs", col2) %>% rlang::sym()
    pif_colname <- paste0("PIF_", col1, "vs", col2) %>% rlang::sym()
    zpif_colname <- paste0("zPIF_", col1, "vs", col2) %>% rlang::sym()

    #print to console to show progress
    print(pif_colname)

    #Calculate PIF
    PIF_out <-
      norm_spread0 %>%
      dplyr::mutate(!!mean_colname := ((!!col1) + (!!col2)) / 2,
                    !!pif_colname := (!!mean_colname) * (!!log2column)
                    )  %>%
      dplyr::ungroup() %>%
      #Scale (z-score) the PIF, for downstream filtering for nominal significance
      dplyr::mutate(dplyr::across(.cols = dplyr::contains("PIF"),
                                  .fns = scale,
                                  .names = "z{col}")) %>%
      dplyr::select(.data$gene_name,
                    .data$gene_ensembl,
                    !!mean_colname,
                    !!log2column,
                    !!padjcolumn,
                    !!pif_colname,
                    !!zpif_colname,
                    !!.data$description) %>%
      dplyr::arrange(!!zpif_colname)
    }



  message(crayon::blue("PIF calculations - Calculating PIF data..."))

  PIF_out1 <- purrr::map2(.x = DE_res_list,
              .f = .calculate_PIF_per_res_obj,
              DE_res_names)

  PIF_out1  <- setNames(PIF_out1, paste0("PIF - ", DE_res_names))

  message(crayon::blue("PIF calculations - Calculating PIF data...COMPLETE"))


if(export_tables == TRUE){

  if(!dir.exists(export_dir)){
    dir.create(export_dir, recursive = TRUE)
    message(crayon::red(paste("Directory created:", export_dir)))
  } else{message(crayon::green(paste(export_dir,"Directory exists")))}

  PIF_out2 <- PIF_out1
  names(PIF_out2) <-
    names(PIF_out2) %>%
    stringr::str_remove_all("PIF - ") %>%
    stringr::str_squish() %>%
    stringr::str_remove_all(" ") %>%
    stringr::str_trunc(31,ellipsis = "")

  openxlsx::write.xlsx(PIF_out2,
                       file = paste0(export_dir,
                                     "PIF tables - SPLIT - ",
                                     top_level_group_individual,
                                     " - ",
                                     format(Sys.time(), "%Y%m%d_%H%M") ,
                                     ".xlsx"),
                       colWidths = "auto")

  message(crayon::green(paste("PIF tables exported to an .xlsx file in the sub-directory:", export_dir)))
}
  message(crayon::green("PIF calculations - END"))
  return(PIF_out1)
}
