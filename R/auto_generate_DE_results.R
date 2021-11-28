#' Auto Generate DE Results
#'
#' Wrapper to DESeq2 that iterates over each element of \code{top_level_groups}
#' and conducts pairwise comparisons, QC plots and other analysis. Can also
#' normalise data over whole dataset (ignoring \code{top_level_groups})
#'
#' This wrapper determines all possible pairwise comparisons, based on the
#' \code{top_level_colname} attribute which is also a column of \code{colData}.
#' Defaults to use a formula design of \code{~Region_Diet}, but could be changed
#' to be \code{~Region + Diet} (see \code{vignette("DESeq2")}) but the pairwise
#' comparisons are generated with \code{results_contrast_factor}, which must be
#' one of the elements of the formula.
#'

#' @param se_data SummarizedExperiment object of experimental data.
#' @param top_level_groups character vector of names of each element of \code{top_level_colname}
#'   to be iterated over.
#' @param top_level_colname non-string. Name of \code{colData} column with top
#'   level factor, which is normally "Tissue_Region".
#' @param sample_colname non-string. Name of \code{colData} column with
#'   sample names (used for filtering).
#' @param samples_to_remove string/s. Names of samples to remove e.g.
#'   \code{c("ARC_01_HCP-HP-RMEI", "ARC_02_HCP-HP-UMEI")}
#' @param rowSums_filter numeric. For intial filter of rowSums across whole dds.
#'   Filters out all samples that have basically 0 expression across all
#'   samples. Saves computing time, but does not change DESeq results as it
#'   filters internally as well. Defaults to 10.
#' @param DESeq2_formula_design = formula in form \code{~Factor} parsed to
#'   \code{DESeq()}.
#' @param results_contrast_factor Name of Factor from
#'   \code{DESeq2_formula_design} that is used to generate pairwise comparisons.
#' @param results_combinations Either \code{NA}, which will automatically
#'   generate all possible pairwise combinatios of the
#'   \code{results_contrast_factor} OR a list of pairwise comparisons with each
#'   element containing 2 strings, numerator and denominator, respectively. See
#'   \code{?DESeq2::results} for requirements. Defaults to \code{NA}.
#' @param use_IHW_filtering logical. Should Independent Hypothesis Weighting be
#'   used when calculating results? \code{TRUE}, the default, indicates yes. See
#'   [IHW::ihw] or \code{vignette("DESeq2")} for more details.
#' @param alpha numeric. What is the p-value threshold to be used for
#'   determining significance. Used in call to [DESeq2::results] and
#'   others.
#' @param gene_annotations Output from call to [annotate_gene_ensembl]
#' @param export_tables logical. \code{TRUE} indicates the normalised counts
#'   tables, both vsd and log2norm, annotated with gene names and descriptions,
#'   are exported to \code{"./outputs/normalised_counts/"}
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created.
#' @param whole_data_normalisation logical. Should DESEQ2 be run for pairwise
#'  comparisons, split by \code{top_level_groups}? or, if FALSE - normalised
#'  output only for whole dataset.
#'
#' @return returns a named, nested list (list of lists) with: \itemize{
#' \item dds_wald_object \item boxplot_cooks_distance \item DESeq2_res_object
#' \item pairwise_plots \item overall_plots \item normalised_data \item PIF
#' \item DE_by_PIF_df}
#'
#' @export

auto_generate_DE_results <-

function(se_data,
         top_level_groups,
         top_level_colname,
         sample_colname,
         samples_to_remove = NA,
         DESeq2_formula_design,
         rowSums_filter = 10, #for dds filtering
         results_contrast_factor,
         results_combinations = NA,
         use_IHW_filtering = TRUE,
         alpha = 0.05,
         gene_annotations,
         export_tables = FALSE,
         export_dir = "./outputs/",
         whole_data_normalisation = FALSE
){


  #check se_data is of class SummarizedExperiment
  if(inherits(se_data, "SummarizedExperiment") == FALSE){
    stop("se_data is not of class 'SummarizedExperiment'")
  }

  #Helpers
  sample_colname_string <- rlang::enquo(sample_colname) %>% rlang::as_label()
  cf <- rlang::enquo(results_contrast_factor)
  tln <- rlang::enquo(top_level_colname)
  a <- alpha


  if(whole_data_normalisation == TRUE){
    message(crayon::green("Whole data normalisation selected. No pairwise results will be generated."))
    se_data0 <- se_data
    top_level_groups <- "Whole_Dataset"
  }


  .f_DE_per_group <-
    function(top_level_group_individual){
      #### 1. subset data
      if(whole_data_normalisation == FALSE){
        #Helpers
        top_level_colname_string <- rlang::enquo(top_level_colname) %>% rlang::as_label()

        message(crayon::black$bgCyan$bold(paste("\n\n ******************* Start of - ",
                                                top_level_group_individual,
                                                "******************* \n\n")))

        #subset out top_level_colname (e.g. Region)
        se_data0 <-
          S4Vectors::subset(se_data,
                            select = SummarizedExperiment::colData(se_data)[[top_level_colname_string]] == top_level_group_individual)

      }

      #subset out samples to remove, if applicable
      if(any(is.na(samples_to_remove)) == FALSE){

        if(any(SummarizedExperiment::colData(se_data0)[[sample_colname_string]] %in% samples_to_remove)){
          se_data0 <-
            SummarizedExperiment::subset(se_data0,
                                         select = !(SummarizedExperiment::colData(se_data0)[[sample_colname_string]] %in% samples_to_remove))


          message(crayon::green(paste("Filtered out samples: ", paste(samples_to_remove, collapse = ", "))))
          message(crayon::green("New dimensions: "))
          print(se_data0)
        }
      }

      ### 2. create dds object and filter counts
      dds <- DESeq2::DESeqDataSet(se_data0,
                                  design = DESeq2_formula_design)

      keep <- rowSums(BiocGenerics::counts(dds)) >= rowSums_filter
      paste("Number of genes with more than", rowSums_filter, "counts accross all samples:  ", length(which(keep == T)))

      dds <- dds[keep,]

      #add a note on dimensions here?


      #### 3.
      message(crayon::blue("Beginning DESeq analysis..."))
      #Run DESeq with Wald - pairwise
      dds_wald <- DESeq2::DESeq(dds, test = "Wald",
                                betaPrior = F,
                                minReplicatesForReplace = 6,
                                parallel = FALSE)#this is default to n = 7; however this replace outlier gene values with a conservative estimate - this is very important for downstream PIF analysis. If n = <5/treatment it would be best to remove gene entirely.

      message(crayon::blue("Completed DESeq analysis."))


      #### 4.
      message(crayon::blue("Plotting cooks distance..."))
      #Check cooks distance boxplot to see if any gene is higher on average
      #boxplot(log10(assays(dds_wald)[["cooks"]]), range=0, las=2)

      boxplot_cooks_distance  <- log10(SummarizedExperiment::assays(dds_wald)[["cooks"]]) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("ensembl") %>%
        tidyr::pivot_longer(where(is.numeric), values_to = "cooks_distance") %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = .data$cooks_distance, y = .data$name))+
        ggplot2::geom_boxplot()+
        ggplot2::theme_classic()
      message(crayon::blue("Cooks distance plot complete."))



      #### 5. Produce results
      ################################################################################ #
      if(whole_data_normalisation == FALSE){
        message(crayon::blue("Generating pairwise DESeq2 results..."))



        res_out <-   .calculate_DESEQ_results(dds_object = dds_wald,
                                              contrast_factor = !!cf,
                                              alpha = a,
                                              use_IHW_filtering = use_IHW_filtering)


        message(crayon::blue("Results generated."))



        ################################################################################## #
        # make MA plots (rough DESeq ones)

        #iterates over all results objects to make plots
        message(crayon::blue("Plotting MA plots (default DESeq2 style for QC)..."))
        MA_pairwise_plots <-
          purrr::map2(.x = res_out, .y = names(res_out),
                      function(x,names){
                        plot_title = paste(top_level_group_individual, names)
                        DESeq2::plotMA(x, main = plot_title)
                        MA_plot <- grDevices::recordPlot()
                      })
        message(crayon::blue("Finished MA plots."))
        ################################################################################## #
        #pvalue histograms - (updated to be x$baseMean >1)

        #iterates over all results objects to make pvalue histograms (consider update with )
        message(crayon::blue("Plotting p value histograms (for QC)..."))
        pvalue_histogram_pairwise_plots <-
          purrr::map2(.x = res_out, .y = names(res_out),
                      function(x,names){
                        graphics::hist(x$pvalue[x$baseMean > 1],
                                       breaks = 0:20/20,
                                       col = "grey50",
                                       border = "white",
                                       main = paste("p values (with baseMean > 1) of ", names))
                        pvalue_histogram <- grDevices::recordPlot()
                      })
        message(crayon::blue("Finished plotting p value histograms."))

      }
      ################################################################################## #
      # Generate normalised counts tables, and plots of overall dataset

      norm_out <- .norm_and_plot(top_level_group_individual,
                                 dds_object = dds_wald,
                                 gene_annotations,
                                 export_tables,
                                 export_dir,
                                 results_contrast_factor = !!cf,
                                 whole_data_normalisation,
                                 top_level_colname = !!tln)
      #accesses two separate parts of output, which already contain the region name in them so are unique
      overall_plots_list <- norm_out[grepl("overall_plots", names(norm_out))][[1]]
      normalised_data_list <- norm_out[grepl("normalised_data", names(norm_out))][[1]]





      if(whole_data_normalisation == FALSE){
        ################################################################################ #
        # Calculate PIF
        PIF_out <-
          calculate_PIF(DE_res_list = res_out, #named list
                        top_level_group_individual,
                        log2_norm_data = normalised_data_list$log2norm, #from within same run
                        colData = SummarizedExperiment::colData(se_data0),
                        results_contrast_factor = !!cf,
                        export_tables,
                        export_dir)

        ################################################################################ #
        # mutate descriptive columns to dataframes
        # filter on adjusted p value to significance threshold
        # adding PIF value to arrange on PIF
        list_results <- res_out
        list_PIF_results <- PIF_out #should be in same order, as res_out was iterated over to generate PIF_out

        DE_annotated <-
          purrr::map2(.x = list_results,
                      .y = list_PIF_results,
                      function(.x, .y, annot, p_thresh){

                        PIF <- .y %>%
                          dplyr::select(.data$gene_ensembl, tidyselect::starts_with("zPIF")) %>%
                          dplyr::rename_with(function(x){
                            stringr::str_trunc(x,4,"right", "")},
                            tidyselect::starts_with("zPIF")
                          )

                        .x %>%
                          BiocGenerics::as.data.frame() %>%
                          tibble::rownames_to_column(var = "gene_ensembl") %>%
                          dplyr::left_join(annot, by = c("gene_ensembl")) %>%
                          dplyr::left_join(PIF, by = c("gene_ensembl")) %>%
                          dplyr::select(.data$gene_ensembl,
                                        .data$gene_name,
                                        .data$description,
                                        tidyselect::everything()) %>%
                          dplyr::filter(.data$padj <= p_thresh) %>%
                          dplyr::arrange(dplyr::desc(abs(.data$zPIF)))#descending, absolute

                      },
                      annot = gene_annotations,
                      p_thresh = alpha)

        ################################################################################ #
        # export DE Tables
        #
        if(export_tables == TRUE){

          if(!dir.exists(export_dir)){
            dir.create(export_dir, recursive = TRUE)
            message(crayon::red(paste("Directory created:", export_dir)))
          } else{message(crayon::green(paste(export_dir,"Directory exists")))}

          DE_annotated2 <- DE_annotated
          names(DE_annotated2) <-
            names(DE_annotated2) %>%
            stringr::str_squish() %>%
            stringr::str_remove_all(" ") %>%
            stringr::str_trunc(31,ellipsis = "")


          openxlsx::write.xlsx(DE_annotated2,
                               file = paste0(export_dir,
                                             "DE genes sorted by PIF - SPLIT by pairwise - ",
                                             top_level_group_individual,
                                             " - ",
                                             format(Sys.time(), "%Y%m%d_%H%M") ,
                                             ".xlsx"),
                               colWidths = "auto")

          message(crayon::green(paste("DE tables exported to an .xlsx file in the sub-directory:", export_dir)))

        }
      }


      ################################################################################## #
      #OUTPUT - results
      # Prepare data for export
      #add's plots into 1 named list for export
      message(crayon::blue("Preparing data for output..."))
      if(whole_data_normalisation == FALSE){
        pairwise_plots0 <- list(MA_plots = MA_pairwise_plots, Pvalue_histogram = pvalue_histogram_pairwise_plots)

        list_out <- list(dds_wald_object = dds_wald,
                         boxplot_cooks_distance = boxplot_cooks_distance,
                         DESeq2_res_object = res_out,
                         pairwise_plots = pairwise_plots0,
                         overall_plots = overall_plots_list,
                         normalised_data = normalised_data_list,
                         PIF = PIF_out,
                         DE_by_PIF_df = DE_annotated)


        #Add to a list with it's own name identifying it by top_level_group_individual (region)
        list_out <- stats::setNames(list(list_out), paste0(top_level_group_individual, "_DESeq2_Output"))
        message(crayon::red("List output succesfully generated."))

        message(crayon::black$bgCyan$bold(paste("\n\n ******************* END of - ",
                                                top_level_group_individual,
                                                "******************* \n\n")))
      } else if(whole_data_normalisation == TRUE){
        list_out <- list(dds_wald_object = dds_wald,
                         boxplot_cooks_distance = boxplot_cooks_distance,
                         overall_plots = overall_plots_list,
                         normalised_data = normalised_data_list)


        #Add to a list with it's own name identifying it by top_level_group_individual (region)
        list_out <- stats::setNames(list(list_out), paste0("Whole_data_normalisation_output"))
        message(crayon::red("List output succesfully generated."))

        message(crayon::black$bgCyan$bold(paste("\n\n ******************* END ******************* \n\n")))

      }
      #print("checking names output:")
      #print(names(list_out))
      return(list_out)
    }

  #Character vector to unnamed list
  top_level_groups <- top_level_groups %>% unique() %>% as.list()

  # map
  out_list <- purrr::map(top_level_groups,
                         .f_DE_per_group)

  #remove_top_level of list
  out_list <- out_list %>% unlist(recursive = FALSE)

  return(out_list)

}
