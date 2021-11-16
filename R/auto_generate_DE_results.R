#' Auto Generate DE Results
#'
#' Wrapper to DESeq2 that is designed to be used iteratively for each element of
#' \code {top_level_filter} e.g. for each Tissue Region separately.
#'
#' This wrapper uses a call to \code{.calculate_DESEQ_results} which determines
#' all possible pairwise comparisons, based on the \code{top_level_name}
#' attribute which is also a column of \code{colData}. Defaults to use a formula
#' design of \code{~Region_Diet}, but could be changed to be \code{~Region +
#' Diet} (see \code{vignette("DESeq2")})but the pairwise comparions are
#' generated with \code{results_contrast_factor}, which must be one of the
#' elements of the formula.
#'
#' @param top_level_filter string. Name of each element of \code{top_level_name}
#'   to be iterated over. Usually this is provided by a call to \code{map()} or
#'   \code{lapply()}.
#' @param se_data SummarizedExperiment object of experimental data.
#' @param top_level_name non-string. Name of \code{colData} column with top
#'   level factor, which is normally Tissue Region. Defaults to \code{Region}.
#' @param column_of_samples non-string. Name of \code{colData} column with
#'   sample names.
#' @param samples_to_remove string/s. Names of samples to remove e.g.
#'   \code{c("ARC_01_HCP-HP-RMEI", "ARC_02_HCP-HP-UMEI")}
#' @param rowSums_filter numeric. For intial filter of rowSums across whole dds.
#'   Filters out all samples that have basically 0 expression across all
#'   samples. Saves computing time, but does not change DESeq results as it
#'   filters internally as well. Defaults to 10.
#' @param DESeq2_formula_design = formula in form \code{~Factor} parsed to
#'   \code{DESeq()}. Defaults to \code{~Region_Diet}.
#' @param results_contrast_factor Name of Factor from
#'   \code{DESeq2_formula_design} that is used to generate pairwise comparisons.
#'   Defaults to \code{Region_Diet}
#' @param results_combinations Either \code{NA}, which will automatically
#'   generate all possible pairwise combinatiosn of the
#'   \code{results_contrast_factor} OR a list of pairwise comparisons with each
#'   element containing 2 strings, numerator and denominator, respectively. See
#'   \code{?DESeq2::results} for requirements. Defaults to \code{NA}.
#' @param use_IHW_filtering logical. Should Independent Hypothesis Weighting be
#'   used when calculating results? \code{TRUE}, the default, indicates yes. See
#'   \code{?IHW::ihw} or \code{vignette("DESeq2")} for more details.
#' @param alpha numeric. What is the p-value threshold to be used for
#'   determining significance. Used in call to \code{DESeq2::results()} and
#'   others.
#' @return returns a named, nested list (list of lists) with dds_wald_object,
#'   list of DESeq result objects and list of pairwise plots (contains MA plots
#'   and p value histograms).
#'
#' @export

auto_generate_DE_results <-

function(top_level_filter, #from lapply e.g. c("ARC", "LAT"), 1 at a time
         se_data,
         top_level_name = Region,
         column_of_samples,
         samples_to_remove = NA,
         DESeq2_formula_design = ~Region_Diet,
         rowSums_filter = 10, #for dds filtering
         results_contrast_factor = Region_Diet,
         results_combinations = NA,
         use_IHW_filtering = TRUE,
         alpha = 0.05
){


  #check se_data is of class SummarizedExperiment
  if(inherits(se_data, "SummarizedExperiment") == FALSE){
    stop("se_data is not of class 'SummarizedExperiment'")
  }

  #### 1. subset data


  #rename column names as strings for subsetting with []
  top_level_name_string <- rlang::enquo(top_level_name) %>% rlang::as_label()
  column_of_samples_string <- rlang::enquo(column_of_samples) %>% rlang::as_label()


  #subset out top_level_name (e.g. Region)
  se_data0 <-
    S4Vectors::subset(se_data,
                      select = SummarizedExperiment::colData(se_data)[[top_level_name_string]] == top_level_filter)



  #subset out samples to remove, if applicable
  if(any(is.na(samples_to_remove)) == FALSE){

    if(any(SummarizedExperiment::colData(se_data0)[[column_of_samples_string]] %in% samples_to_remove)){
      se_data0 <-
        SummarizedExperiment::subset(se_data0,
                                     select = !(SummarizedExperiment::colData(se_data0)[[column_of_samples_string]] %in% samples_to_remove))


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

  boxplot_cooks_distance  <<- log10(SummarizedExperiment::assays(dds_wald)[["cooks"]]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl") %>%
    tidyr::pivot_longer(where(is.numeric), values_to = "cooks_distance") %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = cooks_distance, y = name))+
    ggplot2::geom_boxplot()+
    ggplot2::theme_classic()
  message(crayon::blue("Cooks distance plot complete."))



  #### 5. Produce results
  ################################################################################ #

  message(crayon::blue("Generating pairwise DESeq2 results..."))

  cf <- rlang::enquo(results_contrast_factor)

  a <- alpha

  res_out <-   .calculate_DESEQ_results(dds_object = dds_wald,
                                        contrast_factor = !!cf,
                                        alpha = a,
                                        use_IHW_filtering = use_IHW_filtering)


  message(crayon::blue("Results generated."))
  ################################################################################ #
  #mutate descriptive columns to dataframes and create one long dataset of all results
  list_results <- res_out
  ################################################################################## #

  # make MA plots (rough DESeq ones)

  #iterates over all results objects to make plots
  message(crayon::blue("Plotting MA plots (default DESeq2 style for QC)..."))
  MA_pairwise_plots <-
    purrr::map2(.x = list_results, .y = names(list_results),
                function(x,names){
                  plot_title = paste(top_level_filter, names)
                  DESeq2::plotMA(x, main = plot_title)
                  MA_plot <- recordPlot()
                })
  message(crayon::blue("Finished MA plots."))
  ################################################################################## #
  #pvalue histograms - (updated to be x$baseMean >1)

  #iterates over all results objects to make pvalue histograms (consider update with )
  message(crayon::blue("Plotting p value histograms (for QC)..."))
  pvalue_histogram_pairwise_plots <<-
    purrr::map2(.x = list_results, .y = names(list_results),
                function(x,names){
                  hist(x$pvalue[x$baseMean > 1],
                       breaks = 0:20/20,
                       col = "grey50",
                       border = "white",
                       main = paste("p values (with baseMean > 1) of ", names))
                  pvalue_histogram <- recordPlot()
                })
  message(crayon::blue("Finished plotting p value histograms."))
  ################################################################################## #
  #OUTPUT
  # Prepare data for export
  #add's plots into 1 named list for export
  message(crayon::blue("Preparing data for output..."))
  pairwise_plots0 <- list(MA_plots = MA_pairwise_plots, Pvalue_histogram = pvalue_histogram_pairwise_plots)

  list_out <- list(dds_wald_object = dds_wald, DESeq2_res_object = res_out, pairwise_plots = pairwise_plots0)

  #Add to a list with it's own name identifying it by top_level_filter (region)
  list_out <- setNames(list(list_out), paste0(top_level_filter, "_DESeq2_Output"))
  message(crayon::red("List output succesfully generated. Contains: dds_wald_object [dds after call to DESeq()], list of DESeq result objects and list of pairwise plots."))
  return(list_out)
}
