#' Generate normalised counts and QC plots
#'
#' Internal function to use within \code{auto_generate_DE_results()}.
#'
#' This function generates normalised count tables, using the variance
#' stabilizing transformation technique (called "vsd" throughougt - variance
#' stabilised data), and also using log2 normilisation (log2norm). This function
#' also generates 3 types of quality control plots for the overall dataset.
#'
#' @param top_level_filter inherits from \code{auto_generate_DE_results()}
#' @param dds_object DESeq2 data object, inherits from
#'   \code{auto_generate_DE_results()}
#' @param gene_annotations Output from \code{annotate_gene_ensembl()}
#' @param export_tables logical. \code{TRUE} indicates the normalised counts
#'   tables, both vsd and log2norm, annotated with gene names and descriptions,
#'   are exported to \code{"./outputs/normalised_counts/"})
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created.
#' @param results_contrast_factor inherits from
#'   \code{auto_generate_DE_results()}
#'
#'
#' @return returns a list of 1) heatmap_plot, PCA_plot and Sample distances
#'   heatmap; and 2) normalised vsd and normalised log2norm dataframes.
#'
#' @export
#'

#'
#'
#'
#'
.norm_and_plot <- function(top_level_filter,
                           dds_object,
                           gene_annotations,
                           export_tables = FALSE,
                           export_dir = "./outputs/normalised_counts/",
                           results_contrast_factor){


  results_contrast_factor_string <- rlang::enquo(results_contrast_factor) %>% rlang::as_label()
  ################################################################################## #
  #Generate full, annotated table of results of DESeq2
  #
  #   DE_out$ARC_DESeq2_Output$DESeq2_res_object$`ARC_HCP.HP.RMEI vs ARC_HCP.HP.UMEI`
  #
  #
  #   #change all test to dds_object
  #
  #   test <- fake_output %>%
  #   unlist(recursive = FALSE) %>%
  #     {.[grepl("res_object", names(.))]}
  #
  #   top_level_names <- stringr::str_split(names(test), pattern = "_", n=2, simplify = TRUE)[,1]
  #
  #
  # test$LHA_DESeq2_Output.DESeq2_res_object$`ARC_HCP.HP.RMEI vs ARC_HCP.HP.UMEI`
  #
  #
  # #need a function within a function, as mutli layers of test$...
  #
  #
  # test_output <- purrr::map2(.x =test ,
  #             .y = top_level_names,
  #             .f = function(x, y){
  #
  #               as.data.frame(x) %>%
  #                 dplyr::mutate(top_level_filter = y,
  #                               )
  #             })
  #
  # test_output$LHA_DESeq2_Output.DESeq2_res_object
  #
  ################################################################################## #
  #VSD Normalisation

  #calculate normalised expression
  vsd <- DESeq2::vst(dds_object, blind=FALSE)
  vsd_a <-
    SummarizedExperiment::assay(vsd) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_ensembl") %>%
    dplyr::left_join(gene_annotations, by = c("gene_ensembl")) %>%
    dplyr::select(.data$gene_ensembl,
                  .data$gene_name,
                  .data$description,
                  tidyselect::everything())



  ################################################################################## #
  #log 2 Normalisation

  dds_temp <- BiocGenerics::estimateSizeFactors(dds_object) #this is the 'factor' that it normalises by below


  log2norm <- DESeq2::normTransform(object = dds_temp,
                                    f = log2,#log2 transformation
                                    pc = 1) #pseudocount

  # Annotate
  log2norm_a <-
    SummarizedExperiment::assay(log2norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_ensembl") %>%
    dplyr::left_join(gene_annotations, by = c("gene_ensembl")) %>%
    dplyr::select(.data$gene_ensembl,
                  .data$gene_name,
                  .data$description,
                  tidyselect::everything())



  ################################################################################## #
  #export tables of VST and log2 normalised data, plus the full list of DE output
  if(export_tables == TRUE){

    if(!dir.exists(export_dir)){
      dir.create(export_dir, recursive = TRUE)
      message(crayon::green(paste("Directory created:", export_dir)))
    } else{message(crayon::green(paste(export_dir,"Directory exists")))}

    data.table::fwrite(vsd_a ,
                       file = paste0(export_dir,top_level_filter,"_VST_normalisedcounts_",format(Sys.time(), "%Y%m%d_%H%M"),".txt"), row.names = F)

    data.table::fwrite(log2norm_a,
                       file = paste0(export_dir,top_level_filter,"_log2_of_DESEQ_internal_normalisedcounts_",format(Sys.time(), "%Y%m%d_%H%M"), ".txt"), row.names = F)

    message(crayon::green(paste("Normalised tables exported to the sub-directory:", export_dir)))
    # fwrite(full_dataset_DE,
    #        file = paste0("./outputs/",region,"_DE_output_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
  }


  ################################################################################## #
  #Create PCA plot
  rcf <- rlang::enquo(results_contrast_factor)
  pcaData <- DESeq2::plotPCA(vsd, intgroup=c(results_contrast_factor_string), returnData=TRUE, ntop = 500)

  percentVar <- round(100 * attr(pcaData, "percentVar"))

  PCA_plot <- ggplot2::ggplot(pcaData, ggplot2::aes(.data$PC1, .data$PC2, fill=!!rcf)) +
    ggplot2::geom_point(size=4.2, shape = 21, colour = "black", alpha = 0.8) +
    ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggplot2::coord_fixed()+
    #coord_fixed(xlim = c(-150,150), ylim = c(-120, 120)) +
    ggplot2::theme_classic()+
    ggrepel::geom_text_repel(mapping = ggplot2::aes(label=.data$name), max.overlaps = 10)+
    ggplot2::ggtitle(label = paste(top_level_filter))


  ################################################################################## #
  #Heatmap of count matrix
  topVarGenes <- BiocGenerics::order(genefilter::rowVars( SummarizedExperiment::assay(vsd) ), decreasing=TRUE )[1:30]

  heatmap_input <- SummarizedExperiment::assay(vsd)[ topVarGenes, ]

  heatmap_input <- heatmap_input - BiocGenerics::rowMeans(heatmap_input)

  df1 <- BiocGenerics::as.data.frame(SummarizedExperiment::colData(dds_object)[results_contrast_factor_string])

  heatmap_plot<- pheatmap::pheatmap(heatmap_input,
                                    cluster_rows = F,
                                    annotation_col = df1,
                                    main = paste(top_level_filter, "- Top most variable genes" ))

  ################################################################################## #
  #Heatmap of sample-to-sample distances
  #with dist(t(assay(vsd))) - see vignette

  sampleDists <- stats::dist(BiocGenerics::t(SummarizedExperiment::assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)

  rownames(sampleDistMatrix) <- SummarizedExperiment::colData(vsd)[,results_contrast_factor_string]
  colnames(sampleDistMatrix) <- NULL
  sample_sample_heatmap <- pheatmap::pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists)

  ################################################################################## #
  #OUTPUTS - to lists already generated




  #put plots into a list, append to the list (from global environment) but named
  #plots of overall data
  plots <- list( heatmap = heatmap_plot, PCA = print(PCA_plot), Sample_distances = sample_sample_heatmap)

  #put data into a list, append to the list (from global environment) but named
  norm_data_out <-  list(vsd = vsd_a, log2norm = log2norm_a)

  data_out <- stats::setNames(list(plots, norm_data_out),
                       c(paste0(top_level_filter, "_overall_plots"),
                         paste0(top_level_filter, "_normalised_data")))

  return(data_out)
}
