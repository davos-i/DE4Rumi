#' Generate normalised counts and QC plots
#'
#' Internal function to use within \code{auto_generate_DE_results()}.
#'
#' This function generates normalised count tables, using the variance
#' stabilizing transformation technique (called "vsd" throughougt - variance
#' stabilised data), and also using log2 normilisation (log2norm). This function
#' also generates 3 types of quality control plots for the overall dataset.
#'
#' @param top_level_group_individual inherits from \code{auto_generate_DE_results()}
#' @param dds_object DESeq2 data object, inherits from
#'   \code{auto_generate_DE_results()}
#' @param gene_annotations Output from \code{annotate_gene_ensembl()}
#' @param export_tables logical. \code{TRUE} indicates the normalised counts
#'   tables, both vsd and log2norm, annotated with gene names and descriptions,
#'   are exported to \code{"./outputs/normalised_counts/"})
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this specificed directory doesn't exist, it will be created.
#' @param results_contrast_factor inherits from
#'   \code{auto_generate_DE_results()}
#' @param top_level_colname inherits from \code{auto_generate_DE_results()}.
#'   non-string. Name of \code{colData} column with top level factor, which is
#'   normally Tissue Region.
#' @param whole_data_normalisation inherits from
#'   \code{auto_generate_DE_results()}.
#'
#'
#' @return returns a list of 1) heatmap_plot, PCA_plot and Sample distances
#'   heatmap; and 2) normalised vsd and normalised log2norm dataframes.
#'
#' @export
#'
#'


.norm_and_plot <- function(top_level_group_individual,
                           dds_object,
                           gene_annotations,
                           export_tables = FALSE,
                           export_dir = "./outputs/",
                           results_contrast_factor,
                           top_level_colname = NA,
                           whole_data_normalisation = FALSE){


  #Determines filters for plotting, which require different input if for whole data normalisation from auto_generate_DE_results
if(whole_data_normalisation == FALSE){
  results_contrast_factor_string <- rlang::enquo(results_contrast_factor) %>% rlang::as_label()
  rcf <- rlang::enquo(results_contrast_factor)
} else if(whole_data_normalisation == TRUE){
  rcf <- rlang::enquo(top_level_colname)
  results_contrast_factor_string <- rlang::enquo(top_level_colname) %>% rlang::as_label()
  }

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
      message(crayon::red(paste("Directory created:", export_dir)))
    } else{message(crayon::green(paste(export_dir,"Directory exists")))}

    data.table::fwrite(vsd_a ,
                       file = paste0(export_dir,top_level_group_individual,"_VST_normalisedcounts_",format(Sys.time(), "%Y%m%d_%H%M"),".txt"), row.names = F)

    data.table::fwrite(log2norm_a,
                       file = paste0(export_dir,top_level_group_individual,"_log2_of_DESEQ_internal_normalisedcounts_",format(Sys.time(), "%Y%m%d_%H%M"), ".txt"), row.names = F)

    message(crayon::green(paste("Normalised tables exported to the sub-directory:", export_dir)))
    # fwrite(full_dataset_DE,
    #        file = paste0("./outputs/",region,"_DE_output_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
  }


  ################################################################################## #
  #Create PCA plot

  pcaData <- DESeq2::plotPCA(vsd, intgroup=results_contrast_factor_string, returnData=TRUE, ntop = 500)

  percentVar <- round(100 * attr(pcaData, "percentVar"))

  PCA_plot <- ggplot2::ggplot(pcaData, ggplot2::aes(.data$PC1, .data$PC2, fill=!!rcf)) +
    ggplot2::geom_point(size=4.2, shape = 21, colour = "black", alpha = 0.8) +
    ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggplot2::coord_fixed()+
    #coord_fixed(xlim = c(-150,150), ylim = c(-120, 120)) +
    ggplot2::theme_classic()+
    #ggrepel::geom_text_repel(mapping = ggplot2::aes(label=.data$name), max.overlaps = 10)+
    ggplot2::ggtitle(label = paste(top_level_group_individual))


  ################################################################################## #
  #Heatmap of count matrix
  topVarGenes <- BiocGenerics::order(genefilter::rowVars( SummarizedExperiment::assay(vsd) ), decreasing=TRUE )[1:30]

  heatmap_input <- SummarizedExperiment::assay(vsd)[ topVarGenes, ]

  heatmap_input <- heatmap_input - BiocGenerics::rowMeans(heatmap_input)

  df1 <- BiocGenerics::as.data.frame(SummarizedExperiment::colData(dds_object)[results_contrast_factor_string])

  heatmap_plot<- pheatmap::pheatmap(heatmap_input,
                                    cluster_rows = F,
                                    annotation_col = df1,
                                    main = paste(top_level_group_individual, "- Top most variable genes" ))

  ################################################################################## #
  #Heatmap of sample-to-sample distances
  #with dist(t(assay(vsd))) - see vignette

  sampleDists <- stats::dist(BiocGenerics::t(SummarizedExperiment::assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)

  #rownames(sampleDistMatrix) <- SummarizedExperiment::colData(vsd)[,results_contrast_factor_string]
  rownames(sampleDistMatrix) <- SummarizedExperiment::colData(vsd)[,"sample_names"]
  colnames(sampleDistMatrix) <- NULL
  sample_sample_heatmap <- pheatmap::pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           main = paste(top_level_group_individual, "- Sample-Sample Distances"))

  ################################################################################## #
  #OUTPUTS - to lists already generated




  #put plots into a list, append to the list (from global environment) but named
  #plots of overall data
  plots <- list( heatmap = heatmap_plot, PCA = print(PCA_plot), Sample_distances = sample_sample_heatmap)

  #put data into a list, append to the list (from global environment) but named
  norm_data_out <-  list(vsd = vsd_a, log2norm = log2norm_a)

  data_out <- stats::setNames(list(plots, norm_data_out),
                       c(paste0(top_level_group_individual, "_overall_plots"),
                         paste0(top_level_group_individual, "_normalised_data")))

  return(data_out)
}
