
.norm_and_plot <- function(dds_object,
                           gene_annotations,
                           export_tables = FALSE,
                           top_level_filter){
  top_level_filter_string <- rlang::enquo(top_level_filter) %>% rlang::as_label()
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
    dplyr::select(gene_ensembl, gene_name, description, everything())



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
    dplyr::select(gene_ensembl, gene_name, description, everything())



  ################################################################################## #
  #export tables of VST and log2 normalised data, plus the full list of DE output
  if(export_tables == TRUE){
    data.table::fwrite(vsd_a ,
                       file = paste0("./outputs/normalised_counts/",top_level_filter,"_VST_normalisedcounts_",format(Sys.time(), "%Y%m%d_%H%M"),".csv"), row.names = F)

    data.table::fwrite(log2norm_new_names,
                       file = paste0("./outputs/normalised_counts/",top_level_filter,"_log2_of_DESEQ_internal_normalisedcounts_",format(Sys.time(), "%Y%m%d_%H%M"), ".txt"), row.names = F)

    # fwrite(full_dataset_DE,
    #        file = paste0("./outputs/",region,"_DE_output_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
  }


  ################################################################################## #
  #Create PCA plot
  pcaData <- DESeq2::plotPCA(vsd, intgroup=c("Diet"), returnData=TRUE, ntop = 500)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  PCA_plot <- ggplot2::ggplot(pcaData, ggplot2::aes(PC1, PC2, fill=Diet)) +
    ggplot2::geom_point(size=4.2, shape = 21, colour = "black", alpha = 0.8) +
    ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggplot2::coord_fixed()+
    #coord_fixed(xlim = c(-150,150), ylim = c(-120, 120)) +
    ggplot2::theme_classic()+
    ggplot2::geom_text(ggplot2::aes(label=name))+
    ggplot2::ggtitle(label = paste(top_level_filter_string ))


  ################################################################################## #
  #Heatmap #************** Update with sample distances pre-calculated with dist(t(assay(vsd))) - see vignette
  topVarGenes <- head( order( genefilter::rowVars( SummarizedExperiment::assay(vsd) ), decreasing=TRUE ), 500 )

  heatmap_input <- SummarizedExperiment::assay(vsd)[ topVarGenes, ]

  heatmap_plot<- pheatmap::pheatmap(heatmap_input,
                                    scale = "row",
                                    #cutree_cols = 5,
                                    #cutree_rows = 4,
                                    clustering_method = "complete",
                                    #kmeans_k = 50,
                                    #annotation_row =
                                    main = paste(top_level_filter_string ))



  ################################################################################## #
  #OUTPUTS - to lists already generated




  #put plots into a list, append to the list (from global environment) but named
  #plots of overall data
  plots <- list( heatmap = heatmap_plot, PCA = print(PCA_plot))

  #Plot_out_list <- c(Plot_out_list, setNames(list(plots), paste0(top_level_filter, "_plots")))



  #put data into a list, append to the list (from global environment) but named
  data_out <-  list(vsd = vsd_a, log2norm = log2norm_a)

  #Norm_and_DESEQ_out_list <- c(Norm_and_DESEQ_out_list, setNames(list(data_out), paste0(top_level_filter, "_data")))

  return(list(plots, data_out))
}
