#' Plot genes of interest
#'
#' Creates boxplots showing p values for significantly different (DE) genes between
#' groups.
#'
#' @param norm_exp_data Whole data normalisation output - use VST norm.
#' @param DE_full_out Full output of [auto_generate_DE_results].
#' @param se_data SummarizedExperiment object of experimental data.
#' @param selected_genes Character vector of gene names to plot.
#' @param gene_annotations Output from call to [annotate_gene_ensembl].
#' @param top_level_colname non-string. Name of \code{colData} column with top
#'   level factor, which is normally "Tissue_Region".
#' @param sample_colname non-string. Name of \code{colData} column with
#'   sample names (used for filtering).
#' @param results_contrast_factor Name of Factor from
#'   \code{DESeq2_formula_design} that is used to generate pairwise comparisons
#' @param top_level_groups character vector of names of each element of \code{top_level_colname}
#'   to filter DE data by before plotting.
#' @param plot_groups_colname non-string. Name of \code{colData} column with
#' group data for grouping data when plotting.
#' @param max_genes_per_plot number. This will split the plots so that only n number of genes are displayed per plot. Defaults to 3.
#' @param x_label string. label for x axis. Normally similar to \code{plot_groups_colname}.
#' @param p_thresh P value threshold for plotting significance lines.
#'
#'
#'
#' @return Returns a ggplot.
#'
#' @export
#'
plot_genes_of_interest <-
  function(norm_exp_data,
           DE_full_out,
           se_data,
           gene_annotations,
           selected_genes,
           top_level_colname,
           sample_colname,
           results_contrast_factor,
           top_level_groups,
           plot_groups_colname,
           max_genes_per_plot = 3,
           x_label = "Treatment",
           p_thresh = 0.05){

    tln <- rlang::enquo(top_level_colname)
    tln_name <- tln %>% rlang::as_name()
    rcf <- rlang::enquo(results_contrast_factor) %>% rlang::as_name()
    pgc <- rlang::enquo(plot_groups_colname) %>% rlang::as_name()

    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      stop("Package \"ggpubr\" needed for plot_genes_of_interest function to work. Please install it using install.packages(\"ggpubr\").",
           call. = FALSE)
    }
    if (!requireNamespace("lemon", quietly = TRUE)) {
      stop("Package \"lemon\" needed for plot_genes_of_interest function to work. Please install it using install.packages(\"lemon\").",
           call. = FALSE)
    }
    if (!requireNamespace("viridis", quietly = TRUE)) {
      stop("Package \"viridis\" needed for plot_genes_of_interest function to work. Please install it using install.packages(\"viridis\").",
           call. = FALSE)
    }


    ############################################################# #
    # Access DE tables, and reformat
    # Need FULL table, not just DE
    DE_res_obj <- GET_DESeq2_res_object(DE_full_out, top_level_groups)

    DE_data_list <-
      purrr::map(DE_res_obj,
               function(.x, annot){
                 .x %>%
                   BiocGenerics::as.data.frame() %>%
                   tibble::rownames_to_column(var = "gene_ensembl") %>%
                   dplyr::left_join(annot, by = c("gene_ensembl")) %>%
                   dplyr::select(.data$gene_ensembl,
                                 .data$gene_name,
                                 .data$description,
                                 tidyselect::everything())
                 },
               annot = gene_annotations)

    DE_data <-
      DE_data_list %>%
      dplyr::bind_rows(.id = "comparison") %>%
      dplyr::mutate(pairwise_comparison = stringr::str_split(.data$comparison, pattern = "\\.", n = 2, simplify = TRUE)[,2],
                    top_level = stringr::str_split(.data$comparison, pattern = "_", n = 2, simplify = TRUE)[,1])


    ############################################################# #
    # Access coldata
    colData <-
     se_data %>%
      SummarizedExperiment::colData() %>% as.data.frame()

    #coldata for renaming groups for pvalues
    new_groups <- colData %>%
      dplyr::select(!!rcf,
                    !!pgc) %>%
      dplyr::distinct({{results_contrast_factor}}, .keep_all = TRUE)


    #prepare genes into shorter lists
    genes_to_plot <-
      split(selected_genes,
          ceiling(seq_along(selected_genes) / max_genes_per_plot))



    .plot_genes <- function(genes_to_plot){
      ############################################################# #
      # Prepare normalised data for plotting
      norm_long <-
        norm_exp_data %>%
        dplyr::filter(.data$gene_name %in% genes_to_plot) %>%
        dplyr::mutate(gene_name = factor(.data$gene_name, levels = genes_to_plot)) %>%  #set as factor and relevel
        dplyr::select(.data$gene_name, .data$gene_ensembl, .data$description, tidyselect::contains(top_level_groups)) %>%
        tidyr::pivot_longer(cols = tidyselect::contains("__"),
                            names_to = rlang::enquo(sample_colname) %>% rlang::as_name(),
                            values_to = "Value") %>%
        dplyr::left_join(colData) %>%
        dplyr::mutate(column_to_plot = .data$Value)

      ############################################################# #
      # Create dataframe with maximum values
      # This is joined to following dataframe to help with plotting
      max_for_y.position <- norm_long %>%
        dplyr::group_by(.data$gene_name, {{ top_level_colname }}) %>%
        dplyr::summarise(y.position = max(.data$column_to_plot)) %>%
        dplyr::rename("top_level" = {{ top_level_colname }})


      ############################################################# #
      #Get pvalues data for adding to the ggplot
      data_DE1 <-  DE_data %>%
        dplyr::filter(.data$gene_name %in% genes_to_plot) %>%
        dplyr::mutate(gene_name = factor(.data$gene_name, levels = genes_to_plot)) %>%  #set as factor and relevel for plotting
        dplyr::filter(.data$top_level %in% top_level_groups)

      data_DE_with_groups <- data_DE1 %>%
        dplyr::mutate(group1 = stringr::str_split(.data$pairwise_comparison, pattern = " vs ", simplify = T)[,1],
                      group2 = stringr::str_split(.data$pairwise_comparison, pattern = " vs ", simplify = T)[,2],
                      sort = dplyr::case_when(.data$padj <= p_thresh ~ "sig",
                                       TRUE ~ "ns"),
                      padj = dplyr::case_when(.data$padj < 0.05 ~ format(.data$padj, scientific = T, digits = 1),
                                       TRUE ~ format(round(.data$padj, digits = 3), scientific = F))
        )

      #define new group names
        new_group1 <- data_DE_with_groups %>%
        dplyr::select(!!rcf := .data$group1) %>%
        dplyr::left_join(new_groups) %>%
          dplyr::select({{ plot_groups_colname }})

      new_group2 <- data_DE_with_groups %>%
        dplyr::select(!!rcf := .data$group2) %>%
        dplyr::left_join(new_groups) %>%
        dplyr::select({{ plot_groups_colname }})

      #add groups, yposition max and filter only significant values.
      data_DE <-
      data_DE_with_groups %>%
        dplyr::select(-.data$group1, -.data$group2) %>%
        dplyr::mutate(group1 = new_group1[[1]],
                      group2 = new_group2[[1]]) %>%
        dplyr::left_join(max_for_y.position) %>%
        dplyr::filter(.data$sort == "sig") %>%
        dplyr::select(.data$gene_name, .data$group1, .data$group2, p.adj = .data$padj, .data$y.position, .data$top_level)


      #calculate how many pvalues need ploting per facet, helps to space out in following function
      n_row_data <- data_DE %>%
        dplyr::group_by(.data$gene_name, .data$top_level) %>%
        dplyr::summarise(n_rows = dplyr::n())

      if(nrow(n_row_data) == 0){
        gene_p_data <-
          data_DE %>%
          dplyr::arrange(.data$top_level, .data$gene_name) %>%
          dplyr::mutate(!!tln_name := .data$top_level)
      } else {
        #recalculate y.position based on manual transformation
        gene_p_data <-
          data_DE %>%
          dplyr::left_join(n_row_data) %>%
          dplyr::group_by(.data$top_level, .data$gene_name) %>%
          dplyr::mutate(y.position_2 = (1:.data$n_rows[1])*(.data$y.position/8)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(y.position = .data$y.position+.data$y.position_2) %>%
          dplyr::arrange(.data$top_level, .data$gene_name) %>%
          dplyr::mutate(!!tln_name := .data$top_level) #to match for plotting
      }





      #add a y_max0 value to help with plotting a geom_blank() to increase the y-axis height diffrently for each facet
      data_to_plot <-
        norm_long %>%
        dplyr::left_join(dplyr::select(gene_p_data, .data$gene_name, .data$top_level, .data$y.position), by = "gene_name") %>%
        dplyr::group_by(.data$top_level, .data$gene_name) %>%
        dplyr::mutate(y_max0 = max(.data$y.position)*1.05) %>%
        dplyr::select(-.data$y.position) %>%
        dplyr::distinct()


      p_for_publish<-
        ggplot2::ggplot(data = data_to_plot,
                        mapping = ggplot2::aes(x = {{ plot_groups_colname }},
                                               y = .data$column_to_plot,
                                               colour = {{ plot_groups_colname }})) +
        ggplot2::geom_blank(ggplot2::aes(x= {{ plot_groups_colname }}, y = .data$y_max0))+
        ggplot2::geom_boxplot(ggplot2::aes(fill = {{ plot_groups_colname }}))+
        ggplot2::geom_jitter(ggplot2::aes(fill = {{ plot_groups_colname }}),
                             colour = "black",
                             size = 1,
                             width = 0.3,
                             alpha = 0.3,
                             shape = 21)+
        viridis::scale_color_viridis(discrete=TRUE)+
        viridis::scale_fill_viridis(discrete=TRUE, alpha = 0.5)+
        lemon::facet_rep_grid(cols = ggplot2::vars(!!tln), rows = ggplot2::vars(.data$gene_name), repeat.tick.labels ='y', scale = "free_y")+
        ggplot2::xlab(x_label)+
        ggplot2::ylab("Log2 Expression")+
        #ylim(c(NA,NA))+
        ggpubr::theme_classic2(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6, colour = "black", angle=25, hjust=1),
                       strip.text.x = ggplot2::element_text(size = 8, colour = "black", face = "plain"),
                       strip.text.y = ggplot2::element_text(size = 8, colour = "black", face = "italic"), #or "bold.italic"
                       axis.title = ggplot2::element_text(size = 10, colour = "black", face = "plain"),
                       legend.title = ggplot2::element_text(size = 10, colour = "black", face = "plain"))+
        ggpubr::stat_pvalue_manual(data = gene_p_data ,
                                   y.position = 'y.position',
                                   label = "adj-P = {p.adj}",
                                   #step.increase = 0.12,
                                   label.size = 2.5,
                                   color = "blue",
                                   #step.group.by = rlang::enquo(top_level_colname) %>% rlang::as_name(),
                                   step.group.by = tln_name,
                                   hide.ns=FALSE)

      return(p_for_publish)
    }


    #apply list of genes to function
    plots_out <-
      purrr::map(genes_to_plot,
                 .plot_genes)

    return(plots_out)
  }
