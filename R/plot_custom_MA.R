
#' MA Plots
#'
#' Create labelled MA plots from previously calculated DE and PIF data.
#'
#' These MA plots show the mean expression between Condition 1 and Condition 2
#' along the x axis, and the log2 fold change on the y axis. \cr Points that are
#' both highly abundant and highly different (higher/lower log2 fold change) are
#' considered to be most likely to be biologically relevant, and are represented
#' by the PIF calculation (see \code{calculate_PIF()}). \cr The points represent
#' the averages across multiple individual samples/animals. Therefore, points
#' that appear highly different may not be statistically different due to high
#' variance between the underlying points. Therefore, points that are both
#' Differentially Expressed (statistically different) and PIF are likely the
#' most important points. \cr Labels are generated automatically, and avoid too
#' many overlaps, see \code{?ggrepel::ggrepel}. Increasing the plot size will
#' increase number of labels. In RStudio, this is as simple as clicking on Zoom
#' in the 'Plots' window, or changing the \code{fig.width} and/or
#' \code{fig.length} parameters in a RMarkdown/RNotebook code chunk. \cr\cr
#' These plots are built with ggplot2. See https://ggplot2.tidyverse.org/ for
#' more details on the specs used.
#'
#'
#'
#' @param PIF_data dataframe with at least the following columns: starts_with...
#'   "gene_ensembl", "mean_", "log2FC_" and "zPIF"
#' @param padj_DE_filter Number. Significance threshold for determining if
#'   points are labelled as DE or not.
#' @param PIF_sig Number. Used for filtering the z-scored PIF values, so a
#'   \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas a
#'   value of 1.96 filters for a nominal p value of 0.05.
#' @param n_labels_top_bottom Number. How many labels should be displayed on
#'   both the top and bottom. E.g. if n = 20, then 20 labels will be plotted
#'   above y = 0 and 20 below y = 0.
#' @param label_type a string of either "PIF & DE" or "PIF". If "PIF & DE" is
#'   selected, labels will only plot points that are BOTH PIF and DE, whereas
#'   "PIF" will allow any point that is significantly PIF to be labelled.
#' @param custom_labels character vector of gene names to force labels for e.g.
#'   c("CARTPT")
#' @param point_size Size of points on plot
#' @param point_alpha Transparency of points.
#' @param font_face String. Either "plain", "bold", "italic" or "bold.italic"
#' @param font_colour String of a colour name.
#' @param font_size Number. Size of font for labels.
#' @param label.rectangle logical. Should labels have a rectangle with white
#'   background? See difference between \code{geom_label()} and
#'   \code{geom_text()} by viewing \code{?ggplot2::geom_label}
#' @param plot_title String. Title of the plot.
#' @param xlab String. Title of x axis.
#' @param ylab String. Title of y axis.
#' @param ggtheme A function to a ggplot2 theme. E.g.
#'   \code{ggplot2::theme_classic()}
#' @param palette character vector of colour names. Will label points in same
#'   order as: \code{c("Not PIF or DE", "PIF", "PIF & DE", "DE", "Removed by
#'   DESeq2")}
#' @param export_TIFFs Logical. Should a TIFF image be exported?
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created.
#' @param tiff_width number in mm. Width of .tiff image if exported.
#' @param tiff_height number in mm. Height of .tiff image if exported.
#'
#'
#' @return returns a ggplot object
#'
#' @export
#'
#'
plot_custom_MA  <-
  function(PIF_data,
           padj_DE_filter = 0.05,
           PIF_sig = 2.58, #1.96 or 2.58
           n_labels_top_bottom = 20,
           label_type = "PIF & DE", #or "PIF"
           #label_lower_exp_limit = 3,
           custom_labels = NA,
           point_size = 1.5,
           point_alpha = 0.5,
           font_face = "plain",
           font_colour = "black",
           font_size = 8,
           label.rectangle = FALSE,
           plot_title = "TITLE",
           xlab = "Mean expression",
           ylab = "Log2 fold change",
           ggtheme = ggplot2::theme_classic(),
           palette = c("lightgrey", "orange", "red", "blue", "green"),
           export_TIFFs = FALSE,
           export_dir = "./outputs/MA_plots/",
           tiff_width = 190,
           tiff_height = 140
           ) {


    #remove any attributes, normally remnant attributes from scale function
    PIF_data[] <- lapply(PIF_data, function(x) { attributes(x) <- NULL; x })


    data1 <- PIF_data %>%
      dplyr::rename("zPIF" = tidyselect::starts_with("zPIF"),
                    "padj" = tidyselect::starts_with("padj"),
                    "log2FoldChange" = tidyselect::starts_with("log2FC"),
                    "mean_exp" = tidyselect::starts_with("mean"))

    data1$zPIF[abs(data1$zPIF) < PIF_sig] = NA

    message(crayon::green(paste("Number of significant PIF genes: ",
                                (length(which(!is.na(data1$zPIF))
                                        )
                                 )
                                )
                          )
            )
   # Create Groupings based on data
    data2 <-
      data1 %>%
      dplyr::mutate(Group = dplyr::if_else(is.na(.data$padj), #IF padj = NA then it was removed by DESeq2
                                           true = "Removed by DESeq2",
                                           false = dplyr::if_else(is.na(.data$zPIF), #then, if zPIF is NA...
                                                                  true = dplyr::if_else(.data$padj <= padj_DE_filter, #check if DE
                                                                                        true = "DE",
                                                                                        false = "Not PIF or DE",
                                                                                        missing = "Not PIF or DE"),  #when zPIF & padj == NA
                                                                  false = dplyr::if_else(.data$padj <= padj_DE_filter,
                                                                                         true = "PIF & DE",
                                                                                         false = "PIF",
                                                                                         missing = "PIF"), #when zPIF is significant but padj == NA
                                                                  missing = "ERROR")
      ),
      Group = factor(.data$Group, ordered = T, levels = c("Not PIF or DE", "PIF", "PIF & DE", "DE", "Removed by DESeq2"))
      )



    #reorder so that 'DE' / coloured dots appear above grey dots
    data3 <-
      data2 %>%
      dplyr::mutate(a_zPIF = abs(.data$zPIF)) %>%
      dplyr::arrange(dplyr::desc(.data$padj), .data$a_zPIF)
    #in data - a_zPIF is the absolute z-normalised PIF value; zPIF is not the absolute.


    labs_data0 <-   tidyr::drop_na(data3, .data$gene_name) #removes any rows with missing name (if they exist)


    ############################################################################ #
    if(label_type == "PIF"){
      # for selecting n above the line and n below the line
      labs_data <- plyr::rbind.fill(dplyr::slice_max(labs_data0,
                                              n = n_labels_top_bottom,
                                              order_by = .data$zPIF),
                                    dplyr::slice_min(labs_data0,
                                              n = n_labels_top_bottom,
                                              order_by = .data$zPIF))
    }

    # "top PIF & DE" - Only label above a certain expression level, e.g. above an expression of 3 (filters out strangely skewed data)
    # if(label_type == "top PIF & DE"){
    #   labs_data1 <-labs_data0 %>%
    #     dplyr::filter(.data$Group == "PIF & DE" |(.data$Group == "PIF" & (mean_exp ) > label_lower_exp_limit))
    #
    #   l_max <- dplyr::slice_max(labs_data1, n = n_labels_top_bottom, order_by = .data$zPIF)
    #   l_min <- dplyr::slice_min(labs_data1, n = n_labels_top_bottom, order_by = .data$zPIF)
    #
    #   labs_data <- dplyr::bind_rows(l_max, l_min)
    #
    # }

    if(label_type == "PIF & DE"){
      labs_data1 <-labs_data0 %>%
        dplyr::filter(.data$Group == "PIF & DE" )

      l_max <- dplyr::slice_max(labs_data1, n = n_labels_top_bottom, order_by = .data$zPIF)
      l_min <- dplyr::slice_min(labs_data1, n = n_labels_top_bottom, order_by = .data$zPIF)

      labs_data <- dplyr::bind_rows(l_max, l_min)

    }

    ############################################################################ #
    # Add custom labels
    if(!is.na(custom_labels)){
      labs_data <- labs_data %>%
        dplyr::bind_rows(labs_data0 %>%
                           dplyr::filter(.data$gene_name %in% custom_labels)) %>%
        dplyr::distinct(.data$gene_name, .keep_all = T)
    }


    ############################################################################ #
    # Create ggplot
    p <- ggplot2::ggplot(data3, ggplot2::aes(x = .data$mean_exp, y = .data$log2FoldChange)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$Group),
                 size = point_size,
                 alpha = point_alpha)

    ############################################################################ #
    # Add labels to points
    if(label.rectangle){
      p <- p + ggrepel::geom_label_repel(data = labs_data,
                                         mapping = ggplot2::aes(label = .data$gene_name),
                                         box.padding = ggplot2::unit(0.35, "lines"),
                                         point.padding = ggplot2::unit(0.3, "lines"),
                                         force = 1,
                                         fontface = font_face,
                                         size = font_size/ggplot2::.pt,
                                         color = font_colour,
                                         min.segment.length = 0)
    } else{
      p <- p + ggrepel::geom_text_repel(data = labs_data,
                                        mapping = ggplot2::aes(label = .data$gene_name),
                                        box.padding = ggplot2::unit(0.35, "lines"),
                                        point.padding = ggplot2::unit(0.3, "lines"),
                                        force = 1,
                                        fontface = font_face,
                                        size = font_size/ggplot2::.pt,
                                        color = font_colour,
                                        min.segment.length = 0)
    }


    ############################################################################ #
    # Identify genes that had extreme values and were not included in plots
    # extreme is classed as a log2FC > 8
    # prints names of removed genes under plot title
    out_of_bounds_genes <- data3[which(abs(data3$log2FoldChange) > 8),]$gene_name %>% paste(collapse = ", ")

    scale_limit0 <-
      data3$log2FoldChange %>%
      abs() %>%
      purrr::map(function(x){
        x[which(x < 8)] }) %>%
      unlist() %>%
      max()*1.02

    p <- p +
      ggplot2::scale_x_continuous(n.breaks = 10) +
      ggplot2::scale_y_continuous(n.breaks = 10, limits = c(-scale_limit0, scale_limit0)) +
      ggplot2::labs(x = xlab,
           y = ylab,
           title = plot_title,
           subtitle = paste('Genes out of bounds (absolute log2FC >8): ',
                            out_of_bounds_genes),
           color = "Legend")


    #set graphical params from user input in function above
    p <- ggpubr::ggpar(p, palette = palette, ggtheme = ggtheme)

    ################################################################## #
    # Add customised options for legend position
    p <- p + ggplot2::theme(#legend.position = "none",
      axis.text = ggplot2::element_text(face = 'bold', colour = "black"),
      legend.justification = c(1,0), #which point on legend is the reference c(x,y) between 0 - 1
      legend.position = c(.99, 0.02), #which point to put the legend point c(x,y) between 0 - 1
      legend.box.background = ggplot2::element_rect(color="darkgrey", size=0.5, fill = "transparent"), #colour box around legend, with transaprenty background
      legend.box.margin = ggplot2::margin(2, 2, 2, 2), #margin for legend
      legend.background = ggplot2::element_blank(), #no background to legend, full transparency
      plot.subtitle = ggplot2::element_text(size=6)
    )

    ################################################################## #
    if(export_TIFFs == TRUE){

      if(!dir.exists(export_dir)){
        dir.create(export_dir, recursive = TRUE)
        message(crayon::red(paste("Directory created:", export_dir)))
      } else{message(crayon::green(paste(export_dir,"Directory exists")))}


      filename1 <- paste0(plot_title, "_",
                          format(Sys.time(), "%Y%m%d_%H%M"),
                          ".tiff")
      ggplot2::ggsave(filename = filename1,
           plot = p,
           device = "tiff",
           path = export_dir,
           width = tiff_width, #190 mm wide for ANIMAL
           height = tiff_height,
           units = "mm",
           dpi = 1000
    )

    message(crayon::green(paste0("Saved file: ", filename1)))
    }
    ################################################################## #

    return(p)
  }
