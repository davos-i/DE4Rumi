#' Plot all MA plots automatically
#'
#' Plot all MA plots automatically by iterating over a list of PIF dataframes,
#' naming each plot
#'
#' Wrapper for [plot_custom_MA] that uses sensible defaults and naming,
#' but is customisable.
#'
#'
#' @param PIF_list List of dataframes produced by a call to
#'   [auto_generate_DE_results]
#' @param list_plot_titles A character vector or list of names in the same order
#'   as \code{names(PIF_list)}. If NA (default) names are automatically
#'   generated from  \code{names(PIF_list)}.
#' @param export_TIFFs Logical. Should individual TIFF images be exported for
#'   each MA plot?
#' @param export_pdf Logical. Should all plots be exported to one pdf file?
#' @param pdf_height number. Height of each pdf page in mm.
#' @param pdf_width number. Width of each pdf page in mm.
#' @param export_dir string. Path directory where TIFF and PDF are exported to,
#'   if exported.
#' @param ... Other arguments parsed to [plot_custom_MA]. \cr See
#'   [plot_custom_MA] for details of available arguments.
#'
#' @return Returns a named list of ggplot objects. Also exports plots to .TIFF
#'   or .PDF, if requested.
#'
#' @export
#'
#'
auto_plot_multi_MA <-
  function(PIF_list,
           list_plot_titles = NA,
           export_TIFFs = FALSE,
           export_pdf = FALSE,
           pdf_height = 140,
           pdf_width = 220,
           export_dir = "./outputs/MA_plots/",
           ...){


    if(any(is.na(list_plot_titles))){
      list_of_names <- names(PIF_list) %>%
      stringr::str_split(pattern = "-", n = 2, simplify = TRUE) %>%
      magrittr::extract(,2) %>%
      stringr::str_trim()

    list_plot_titles <- paste("MA Plot of", list_of_names)

    } else{message(crayon::red("Using custom titles from list_plot_titles"))}


    .plot_custom_MA <-
      function(.x, .y, ...){
        plot_custom_MA(PIF_data = .x,
                    plot_title = .y,
                    export_TIFFs = export_TIFFs,
                    export_dir = export_dir,
                    ...)
    }


    plots_out <- purrr::map2(.x = PIF_list,
                .y = list_plot_titles,
                .f = .plot_custom_MA,
                ...)

   plots_out <- setNames(plots_out, list_plot_titles)

   if(export_pdf == TRUE){
     message(crayon::green("Plotting MA plots to pdf file..."))

     ggplot2::ggsave(
       filename = paste0("MA_PLOTS_",format(Sys.time(), "%Y%m%d_%H%M"), ".pdf"),
       plot = gridExtra::marrangeGrob(plots_out, nrow=1, ncol=1),
       path = export_dir,
       width = pdf_width, height = pdf_height, units = "mm"
     )
     message(crayon::green(paste0("MA plot pdf exported to: ", export_dir)))
   }

   return(plots_out)
   }
