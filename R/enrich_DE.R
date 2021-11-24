#'Functional Enrichment of DE lists
#'
#'Identify GO, REACTOME and KEGG pathways for DE lists via gProfiler's API
#'
#'\code{enrich_DE()} takes the DE tables from a call to
#'\code{auto_generate_DE_results} as input to a call to
#'\code{gprofiler2::gost()}. Various custom options include: \itemize{\item
#'Filtering by adjusted p value or log2 fold change \item Split the DE list into
#'positive and negative fold change, to enrich separately. \item Remove parent terms from list
#'to de-clutter redudant terms. }
#'
#'@param list_DE_obj list of DESeq2 results objects, normally accessed via a
#'  call to \code{auto_generate_DE_results()}.
#'@param sig_thresh number. Threshold to filter results by adjusted P value.
#'  Defaults to 0.05.
#'@param log2_thresh number. Optional filter - if there is too many DE genes,
#'  then list is filtered by log2 fold change. Will keep only DE genes with a FC
#'  above this threshold. Uses absolute values (includes + and -).
#'@param split_by_pos_neg logical. Should the query be run separately for genes
#'  with a + fold change, compared to genes with a - fold change?
#'@param gene_annotations Output from call to \code{annotate_gene_ensembl()}
#'@param organism string. From gprofiler2: Organism names are constructed by
#'  concatenating the first letter of the name and the family. For sheep data:
#'  "oaries". For cattle data:"btaurus". Default is "oaries".
#'@param base_URL string. URL for the version of g:Profiler to use. This
#'  defaults to e100 as it matches the galaxy workflow, but also e103 broke a
#'  lot of names and is yet unresolved. Provided in the format:
#'  "http://biit.cs.ut.ee/gprofiler_archive3/e100_eg47_p14". \cr\cr Note that it
#'  is "http://" NOT "https://". \cr\cr Providing NA will default to g:Profiler2
#'  default of "http://biit.cs.ut.ee/gprofiler". \cr\cr All available archives
#'  at: https://biit.cs.ut.ee/gprofiler/page/archives
#' @param export_tables logical. \code{TRUE} indicates the normalised counts
#'   tables, both vsd and log2norm, annotated with gene names and descriptions,
#'   are exported to \code{"./outputs/normalised_counts/"}
#' @param export_dir string. "./" indicates relative to working directory. IF
#'   this directory doesn't exist, it will be created.
#'
#'@return Returns a list of tables of same length as \code{list_DE_obj}, each with
#'  g:Profiler results
#'
#'@export

enrich_DE <-
  function( list_DE_obj,
            sig_thresh = 0.05,
            log2_thresh = 0,
            split_by_pos_neg = FALSE,
            gene_annotations,
            organism,
            base_URL = "http://biit.cs.ut.ee/gprofiler_archive3/e100_eg47_p14",
            export_tables = FALSE,
            export_dir = "./outputs/"){

    message(crayon::green("Filtering data..."))
    DE_lists <-
      purrr::map(list_DE_obj,
                 function(x, sig_thresh, log2_thresh){
                   data <- x %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column("gene_ensembl") %>%
                     dplyr::filter(.data$padj <= sig_thresh & abs(.data$log2FoldChange)>log2_thresh) %>%
                     dplyr::arrange(.data$padj) %>%
                     dplyr::mutate(log2FC_positive_negative =
                                     dplyr::case_when(.data$log2FoldChange < 0 ~ "Neg",
                                                      .data$log2FoldChange > 0 ~ "Pos"))
                   if(split_by_pos_neg == TRUE){
                     #split by Pos Neg
                     return(data %>% split(magrittr::use_series(data, "log2FC_positive_negative")))
                   } else{
                     return(data)
                   }

                 },
                 sig_thresh,
                 log2_thresh
      )

    #collapse list
    if(split_by_pos_neg == TRUE){
      DE_lists <- DE_lists %>% unlist(recursive = FALSE)
    }
    message(crayon::green("Filtering data... COMPLETE"))

    #select only the 'ensembl' column from each dataframe in the list
    ensembl_named_list <- purrr::map(DE_lists, 'gene_ensembl')

    if(is.na(base_URL)){
      gprofiler2::set_base_url("http://biit.cs.ut.ee/gprofiler")
      message(paste("Using default g:Profiler URL: ", gprofiler2::get_base_url()))
    } else {
      gprofiler2::set_base_url(base_URL)
      message(paste("g:Profiler Version URL: ", gprofiler2::get_base_url()))
    }

    #run g:Profiler - Returns one large dataframe, if multi_query is set to TRUE
    #then the intersection column is not produced.
    message(crayon::green("Running gprofiler2::gost()..."))

    gost_result1 <- gprofiler2::gost(query = ensembl_named_list,
                                      organism = organism,
                                      multi_query = FALSE,
                                      evcodes = T)
    message(crayon::green("Running gprofiler2::gost()... COMPLETE"))

    ############################################################# #
    #### annotate the 'intersection' column
    message(crayon::green("Reformatting dataframe..."))
    gost_result_new <-
      gost_result1$result %>% #df or results
      #this function takes a column in a dataframe that is a list, and collapses them into a printed list. Applied to each row (see specific use case of pmap).
      purrr::pmap_dfr(
        function(...){
          current <- tibble::tibble(...)
          parents_new <- current$parents %>% stringr::str_c(collapse = ", ")
          current[1,] %>% dplyr::mutate(parents = parents_new)
        }
      ) %>%
      #Re-name the genes in the intersection to be HGNC names (or other annotation).
      #Intersection are the genes from DE list that matched the GO term
      purrr::pmap_dfr(function(...){
        current_iteration <- data.frame(...) #pmap takes '...' (all arguments - in this case all column names), this line then puts the current 'row' data together
        #Then, the column name of interest can be selected, as all columns are available in 'current_iteration'
        current_iteration$intersection_gene_names <-
          current_iteration$intersection %>%  #take the ensembl IDS that were matched by the GO term
          stringr::str_split(pattern = ",") %>% #split into individual items
          as.data.frame(col.names = "gene_ensembl") %>%
          dplyr::left_join(gene_annotations, by = "gene_ensembl") %>%
          magrittr::use_series("gene_name") %>%
          stringr::str_c(collapse = ",")

        #Truncate long strings within each 'cell' (for export to excel later)
        current_iteration %>%
          dplyr::mutate(dplyr::across(where(is.character),
                                      function(x){
                                        stringr::str_trunc(x, width = 20000, side = "right")
                                      })
          )
      })

    message(crayon::green("Reformatting dataframe... COMPLETE"))
    ############################################################# #
    # Split by query column
    message(crayon::green("Split dataframe... "))

    gost_result_split <- split(gost_result_new, gost_result_new$query)

    message(crayon::green("Split dataframe... COMPLETE"))


    ############################################################# #
    # rename to truncate to fit excel worksheet requirements
    gost_result_split2 <- gost_result_split

    names1 <- names(gost_result_split2) %>%
      stringr::str_split(pattern = "\\.", n = 2, simplify = TRUE) %>%
      magrittr::extract(,2) %>%
      stringr::str_squish() %>%
      stringr::str_remove_all(" ") %>%
      stringr::str_trunc(27,ellipsis = "")

    #truncated but with Pos or neg at start.
    if(split_by_pos_neg == TRUE){
      pos_neg_names <- names(gost_result_split2) %>%
        stringr::str_extract(pattern = "\\.(.{3})$") %>%
        stringr::str_remove(pattern = "\\.")

      names(gost_result_split2) <-
        stringr::str_c(pos_neg_names, names1, sep = ".")
    } else{
      names(gost_result_split2) <- names1
      }

    ############################################################# #
    # Export
    if(export_tables == TRUE){

      if(!dir.exists(export_dir)){
        dir.create(export_dir, recursive = TRUE)
        message(crayon::red(paste("Directory created:", export_dir)))
      } else{message(crayon::green(paste(export_dir,"Directory exists")))}


      openxlsx::write.xlsx(gost_result_split2,
                           file = paste0(export_dir,"gProfiler functional enrichment tables - SPLIT - ",
                                         format(Sys.time(), "%Y%m%d_%H%M"),
                                         ".xlsx"),
                           colWidths = "auto")

      message(crayon::green(paste("gProfiler functional enrichment tables exported to an .xlsx file in the sub-directory:", export_dir)))

    }

    return(gost_result_split)
  }

