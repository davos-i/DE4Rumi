% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCIT.R
\name{PCIT_prepare_data}
\alias{PCIT_prepare_data}
\title{Prepare PCIT}
\usage{
PCIT_prepare_data(
  data_type_selection = c("DE", "PIF", "RIF"),
  norm_exp_data,
  DE_full_out,
  RIF_output_data,
  gene_annotations,
  search_top_level = NULL,
  filter_norm_data_columns = TRUE,
  DEn,
  PIF_sig,
  PIFn,
  RIF_sig,
  RIFn,
  export_tables = FALSE,
  export_venn_plot = FALSE,
  export_dir = "./outputs/"
)
}
\arguments{
\item{data_type_selection}{A character vector of one or any of: "DE", "PIF",
"RIF".}

\item{norm_exp_data}{Whole data normalisation output - use VST norm.}

\item{DE_full_out}{Full output of \link{auto_generate_DE_results}.}

\item{RIF_output_data}{Full output of \link{calculate_RIF}, if "RIF" included in
\code{data_type_selection}.}

\item{gene_annotations}{Output from call to \link{annotate_gene_ensembl}}

\item{search_top_level}{character vector of top_level_groups (see
\link{auto_generate_DE_results}) to search for, returning only tables from these
groups. e.g. \code{c("ARC","LIV")}. If not provided it will return all top
level groups in output.}

\item{filter_norm_data_columns}{logical. If \code{TRUE}, the
\code{norm_exp_data} will be filtered to only include samples that are
included in the top_level_groups selected by \code{search_top_level}. If
\code{TRUE} and {search_top_level} is not provided or NULL, it will return
only columns for top_level_groups that are within the DE_full_out.
Sometimes a user may want to select genes on criteria from only a certain
group, but still calculate PCIT across all normalised data.}

\item{DEn}{Number of values to return, ranked on adjusted p value. E.g. if
\code{DEn = 20}, then the 20 genes with the lowest p value are returned
(there may be less than this due to prior filtering for significance).}

\item{PIF_sig}{Number. Threshold used for filtering the z-scored PIF values,
so a \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas
a value of 1.96 filters for a nominal p value of 0.05.}

\item{PIFn}{Number of values to return for both the max PIF and min PIF. E.g.
if \code{PIFn = 20}, then the 20 highest and 20 lowest are returned (if
that many are present after filtering by PIF_sig)}

\item{RIF_sig}{Number. Threshold used for filtering the z-scored PIF values,
so a \code{PIF_sig} of 2.58 filters for a nominal P value of 0.01, whereas
a value of 1.96 filters for a nominal p value of 0.05.}

\item{RIFn}{Number of values to return for both the max PIF and min PIF. E.g.
if \code{PIFn = 20}, then the 20 highest and 20 lowest are returned (if
that many are present after filtering by PIF_sig)}

\item{export_tables}{logical. \code{TRUE} indicates the tables of selected
genes and filtered normalised counts tables are exported to
\code{export_dir}. This is exported in formats useful for input to
Cytoscape software, or for running PCIT using FORTRAN code (or windows
executable) instead of in R.}

\item{export_venn_plot}{logical. \code{TRUE} indicates a venn diagram
summarising the catagories of selected genes is exported to
\code{export_dir}. This can only be exported, not viewed interactively
within RStudio.}

\item{export_dir}{string. "./" indicates relative to working directory. IF
this directory doesn't exist, it will be created. Defaults to "./outputs/".}
}
\value{
Returns a list of dataframes of filtered output for all selected data
types, as well as a normalised counts table filtered to include only
selected genes, for use in PCIT.
}
\description{
Function to prepare various data sources for selection for a Partical
Corellation Information Theory (PCIT) analysis. It is not necessarily useful
to use all DE, significant PIF or significant RIF values if there are many
hundreds or thousands as input to PCIT. This function will filter for
significance (PIF and RIF) and return a number of genes, based on user input.
}
