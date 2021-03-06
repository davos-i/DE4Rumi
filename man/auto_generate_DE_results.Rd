% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auto_generate_DE_results.R
\name{auto_generate_DE_results}
\alias{auto_generate_DE_results}
\title{Auto Generate DE Results}
\usage{
auto_generate_DE_results(
  se_data,
  top_level_groups,
  top_level_colname,
  sample_colname,
  samples_to_remove = NA,
  DESeq2_formula_design,
  rowSums_filter = 10,
  results_contrast_factor,
  results_combinations = NA,
  use_IHW_filtering = TRUE,
  alpha = 0.05,
  gene_annotations,
  export_tables = FALSE,
  export_dir = "./outputs/",
  whole_data_normalisation = FALSE
)
}
\arguments{
\item{se_data}{SummarizedExperiment object of experimental data.}

\item{top_level_groups}{character vector of names of each element of \code{top_level_colname}
to be iterated over.}

\item{top_level_colname}{non-string. Name of \code{colData} column with top
level factor, which is normally "Tissue_Region".}

\item{sample_colname}{non-string. Name of \code{colData} column with
sample names (used for filtering).}

\item{samples_to_remove}{string/s. Names of samples to remove e.g.
\code{c("ARC_01_HCP-HP-RMEI", "ARC_02_HCP-HP-UMEI")}}

\item{DESeq2_formula_design}{= formula in form \code{~Factor} parsed to
\code{DESeq()}.}

\item{rowSums_filter}{numeric. For intial filter of rowSums across whole dds.
Filters out all samples that have basically 0 expression across all
samples. Saves computing time, but does not change DESeq results as it
filters internally as well. Defaults to 10.}

\item{results_contrast_factor}{Name of Factor from
\code{DESeq2_formula_design} that is used to generate pairwise comparisons.}

\item{results_combinations}{Either \code{NA}, which will automatically
generate all possible pairwise combinatios of the
\code{results_contrast_factor} OR a list of pairwise comparisons with each
element containing 2 strings, numerator and denominator, respectively. Should
match format of output from \link{make_pairwise_combinations}. Defaults to \code{NA}.}

\item{use_IHW_filtering}{logical. Should Independent Hypothesis Weighting be
used when calculating results? \code{TRUE}, the default, indicates yes. See
\link[IHW:ihw.default]{IHW::ihw} or \code{vignette("DESeq2")} for more details.}

\item{alpha}{numeric. What is the p-value threshold to be used for
determining significance. Used in call to \link[DESeq2:results]{DESeq2::results} and
others.}

\item{gene_annotations}{Output from call to \link{annotate_gene_ensembl}}

\item{export_tables}{logical. \code{TRUE} indicates the normalised counts
tables, both vsd and log2norm, annotated with gene names and descriptions,
are exported to \code{export_dir}}

\item{export_dir}{string. "./" indicates relative to working directory. IF
this directory doesn't exist, it will be created.}

\item{whole_data_normalisation}{logical. Should DESEQ2 be run for pairwise
comparisons, split by \code{top_level_groups}? or, if FALSE - normalised
output only for whole dataset.}
}
\value{
returns a named, nested list (list of lists) with: \itemize{
\item dds_wald_object \item boxplot_cooks_distance \item DESeq2_res_object
\item pairwise_plots \item overall_plots \item normalised_data \item PIF
\item DE_by_PIF_df}
}
\description{
Wrapper to DESeq2 that iterates over each element of \code{top_level_groups}
and conducts pairwise comparisons, QC plots and other analysis. Can also
normalise data over whole dataset (ignoring \code{top_level_groups})
}
\details{
This wrapper determines all possible pairwise comparisons, based on the
\code{top_level_colname} attribute which is also a column of \code{colData}.
Defaults to use a formula design of \code{~Region_Diet}, but could be changed
to be \code{~Region + Diet} (see \code{vignette("DESeq2")}) but the pairwise
comparisons are generated with \code{results_contrast_factor}, which must be
one of the elements of the formula.
}
