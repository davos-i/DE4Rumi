% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genes_of_interest.R
\name{plot_genes_of_interest}
\alias{plot_genes_of_interest}
\title{Plot genes of interest}
\usage{
plot_genes_of_interest(
  norm_exp_data,
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
  p_thresh = 0.05
)
}
\arguments{
\item{norm_exp_data}{Whole data normalisation output - use VST norm.}

\item{DE_full_out}{Full output of \link{auto_generate_DE_results}.}

\item{se_data}{SummarizedExperiment object of experimental data.}

\item{gene_annotations}{Output from call to \link{annotate_gene_ensembl}.}

\item{selected_genes}{Character vector of gene names to plot.}

\item{top_level_colname}{non-string. Name of \code{colData} column with top
level factor, which is normally "Tissue_Region".}

\item{sample_colname}{non-string. Name of \code{colData} column with
sample names (used for filtering).}

\item{results_contrast_factor}{Name of Factor from
\code{DESeq2_formula_design} that is used to generate pairwise comparisons}

\item{top_level_groups}{character vector of names of each element of \code{top_level_colname}
to filter DE data by before plotting.}

\item{plot_groups_colname}{non-string. Name of \code{colData} column with
group data for grouping data when plotting.}

\item{max_genes_per_plot}{number. This will split the plots so that only n number of genes are displayed per plot. Defaults to 3.}

\item{x_label}{string. label for x axis. Normally similar to \code{plot_groups_colname}.}

\item{p_thresh}{P value threshold for plotting significance lines.}
}
\value{
Returns a ggplot.
}
\description{
Creates boxplots showing p values for significantly different (DE) genes between
groups.
}
