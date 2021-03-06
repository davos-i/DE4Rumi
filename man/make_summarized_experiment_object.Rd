% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESEQ_helpers.R
\name{make_summarized_experiment_object}
\alias{make_summarized_experiment_object}
\title{Make SE object}
\usage{
make_summarized_experiment_object(counts_data, gene_annotations, colData)
}
\arguments{
\item{counts_data}{Dataframe of gene expression counts data.}

\item{gene_annotations}{Dataframe with columns "gene_ensembl", "gene_name" and
"description", returns from a call to \code{\link{annotate_gene_ensembl}}}

\item{colData}{Dataframe of column annotation data. each column should be
metadata about the sample/animal.}
}
\value{
Returns a SummarizedExperiment object of experimental data.
}
\description{
Takes gene expression count data, gene annotations and column data and returns
a SummarizedExperiment object to use with
\code{\link{auto_generate_DE_results}}
}
\details{
\code{counts_data} and \code{colData} inputs to this function should be
pre-checked using \code{\link{check_count_matrix}}
}
