% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESEQ_helpers.R
\name{make_pairwise_combinations}
\alias{make_pairwise_combinations}
\title{Make pairwise combinations}
\usage{
make_pairwise_combinations(coldata, contrast_factor)
}
\arguments{
\item{coldata}{either a dataframe or a call to \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment::colData}
on a DESeqDataSet}

\item{contrast_factor}{string. Name of column containing treatments to contrast}
}
\description{
Creates a list of all possible unique
combinations of treatment comparisons
}
\details{
This takes the colData from a dds object and finds the column name defined by
\code{contrast_factor} to generate all unique treatment names as input to a call
to \link[utils:combn]{utils::combn} to get all unique pairwise comparisons.
Returns a list of characters of length 2. \cr
e.g. \code{[1] "LIV_HCP-HP-UMEI" "LIV_LCP-LP-UMEI"}
}
