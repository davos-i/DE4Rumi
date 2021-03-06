% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESEQ_helpers.R
\name{check_count_matrix}
\alias{check_count_matrix}
\title{Check formatting of supplied count matrix}
\usage{
check_count_matrix(count_data, colData, column_with_col_names)
}
\arguments{
\item{count_data}{Dataframe of counts data to check.}

\item{colData}{Dataframe of column annotation data. each column should be
metadata about the sample/animal. Particularly include all treatment info
and animal ID in columns.}

\item{column_with_col_names}{non-string. Name of the column in \code{colData}
that matches names of columns in \code{count_data}. This column hould not
include and entry called "gene_ensembl", as this is automatically detected.}
}
\value{
Returns a dataframe of \code{count_data} with columns ordered to
match names in the \code{column_with_col_names} column of \code{colData}.
}
\description{
\code{check_count_matrix()} returns messages to the console to determine if
count matrix has been formatted correctly for use in downstream analysis, and
to see if column names match colData (metadata for each column).
}
\details{
The count matrix, normally as outputted from featureCounts, needs to be
formatted with:\itemize{ \item Column headings that match metadata \item No
rownames \item First column titled gene_ensembl and contains all gene ensembl
names }
This function also prints the dimensions and predicted treatments to
console.
}
