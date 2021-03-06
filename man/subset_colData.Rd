% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESEQ_helpers.R
\name{subset_colData}
\alias{subset_colData}
\title{Subset colData automatically}
\usage{
subset_colData(count_data, colData, column_with_col_names)
}
\arguments{
\item{count_data}{Dataframe of counts data to check.}

\item{colData}{Dataframe of column annotation data. each column should be
metadata about the sample/animal. Particularly include all treatment info
and animal ID in columns.}

\item{column_with_col_names}{non-string. Name of the column in \code{colData}
that matches names of columns in \code{count_data}.}
}
\value{
Returns the filtered dataframe (to be used downstream).
}
\description{
Run if \link{check_count_matrix} warns that there are entries in colData
that do not match columns in counts_matrix.
}
\details{
Matches the entries in the \code{column_with_col_names} in \code{colData} to
the columns of \code{count_data}, and removes them from \code{colData}.
}
