% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncbitax.R
\name{read.ncbi}
\alias{read.ncbi}
\title{Load NCBI Taxonomy Tree}
\usage{
read.ncbi(dir, use.cache = TRUE, clear.cache = FALSE, env.cache = TRUE)
}
\arguments{
\item{dir}{Optional directory containing the unpackaged NCBI taxonomy dump
files names.dmp, nodes.dmp and merged.dmp. If not provided, files will be
downloaded using \code{\link{download.ncbitax}} automatically.}

\item{use.cache}{Since downloading and parsing the dumps from NCBI is time
consuming, read.ncbi will cache the final data.table in the file
"\code{ncbitax.Rds}" in the CWD. If this parameter is set to \code{FALSE},
an existing cache file will be ignored.}

\item{clear.cache}{Delete an existing "ncbitax.Rds" before proceeding to
prepare data structure if set to \code{TRUE}.}

\item{env.cache}{If \code{TRUE} (default), \code{read.ncbi} will try to
return the taxonomy object cached from its last invocation.}
}
\value{
Returns a data.table containing NCBI taxonomy. The table is also
  cached internally for use with the functions within this package.
}
\description{
Prepares the  NCBI taxonomy for use with the other functions in \link{taxme}.
}
\examples{

ncbi <- read.ncbi()

}
\seealso{
\code{\link{ncbi.lineage}}
}
