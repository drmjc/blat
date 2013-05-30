\name{midpoint.bed}
\alias{midpoint.bed}
\title{Determine the midpoint of an interval}
\usage{
  midpoint.bed(x)
}
\arguments{
  \item{x}{a \code{data.frame} representation of a BED
  file}
}
\value{
  a \code{data.frame} like \code{x}, with chromStart and
  chromEnd updated to reflect the midpoint.
}
\description{
  Determine the midpoint of an interval
}
\examples{
bed <- data.frame(chrom="chr1", chromStart=c(100,200,300), chromEnd=c(110, 250, 500))
midpoint.bed(bed)
}
\author{
  Mark Cowley, 2011-11-29
}

