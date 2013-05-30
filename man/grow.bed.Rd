\name{grow.bed}
\alias{grow.bed}
\title{grow the intervals in a bed file}
\usage{
  grow.bed(x, upstream = 0, downstream = 0)
}
\arguments{
  \item{x}{a \code{data.frame} from import.ucsc.blat}

  \item{upstream}{numeric(1) indicating the nuber of bases
  to grow in the upstream direction. default=0}

  \item{downstream}{numeric(1) indicating the nuber of
  bases to grow in the downstream direction. default=0}
}
\value{
  a \code{data.frame}
}
\description{
  The galaxy \code{flank} tool returns only the flanks, not
  the original regions. \code{grow} aims to provide the
  functionality of \code{flank},\ but also return the
  original regions. ie, the seed regions will grow. Since
  the regions grow, an \code{offset} is not supported.
}
\section{TODO}{
  Add a check to make sure the regions don't extend beyond
  the end of a chromosome.
}
\author{
  Mark Cowley, 2011-11-29
}

