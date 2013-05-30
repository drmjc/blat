\name{make.unique.psl}
\alias{make.unique.psl}
\title{Convert a PSL alignment data.frame into a minimal 1:1 alignment file
where each query has at most one target}
\usage{
  make.unique.psl(psl, quick = FALSE, ambig.thresh = 0.95)
}
\arguments{
  \item{psl}{a psl data.frame; see import.psl}

  \item{quick}{if \code{TRUE}, then only the best hit for
  each transcript is returned with no regard for ties. If
  there is a score tie for a transcript then the first one
  (alphabetically) will be returned.}

  \item{ambig.thresh}{undocumented}
}
\value{
  a \code{data.frame} with 1 row per query containing the
  best alignment for each query.
}
\description{
  Convert a PSL alignment data.frame into a minimal 1:1
  alignment file where each query has at most one target
}
\section{Ambiguity Rules}{
  If there's only one alignment for a query, then no change
  If there's > 1 alignment for a query, work out the
  pslScore's for each alignment, and choose the best one.
  If there is a tie for the best alignment, then the one
  closest to the top of the \code{data.frame} is chosen. In
  the case of a tie, warning messages are displaying
  indicating to the user whether the ties were all on the
  same chromosome or not.
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
make.unique.psl(psl)
}
\author{
  Mark Cowley, 12 April 2006
}

