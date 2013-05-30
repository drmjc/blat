\name{import.psl}
\alias{import.psl}
\title{Import a psl formatted alignment file of queries to targets}
\usage{
  import.psl(x, score = TRUE, target2chromosome = FALSE)
}
\arguments{
  \item{x}{the filename of a psl alignment, with or without
  a pslHeader.}

  \item{score}{logical: if \code{TRUE}, calculate the
  alignment score for all genomic alignments}

  \item{target2chromosome}{If psl contains alignments to a
  chromosome, then you can convert the FastA header
  chromosome name to the human readable form of
  \dQuote{chr11} for example.}
}
\value{
  a \code{data.frame} representation of a psl file
}
\description{
  Import a psl formatted alignment file of queries to
  targets
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
head(psl)

psl <- import.psl(f, score=TRUE)
head(psl)
}
\author{
  Mark Cowley, 12 April 2006
}

