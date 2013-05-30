\name{psl2bed}
\alias{psl2bed}
\title{Convert a psl alignment from BLAT into a BED alignment.}
\usage{
  psl2bed(psl, col = "0,0,0", calc.score = FALSE)
}
\arguments{
  \item{psl}{A psl alignment object}

  \item{col}{the colour, in the form of a character vector.
  default: \dQuote{0,0,0}}

  \item{calc.score}{logical. if \code{TRUE}, then the BED
  score (column 5) will be calculated.  see
  \link{pslScore}. if \code{FALSE}, then the Score (column
  5) will be the psl$score, if present or 1000.}
}
\value{
  a BED-formatted \code{data.frame}.
}
\description{
  Convert a psl alignment from BLAT into a BED alignment.
}
\section{TODO}{
  thickStart and thickEnd are currently just the same as "T
  start" and "T end"
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
head(psl2bed(psl))
head(psl2bed(psl, calc.score=T))

psl <- import.psl(f, score=TRUE)
head(psl2bed(psl))
head(psl2bed(psl, calc.score=T))
}
\author{
  Mark Cowley, 5 Jan 2007
}

