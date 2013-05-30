\name{psl.besthit}
\alias{psl.besthit}
\title{obtain the best hit for each query}
\usage{
  psl.besthit(psl, is.sorted = FALSE)
}
\arguments{
  \item{is.sorted}{logical. default=\code{FALSE}}

  \item{psl}{a psl object}
}
\description{
  Like google's "i'm feeling lucky", it will get the best
  hit for each query
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
head(psl.besthit(psl))

f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=TRUE) # this should already sort the psl object.
head(psl.besthit(psl))
}
\author{
  Mark Cowley, 3 August 2006
}

