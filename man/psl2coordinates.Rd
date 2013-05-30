\name{psl2coordinates}
\alias{psl2coordinates}
\title{convert psl file to UCSC-style coordinates}
\usage{
  psl2coordinates(psl)
}
\arguments{
  \item{psl}{a psl object}
}
\description{
  convert psl file to UCSC-style coordinates
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
head(psl2coordinates(psl))
}
\author{
  Mark Cowley, 9 August 2006
}

