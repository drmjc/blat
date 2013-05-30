\name{pslCalcMilliBad}
\alias{pslCalcMilliBad}
\title{internal function used to determing a psl alignment score.}
\usage{
  pslCalcMilliBad(psl, isMrna = !pslIsProtein(psl))
}
\arguments{
  \item{psl}{a PSL object}

  \item{isMrna}{logical: protein scores are 3x size of
  mrna/nucleotide scores.}
}
\value{
  integer. largely undocumented
}
\description{
  internal function used to determing a psl alignment
  score.
}
\examples{
## Note the order within psl differs if import.psl(score=TRUE)
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)[1:5,]
pslCalcMilliBad(psl[3,])
psl[3,]

psl <- import.psl(f, score=TRUE)[1:5,]
pslCalcMilliBad(psl[3,])
psl[3,]
}
\author{
  Mark Cowley, 12 April 2006
}

