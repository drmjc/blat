\name{pslPercentIdentity}
\alias{pslPercentIdentity}
\title{calculate the \% identity of a psl alignment}
\usage{
  pslPercentIdentity(psl, isMrna = !pslIsProtein(psl))
}
\arguments{
  \item{psl}{a PSL object}

  \item{isMrna}{logical: protein scores are 3x size of
  mrna/nucleotide scores.}
}
\description{
  Many years later, this version looks more robust than
  pslPercentIdentity2 and pslPercentIdentity3
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=TRUE)[1:5,]
pslPercentIdentity(psl, FALSE)
}
\author{
  Mark Cowley, 12 April 2006
}

