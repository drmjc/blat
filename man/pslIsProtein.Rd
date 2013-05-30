\name{pslIsProtein}
\alias{pslIsProtein}
\title{is the psl alignment from a protein?}
\usage{
  pslIsProtein(psl)
}
\arguments{
  \item{psl}{a PSL object}
}
\value{
  logical
}
\description{
  fixed the bug on the - strand where psl$"T end" - (_ +
  3*_)
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)[1:5,]
pslIsProtein(psl)

psl <- import.psl(f, score=TRUE)[1:5,]
pslIsProtein(psl)
}
\author{
  Mark Cowley, 12 April 2006
}

