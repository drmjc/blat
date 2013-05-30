\name{order.psl}
\alias{order.psl}
\title{determine the order of rows in a psl object, based on grouping
query, then score, then target name.}
\usage{
  order.psl(psl)
}
\arguments{
  \item{psl}{a psl object}
}
\description{
  determine the order of rows in a psl object, based on
  grouping query, then score, then target name.
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
head(order.psl(psl))

f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=TRUE) # this should already sort the psl object.
head(order.psl(psl))
}
\author{
  Mark Cowley, 3 August 2006
}

