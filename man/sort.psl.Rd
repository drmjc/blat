\name{sort.psl}
\alias{sort.psl}
\title{sort a psl file such that each query is together, then by score
then by target name}
\usage{
  sort.psl(psl)
}
\arguments{
  \item{psl}{a psl object}
}
\description{
  sort a psl file such that each query is together, then by
  score then by target name
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
head(sort.psl(psl))

f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=TRUE)
head(sort.psl(psl))
}
\author{
  Mark Cowley, 3 August 2006
}

