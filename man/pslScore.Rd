\name{pslScore}
\alias{pslScore}
\title{calculate the score from a psl alignment}
\usage{
  pslScore(psl, isMrna = !pslIsProtein(psl[1, ]))
}
\arguments{
  \item{psl}{a PSL object}

  \item{isMrna}{logical: protein scores are 3x size of
  mrna/nucleotide scores.}
}
\value{
  vector psl alignment score
}
\description{
  This calculates the score for a PSL alignment, based on C
  code from Jim Kent, see comment in pslScore.R. This has
  been optimised, breaking the problem into smaller chunks.
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
psl$score <- pslScore(psl, FALSE)
head(psl)

# or the simpler appraoch
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=TRUE)
head(psl)
}
\author{
  Mark Cowley, 12 April 2006
}

