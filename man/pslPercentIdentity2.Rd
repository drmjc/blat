\name{pslPercentIdentity2}
\alias{pslPercentIdentity2}
\title{calculate the \% identity of a psl alignment V2}
\usage{
  pslPercentIdentity2(psl, isMrna = !pslIsProtein(psl))
}
\arguments{
  \item{psl}{a PSL object}

  \item{isMrna}{logical: protein scores are 3x size of
  mrna/nucleotide scores.}
}
\description{
  calculate the \% identity of a psl alignment V2
}
\author{
  Mark Cowley, 12 April 2006
}

