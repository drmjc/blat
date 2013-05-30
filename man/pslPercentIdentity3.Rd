\name{pslPercentIdentity3}
\alias{pslPercentIdentity3}
\title{calculate the \% identity of a psl alignment V3}
\usage{
  pslPercentIdentity3(psl, isMrna = !pslIsProtein(psl))
}
\arguments{
  \item{psl}{a PSL object}

  \item{isMrna}{logical: protein scores are 3x size of
  mrna/nucleotide scores.}
}
\description{
  should be slower than V2
}
\author{
  Mark Cowley, 12 April 2006
}

