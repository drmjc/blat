\name{alignment2psl}
\alias{alignment2psl}
\title{Take an alignment of targets to chromosome(s) where the minimum
requirements are the query name, target name (ie chromosome),
start, stop and strand of the alignment.}
\usage{
  alignment2psl(query = x[, "query"], target = x[, "chr"],
    start = x[, "start"], stop = x[, "stop"],
    strand = x[, "strand"], score = NULL, x)
}
\arguments{
  \item{query}{a character vector of query names}

  \item{target}{a vector of chromosome numbers; either
  chr1, chrX or 1, "X" etc...}

  \item{start}{numerical vectors of start of the alignment}

  \item{stop}{numerical vectors of stop of the alignment}

  \item{strand}{"+" or "-" for each alignment}

  \item{score}{alignment score. default=\code{NULL}}

  \item{x}{the alignment object}
}
\value{
  a psl data.frame where it is assumed that the alignment
  is a perfect one between the query and the target, with
  only one alignment block.
}
\description{
  Take an alignment of targets to chromosome(s) where the
  minimum requirements are the query name, target name (ie
  chromosome), start, stop and strand of the alignment.
}
\author{
  Mark Cowley, 24 August 2006
}

