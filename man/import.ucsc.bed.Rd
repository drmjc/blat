\name{import.ucsc.bed}
\alias{import.ucsc.bed}
\title{import a UCSC-format BED file}
\usage{
  import.ucsc.bed(x)
}
\arguments{
  \item{x}{the path to a bed file}
}
\value{
  a \code{data.frame} representation of the BED file. Any
  browser or track lines will be lost upon import.
}
\description{
  import a UCSC-format BED file
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "bed9_with_hdr.bed")
head(import.ucsc.bed(f))

f <- file.path(system.file(package="blat"), "examples", "bed9_no_hdr.bed")
head(import.ucsc.bed(f))

f <- file.path(system.file(package="blat"), "examples", "bed6_with_hdr.bed")
head(import.ucsc.bed(f))
f <- file.path(system.file(package="blat"), "examples", "bed6_no_hdr.bed")
head(import.ucsc.bed(f))

f <- file.path(system.file(package="blat"), "examples", "bed4_no_hdr.bed")
head(import.ucsc.bed(f))

f <- file.path(system.file(package="blat"), "examples", "bed3_no_hdr.bed")
head(import.ucsc.bed(f))
}
\author{
  Mark Cowley, 2011-11-28
}

