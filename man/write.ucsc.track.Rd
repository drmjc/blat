\name{write.psl.track}
\alias{write.bed.track}
\alias{write.psl.track}
\alias{write.ucsc.track}
\title{Write a psl object to a UCSC-style track}
\usage{
  write.psl.track(..., header = FALSE)

  write.bed.track(..., header = TRUE)

  write.ucsc.track(psl, file, name = "localBLAT",
    description = "local BLAT alignment", splitchr = FALSE,
    pack.tracks = c("knownGene", "mgcGenes", "mrna", "multiz17way"),
    dense.tracks = c("ruler", "refGene", "xenoRefGene", "ensGene", "intronEst", "snp126", "rmsk", "cytoBand", "stsMapMouseNew", "miRNA"),
    colour = c(0, 0, 0), useScore = FALSE, method = "psl",
    header = TRUE)
}
\arguments{
  \item{\dots}{arguments passed to write.ucsc.track}

  \item{psl}{a PSL object}

  \item{file}{the output file name}

  \item{name}{the track name}

  \item{description}{the track description}

  \item{splitchr}{logical. default=\code{FALSE}}

  \item{pack.tracks}{character vector of UCSC track names
  to be displaed in 'pack' form}

  \item{dense.tracks}{character vector of UCSC track names
  to be displaed in 'dense' form}

  \item{colour}{numeric vector[3]}

  \item{useScore}{logical}

  \item{method}{one of \dQuote{psl}, \dQuote{bed},}

  \item{header}{logical. default=\code{TRUE}}
}
\value{
  nothing. writes files.
}
\description{
  Write a psl object to a UCSC-style track
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
out <- tempfile(fileext=".psl")
write.psl.track(psl, file=out)
cat("wrote to: ", out, "\\n")
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)
out <- tempfile(fileext=".bed")
write.bed.track(psl, file=out)
cat("wrote to: ", out, "\\n")
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=FALSE)

out <- tempfile(fileext=".psl")
write.ucsc.track(psl, file=out, method="psl")
cat("wrote to: ", out, "\\n")

out <- tempfile(fileext=".bed")
write.ucsc.track(psl, file=out, method="bed")
cat("wrote to: ", out, "\\n")
}
\author{
  Mark Cowley, 2013-05-30
}

