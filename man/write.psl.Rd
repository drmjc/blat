\name{write.psl}
\alias{write.psl}
\title{Write a psl object out to a psl file.}
\usage{
  write.psl(psl, out, header = TRUE)
}
\arguments{
  \item{psl}{a psl \code{data.frame} from import.psl or
  new.psl (with or without score column)}

  \item{out}{the outfilename or con (see write,
  write.delim)}

  \item{header}{logical: write the psl header?}
}
\description{
  Useful if you want to use the kent utilities on the psl
  object that you've been modifying in R. Note this ignores
  the psl$score, as per the PSL file format.
}
\examples{
f <- file.path(system.file(package="blat"), "examples", "test.psl")
psl <- import.psl(f, score=TRUE)

out <- tempfile(fileext=".psl")
write.psl(psl, out, header=FALSE)
cat("wrote to: ", out, "\\n")

out <- tempfile(fileext=".psl")
write.psl(psl, out, header=TRUE)
cat("wrote to: ", out, "\\n")
}
\author{
  Mark Cowley, 13 March 2007
}

