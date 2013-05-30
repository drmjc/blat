\name{new.psl}
\alias{new.psl}
\title{Make a new empty psl table.}
\usage{
  new.psl(nrow = 0, score = TRUE, alignment.length = 10)
}
\arguments{
  \item{nrow}{how many rows would you like? If > 0 then the
  appropriate default values are assigned to each column}

  \item{score}{logical: if \code{TRUE}, add a column to the
  end called "score"}

  \item{alignment.length}{A vector of alignment lengths
  which inform the function as to which data columns to
  change. This is recycled if necessary.}
}
\value{
  a data.frame of the same specs as import.psl()
}
\description{
  Make a new empty psl table.
}
\author{
  Mark Cowley, 24 August 2006
}
\seealso{
  \code{\link{import.psl}}
}

