#' determine the order of rows in a psl object, based on grouping
#' query, then score, then target name.
#'
#' @inheritParams sort.psl
#' 
#' @author Mark Cowley, 3 August 2006
#' @export
#' 
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' head(order.psl(psl))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=TRUE) # this should already sort the psl object.
#' head(order.psl(psl))
order.psl <- function(psl) {
	if(! "score" %in% colnames(psl)) psl$score <- pslScore(psl)
	res <- order(psl$"Q name", -psl$score, psl$"T name")
	
	return(res)
}

#' sort a psl file such that each query is together, then by score
#' then by target name
#'
#' @param psl a psl object
#' 
#' @author Mark Cowley, 3 August 2006
#' @export
#' 
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' head(sort.psl(psl))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=TRUE)
#' head(sort.psl(psl))
sort.psl <- function(psl) {
	if(! "score" %in% colnames(psl)) {
		message("adding score to psl data")
		psl$score <- pslScore(psl)
	}
	psl <- psl[order.psl(psl), ]
	rownames(psl) <- 1:nrow(psl)
	return( psl )
}

#' obtain the best hit for each query
#' 
#' Like google's "i'm feeling lucky", it will get the
#' best hit for each query
#' @inheritParams sort.psl
#' @param is.sorted logical. default=\code{FALSE}
#' @author Mark Cowley, 3 August 2006
#' @export
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' head(psl.besthit(psl))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=TRUE) # this should already sort the psl object.
#' head(psl.besthit(psl))
psl.besthit <- function(psl, is.sorted=FALSE) {
	if(! "score" %in% colnames(psl)) {
		message("adding score to psl data")
		psl$score <- pslScore(psl)
	}
	if(!is.sorted)
		psl <- sort.psl(psl)

	psl <- psl[match(unique(psl$"Q name"), psl$"Q name"), ]
	return( psl )
}

#' convert psl file to UCSC-style coordinates
#'
#' @inheritParams sort.psl
#' @author Mark Cowley, 9 August 2006
#' @export
#' @importFrom mjcbase rowapply trim
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' head(psl2coordinates(psl))
psl2coordinates <- function(psl) {
	rowapply(psl, function(x) paste(trim(x[[14]]), ":", trim(x[[16]]), "-", trim(x[[17]]), sep=""))
}
