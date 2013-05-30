#' Import a psl formatted alignment file of queries to targets
#' 
#' @param x the filename of a psl alignment, with or without a pslHeader.
#' @param score logical: if \code{TRUE}, calculate the alignment score
#'  for all genomic alignments
#' @param target2chromosome If psl contains alignments to a chromosome, then
#'   		  you can convert the FastA header chromosome name to the human
#'   		  readable form of \dQuote{chr11} for example.
#' @return a \code{data.frame} representation of a psl file
#' @author Mark Cowley, 12 April 2006
#' @export
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' head(psl)
#' 
#' psl <- import.psl(f, score=TRUE)
#' head(psl)
#' 
import.psl <- function(x, score=TRUE, target2chromosome=FALSE) {
	#
	# if there's a psl header, read it, then skip over it
	# when importing the psl data as a table.
	#
	colnames <- colnames( new.psl(score=FALSE) )
	skip <- 0

	#
	# check for the existence of a header, and skip it if it exists.
	#
	tmp <- readLines(x, 10)
	if( any(grepl("-{2,}", tmp)) ) { # then there is a header
		skip <- grep("-{2,}", tmp)
		for(i in skip:length(tmp)) {
			if( nchar(tmp[i]) == 0 )
				skip <- skip + 1
			else
				break
		}
	}

	psl <- read.delim(x, as.is=TRUE, skip=skip, header=FALSE, comment.char='')
	colnames(psl) <- colnames

	if( score ) {
		cat("calculating psl scores\n")
		psl$score <- pslScore(psl)
		psl <- sort.psl(psl)
	}

	if( target2chromosome ) {
		tmp <- paste0("chr", accession2chr( psl$"T name" ))
		if( any(tmp == "chrNA") )
			tmp[tmp == "chrNA"] <- psl$"T name"[tmp == "chrNA"]
		psl$"T name" <- tmp
	}

	return( psl )
}
