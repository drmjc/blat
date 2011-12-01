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
#' 
import.psl <- function(x, score=T, target2chromosome=F) {
	#
	# if there's a psl header, read it, then skip over it
	# when importing the psl data as a table.
	#
	colnames <- colnames( new.psl(score=F) )
	skip <- 0

	#
	# check for the existence of a header, and skip it if it exists.
	#
	tmp <- readLines(x, 10)
	if( grepT("-{2,}", tmp) ) { # then there is a header
		skip <- grep("-{2,}", tmp)
		for(i in skip:length(tmp)) {
			if( nchar(tmp[i]) == 0 )
				skip <- skip + 1
			else
				break
		}
	}

	psl <- read.delim(x, as.is=T, skip=skip, header=F, comment.char='')
	colnames(psl) <- colnames

	if( score ) {
		cat("calculating psl scores\n")
		psl$score <- pslScore(psl)
		psl <- order.psl(psl)
	}

	if( target2chromosome ) {
		tmp <- p("chr", accession2chr( psl$"T name" ))
		if( any(tmp == "chrNA") )
			tmp[tmp == "chrNA"] <- psl$"T name"[tmp == "chrNA"]
		psl$"T name" <- tmp
	}

	return( psl )
}


#' import a UCSC-format BED file
#'
#' @param x the path to a bed file
#' @return a \code{data.frame} representation of the BED file. Any browser
#' or track lines will be lost upon import.
#' @author Mark Cowley, 2011-11-28
#' @export
import.ucsc.bed <- function(x) {
	skip <- 0
	header <- readLines(x, 2)
	skip <- 0
	if( any(grepl("^track", header)) ) skip <- skip + 1
	if( any(grepl("^browser", header)) ) skip <- skip + 1

	tmp <- read.delim(x, as.is=T, header=F, skip=skip)
	idx <- grep("^track", tmp[,1])
	if( length(idx) > 0 )
		tmp <- tmp[-idx,]

	COLNAMES <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")

	colnames(tmp) <- COLNAMES[1:ncol(tmp)]
	return(tmp)
}
