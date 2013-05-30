#' Determine the midpoint of an interval
#'
#' @param x a \code{data.frame} representation of a BED file
#' @return a \code{data.frame} like \code{x}, with chromStart and chromEnd
#'  updated to reflect the midpoint.
#' 
#' @author Mark Cowley, 2011-11-29
#' @export
#' 
#' @examples
#' bed <- data.frame(chrom="chr1", chromStart=c(100,200,300), chromEnd=c(110, 250, 500))
#' midpoint.bed(bed)
midpoint.bed <- function(x) {
	a <- x$chromStart + (x$chromEnd-x$chromStart)/2
	a <- floor(a)
	x$chromStart <- as.integer(a)
	x$chromEnd <- as.integer(a)
	
	x
}
