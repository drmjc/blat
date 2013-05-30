#' grow the intervals in a bed file
#'
#' The galaxy \code{flank} tool returns only the flanks, not the original
#' regions. \code{grow} aims to provide the functionality of \code{flank},\
#' but also return the original regions. ie, the seed regions will grow.
#' Since the regions grow, an \code{offset} is not supported.
#' 
#' @section TODO:
#' Add a check to make sure the regions don't extend beyond the end of a
#' chromosome.
#' 
#' @param x a \code{data.frame} from import.ucsc.blat
#' @param upstream numeric(1) indicating the nuber of bases to grow
#'  in the upstream direction. default=0
#' @param downstream numeric(1) indicating the nuber of bases to grow
#'  in the downstream direction. default=0
#' @return a \code{data.frame}
#' @author Mark Cowley, 2011-11-29
#' @export
grow.bed <- function(x, upstream=0, downstream=0) {
	if("strand" %in% colnames(x)) strand <- x$strand
	else strand <- rep("+", nrow(x))
	
	x$chromStart <- as.integer(x$chromStart - ifelse(strand =="+", rep(upstream,nrow(x)), rep(downstream,nrow(x))))
	x$chromEnd <- as.integer(x$chromEnd + ifelse(strand =="+", rep(downstream,nrow(x)), rep(upstream,nrow(x))))
	x$chromStart <- as.integer(pmax(0, x$chromStart))
	# @TODO make sure that the chromosomes haven't become too long!
	
	x
}
