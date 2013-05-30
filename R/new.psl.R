#' Make a new empty psl table.
#' 
#' @param nrow how many rows would you like? If > 0 then the appropriate
#'   default values are assigned to each column
#' @param score logical: if \code{TRUE}, add a column to the end called "score"
#' @param alignment.length A vector of alignment lengths which inform the
#'   function as to which data columns to change. This is recycled if
#'   necessary.
#' @return a data.frame of the same specs as import.psl()
#' @author Mark Cowley, 24 August 2006
#' @export
#' @seealso \code{\link{import.psl}}
#' @importFrom mjcbase "colclasses<-" recycle
new.psl <- function(nrow=0, score=TRUE, alignment.length=10) {
	res <- as.data.frame(matrix(NA, nrow=nrow, ncol=22))
	colclasses(res) <- c("integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "character", "character", "integer", "integer", "integer", "character", "integer", "integer", "integer", "integer", "character", "character", "character", "numeric")
	colnames(res) <- c("match", "mis-match", "rep-match", "N's", "Q gapcount", "Q gapbases", "T gapcount", "T gapbases", "strand", "Q name", "Q size", "Q start", "Q end", "T name", "T size", "T start", "T end", "blockcount", "blockSizes", "qStarts", "tStarts", "score")

	if( nrow > 0 ) {
		alignment.length <- recycle(alignment.length, nrow)

		res$match <- alignment.length
		for(col in c("mis-match", "rep-match", "N's", "Q gapcount", "Q gapbases", "T gapcount", "T gapbases"))
			res[,col] <- 0
		res$"strand" <- "+"
		res$"Q size" <- alignment.length
		res$"Q start" <- 0
		res$"Q end" <- alignment.length
		res$"T size" <- 1234567890
		res$"T start" <- 0
		res$"T end" <- alignment.length
		res$"blockcount" <- 1
		res$"blockSizes" <- paste(sep="", alignment.length, ",")
		res$"qStarts" <- "0,"
		res$"tStarts" <- "0,"
		res$"score" <- alignment.length
	}

	if( !score ) {
		res <- res[, 1:21]
	}

	return( res )
}
