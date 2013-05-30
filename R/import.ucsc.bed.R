#' import a UCSC-format BED file
#'
#' @param x the path to a bed file
#' @return a \code{data.frame} representation of the BED file. Any browser
#' or track lines will be lost upon import.
#' @author Mark Cowley, 2011-11-28
#' @export
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "bed9_with_hdr.bed")
#' head(import.ucsc.bed(f))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "bed9_no_hdr.bed")
#' head(import.ucsc.bed(f))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "bed6_with_hdr.bed")
#' head(import.ucsc.bed(f))
#' f <- file.path(system.file(package="blat"), "examples", "bed6_no_hdr.bed")
#' head(import.ucsc.bed(f))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "bed4_no_hdr.bed")
#' head(import.ucsc.bed(f))
#' 
#' f <- file.path(system.file(package="blat"), "examples", "bed3_no_hdr.bed")
#' head(import.ucsc.bed(f))
#'
import.ucsc.bed <- function(x) {
	skip <- 0
	header <- readLines(x, 10, warn=FALSE)
	skip <- 0
	if( (s <- sum(grepl("^track", header))) > 0 ) skip <- skip + s
	if( (s <- sum(grepl("^browser", header))) > 0 ) skip <- skip + s

	tmp <- read.table(x, as.is=TRUE, header=FALSE, skip=skip, sep="")
	# idx <- grep("^track", tmp[,1])
	# if( length(idx) > 0 )
	# 	tmp <- tmp[-idx,]

	COLNAMES <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")

	colnames(tmp) <- COLNAMES[1:ncol(tmp)]
	return(tmp)
}
# CHANGELOG
# 2013-05-30: improved robustness, added examples & tested.
