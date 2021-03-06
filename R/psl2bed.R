#' Convert a psl alignment from BLAT into a BED alignment.
#'
#' @section TODO:
#' thickStart and thickEnd are currently just the same as "T start" and "T end"
#' 
#' @param psl A psl alignment object
#' @param col the colour, in the form of a character vector. default: \dQuote{0,0,0}
#' @param calc.score logical. if \code{TRUE}, then the BED score (column 5) will be calculated.
#'  see \link{pslScore}. if \code{FALSE}, then the Score (column 5) will be the psl$score, if present
#' or 1000.
#' 
#' @return a BED-formatted \code{data.frame}.
#' 
#' @author Mark Cowley, 5 Jan 2007
#' @export
#' @importFrom mjcbase "colclasses<-" recycle
#' 
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' head(psl2bed(psl))
#' head(psl2bed(psl, calc.score=T))
#' 
#' psl <- import.psl(f, score=TRUE)
#' head(psl2bed(psl))
#' head(psl2bed(psl, calc.score=T))
psl2bed <- function(psl, col="0,0,0", calc.score=FALSE) {
    if( ! "score" %in% colnames(psl) ) {
        if( calc.score)
            psl$score <- pslScore(psl)
        else
            psl$score <- 1000
    }
    maxscore <- max(psl$score)
    psl$score <- round(psl$score / maxscore * 1000, 0)

    thickStart <- as.numeric( sapply(strsplit(psl$tStarts, split=","), "[[", 1) )
    thickEnd <- as.numeric( sapply(strsplit(psl$tStarts, split=","), function(x) x[[length(x)]]) ) + # get the last element of tStarts
                as.numeric( sapply(strsplit(psl$blockSizes, split=","), function(x) x[[length(x)]]) )

    toUCSCcsv <- function(x) {
        paste(paste(x, collapse=","), ",", sep="")
    }
    #
    # uncsv tStarts;
    # subtract "T start" from each numbers;
    # re-csv these corrected numbers.
    #
    tStarts <- lapply( strsplit(psl$tStarts, split=","), as.numeric ) # == tStarts as numeric vectors in a list.
    for(i in 1:length(tStarts))
        tStarts[[i]] <- tStarts[[i]] - psl$"T start"[[i]]
    # ^^^ == tStarts corrected by THE T Start
    tStarts <- sapply(tStarts, toUCSCcsv) # re-csv them!

    bed <- cbind(
        chrom=psl$"T name",
        chromStart=psl$"T start",
        chromEnd=psl$"T end",
        name=psl$"Q name",
        score=psl$score,
        strand=psl$strand,
        thickStart=thickStart,#psl$"T start",
        thickEnd=thickEnd,#psl$"T end",
        itemRgb=recycle(col, nrow(psl))[1:nrow(psl)],
        blockCount=psl$blockcount,
        blockSizes=psl$blockSizes,
        blockStarts=tStarts)#psl$tStarts)

    bed <- as.data.frame(bed)
    colclasses(bed) <- c("character", "numeric", "numeric", # chrom, cS, cE
                        "character", "numeric", "character", # name, score, strand
                        "numeric", "numeric", # thickStart, thickEnd
                        "character", # itemRgb
                        "numeric", "character", "character") # bCount, bSizes, bStarts

    return( bed )
}
