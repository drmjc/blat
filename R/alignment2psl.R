# Take an alignment of targets to chromosome(s) where the minimum
# requirements are the query name, target name (ie chromosome),
# start, stop and strand of the alignment.
#
# Parameters:
# 	query: a character vector of query names
# 	target: a vector of chromosome numbers; either chr1, chrX or 1, "X" etc...
# 	start/stop: numerical vectors of start/stop of the alignment
# 	strand: "+" or "-" for each alignment
#
# Value:
#  a psl data.frame where it is assumed that the alignment is a perfect one between
#  the query and the target, with only one alignment block.
#
# Mark Cowley, 24 August 2006
#
alignment2psl <- function(query=x[,"query"], target=x[,"chr"], start=x[,"start"], stop=x[,"stop"], strand=x[,"strand"], score=NULL, x) {
	match <- stop - start

    if( any(target == "MT") )
        target[target=="MT"] <- "M"


    if( !all(strand %in% c("-", "+")) ) {
        if( is.numeric(strand) ) {
            strand[strand < 0] <- "-"
            strand[strand >= 0] <- "+"
        }
        else {
            stop("strand should be \"+\" or \"-\"\n")
        }
    }

    #
    # convert chromosome numbers (or X Y M) into
    # chr1 or chrX for eg.
    #
    # if the target is a contig eg NT_039877, then don't
    # prepend the 'chr'.
    #
	if( any(target %in% c(1:30,"X","Y","M")) ) {
        idx <- which( target %in% c(1:30,"X","Y","M") )
        target[idx] <- p("chr", target[idx])
	}

	psl <- new.psl(length(query), !is.null(score), match)
	psl$"Q name" <- query
	psl$"T name" <- target
	psl$"T start" <- start
	psl$"T end" <- stop
	psl$"strand" <- strand
	psl$"tStarts" <- p(start, ",")
	psl$"score" <- score

	return( psl )
}
