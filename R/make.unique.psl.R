#' Convert a PSL alignment data.frame into a minimal 1:1 alignment file
#' where each query has at most one target
#'
#' @section Ambiguity Rules:
#'     If there's only one alignment for a query, then no change
#'     If there's > 1 alignment for a query, work out the pslScore's
#'      for each alignment, and choose the best one. If there is a tie
#'      for the best alignment, then the one closest to the top of the
#'      \code{data.frame} is chosen. In the case of a tie, warning messages are
#'      displaying indicating to the user whether the ties were all on
#'      the same chromosome or not.
#'
#' @param psl a psl data.frame; see import.psl
#' @param quick if \code{TRUE}, then only the best hit for each transcript is returned
#'            with no regard for ties. If there is a score tie for a transcript
#'            then the first one (alphabetically) will be returned.
#' @param ambig.thresh undocumented
#'
#' @return a \code{data.frame} with 1 row per query containing the best alignment for
#'     each query.
#'
#' @author Mark Cowley, 12 April 2006
#' @export
#' 
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' make.unique.psl(psl)
#' 
make.unique.psl <- function(psl, quick=FALSE, ambig.thresh=0.95) {
    queries <- unique(psl[,"Q name"])

    if( quick ) {
        #
        # just return the best hit for each transcript, completely ignoring ties.
        #
        if( !"score" %in% colnames(psl) )
            psl$score <- pslScore(psl)

        psl <- psl[order(psl$"Q name", -psl$score), ]
        psl <- psl[match(queries, psl$"Q name"), ]

        return( psl )
    }
    multihits <- unique( psl[which(duplicated(psl[, "Q name"])), "Q name"] )

    #
    # store the hits for which there was only one alignment, and remove
    # them from psl object to make the ensuing searches quicker.
    #
    res <- psl[match(setdiff(queries, multihits), psl[,"Q name"]), ]
    psl <- psl[-match(setdiff(queries, multihits), psl[,"Q name"]), ]
    # find the first hit for each of these multihits

    rm.idx <- NULL
    ambig.queries <- NULL

    multihits.rowindices <- c(match(multihits, psl[, "Q name"]), nrow(psl)+1) ## add element so loop algo is easier
    for( i in 1:length(multihits) ) {
        rows <- c( multihits.rowindices[i]:(multihits.rowindices[i+1]-1) )

        # get scores and reorder.
        if( "score" %in% colnames(psl) )
            scores <- psl$score[rows]
        else
            scores <- pslScore(psl[rows,])
        rows <- rows[order(scores, decreasing=TRUE)]
        scores <- scores[order(scores)]

        #
        # section to handle ties:
        #
#
# How to handle ties and ambiguous alignments?
# 1) as per UCSC (96% identity criteria)
# 2) a score that incorporates pslScore and a measure of ambiguity.
#
# 1)   http://genome.ucsc.edu/cgi-bin/hgTables: all_mrna Track Description
#
#    GenBank mouse mRNAs were aligned against the genome using the blat program.
#    When a single mRNA aligned in multiple places, the alignment having the
#    highest base identity was found. Only alignments having a base identity
#    level within 0.5% of the best and at least 96% base identity with the
#    genomic sequence were kept.
#
#
# 2)
# A possible ambiguity score:
# (1 - nambig/nalignments ) * pslScore / len(Query) * 100
# so if perfect, unambiguous alignment:
#     (1 - 0/1) * score / score * 100 = 100
# if a good, unambiguous alignment:
#    (1 - 0/1) * score / length * 100 = 90+%
# poor, unambiguous alignment:
#    (1 - 0/1) * score / length * 100 =~ 50%
# good, ambiguous alignment:
#    (1 - 5/10) * score / length * 100 =~ 50%
#
        if(scores[1] >= scores[2] * 1/ambig.thresh) {
            #
            # there's a tie of at least 2 hits for the best alignment score
            #
            tiedrows <- rows[ scores==scores[1] ]
            ambig.queries <- c(ambig.queries, rows[ scores[1] ])

##             if( alleq(psl[tiedrows, "T name"]) )
##                 # then they're on the same chromosome
## ##                 warning( paste0(multihits[i], " has ", length(tiedrows), " hits to same chromosome") )
##                 dummy <- 1
##             else
##                 # then they're on the different chromosomes
##                 warning( paste0(multihits[i], " has ", length(tiedrows), " hits to different chromosomes") )



        }
        rm.idx <- c(rm.idx, rows[2:length(rows)])
    }

    res <- rbind(res, psl[-rm.idx,])
    # queries contains the accession no's in original order
    res <- res[match(queries, res[,"Q name"]),]

    return( res )
}

