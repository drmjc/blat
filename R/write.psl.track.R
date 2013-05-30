#' @param \dots arguments passed to write.ucsc.track
#' @export
#' @rdname write.ucsc.track
#' @aliases write.psl.track
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' out <- tempfile(fileext=".psl")
#' write.psl.track(psl, file=out)
#' cat("wrote to: ", out, "\n")
write.psl.track <- function(..., header=FALSE) {
    write.ucsc.track(..., header=header, method="psl")
}

#' @export
#' @rdname write.ucsc.track
#' @aliases write.bed.track
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' out <- tempfile(fileext=".bed")
#' write.bed.track(psl, file=out)
#' cat("wrote to: ", out, "\n")
write.bed.track <- function(..., header=TRUE) {
    write.ucsc.track(..., header=header, method="bed")
}

#' Write a psl object to a UCSC-style track
#'
#' @param psl a PSL object
#' @param file the output file name
#' @param name the track name
#' @param description the track description
#' @param splitchr logical. default=\code{FALSE}
#' @param pack.tracks character vector of UCSC track names to be displaed in 'pack' form
#' @param dense.tracks character vector of UCSC track names to be displaed in 'dense' form
#' @param colour numeric vector[3]
#' @param useScore logical
#' @param method one of \dQuote{psl}, \dQuote{bed}, 
#' @param header logical. default=\code{TRUE}
#' 
#' @return nothing. writes files.
#' 
#' @author Mark Cowley, 2013-05-30
#' @export
#' @importFrom mjcbase write.delim sortnum
#' @rdname write.ucsc.track
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' 
#' out <- tempfile(fileext=".psl")
#' write.ucsc.track(psl, file=out, method="psl")
#' cat("wrote to: ", out, "\n")
#' 
#' out <- tempfile(fileext=".bed")
#' write.ucsc.track(psl, file=out, method="bed")
#' cat("wrote to: ", out, "\n")
#' 
write.ucsc.track <- function( psl, file, name="localBLAT", description="local BLAT alignment", splitchr=FALSE,
                             pack.tracks=c("knownGene", "mgcGenes", "mrna", "multiz17way" ),
                             dense.tracks=c("ruler", "refGene", "xenoRefGene", "ensGene", "intronEst", "snp126", "rmsk", "cytoBand", "stsMapMouseNew", "miRNA"),
                             colour=c(0,0,0),
                             useScore=FALSE,
                             method="psl",
                             header=TRUE ) {

    #
    # reorder the psl to genomic order
    #

    psl <- psl[order(psl$"T name", psl$"T start"), ]

    if( splitchr ) {

        FILE <- file
        chrs <- grep("chr", sortnum( unique(psl$"T name") ), value=TRUE)
        for(i in 1:length(chrs)) {
            file <- paste0(FILE, ".", chrs[i], ".track" )
            x <- psl[psl$"T name" == chrs[i], ]

            write.psl.track(x, file, paste(name, chrs[i], sep="."), description, splitchr=FALSE, method=method, header=header)
        }

        #
        # are there any alignments to contigs??
        #
        chrs <- setdiff(sortnum( unique(psl$"T name") ), chrs)
        if( length(chrs) > 0 ) {
            file <- paste0(FILE, ".", "chrUn_random", ".track")
            x <- psl[psl$"T name" %in% chrs, ]

            write.psl.track(x, file, paste0(name, ".ChrUn_random"), description, splitchr=FALSE,pack.tracks=pack.tracks,
                            dense.tracks=dense.tracks, colour=colour, useScore=useScore, method=method, header=header)
        }
    }
    else {
        OUT <- file(file, "w")
        chr <- unique(psl$"T name")[1]
        if( header ) {
            write(paste0("browser position ", psl2coordinates(psl[1,])), OUT)
            write(paste0("track name=", name, " description=\"", description, "\" color=", paste(colour, collapse=","),
                    " visibility=2 useScore=", as.numeric(useScore)), OUT)

            for(track in pack.tracks)
                write(paste0("browser pack ", track), OUT)
            for(track in dense.tracks)
                write(paste0("browser dense ", track), OUT)
##         #
##         # the UCSC genome tracks have +- for '+ strand' and -+ for '- strand'
##         #
##         pos.strand <- grep("\\+", psl$strand)
##         psl$strand[pos.strand] <- "+-"
##         psl$strand[-pos.strand] <- "-+"
        }

        psl$score <- NULL

        if( method %in% c("bed", "BED", "Bed") )
            psl <- psl2bed( psl )

        write.delim(psl, OUT, col.names=FALSE)

        close(OUT)
    }
}


#' Write a psl object out to a psl file.
#' 
#' Useful if you want to use the kent utilities on the
#' psl object that you've been modifying in R. Note this ignores the psl$score, as
#' per the PSL file format.
#' 
#' @param psl a psl \code{data.frame} from import.psl or new.psl (with or without score column)
#' @param out the outfilename or con (see write, write.delim)
#' @param header logical: write the psl header?
#' 
#' @author Mark Cowley, 13 March 2007
#' @export
#' @importFrom mjcbase write.delim
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=TRUE)
#' 
#' out <- tempfile(fileext=".psl")
#' write.psl(psl, out, header=FALSE)
#' cat("wrote to: ", out, "\n")
#' 
#' out <- tempfile(fileext=".psl")
#' write.psl(psl, out, header=TRUE)
#' cat("wrote to: ", out, "\n")
#' 
write.psl <- function(psl, out, header=TRUE) {
    if( "score" %in% colnames(psl))
        psl[, "score"] <- NULL

    theheader <-
"psLayout version 3

match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T    T       block   blockSizes      qStarts  tStarts
        match   match           count   bases   count   bases           name            size    start   end     name            size    startend     count
---------------------------------------------------------------------------------------------------------------------------------------------------------------"
    if( header )
        write(theheader, out)

    write.delim(psl, out, append=header, col.names=FALSE)
}
