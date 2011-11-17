write.psl.track <- function(..., header=F) {
    write.ucsc.track(..., header=header, method="psl")
}

write.bed.track <- function(..., header=T) {
    write.ucsc.track(..., header=header, method="bed")
}

write.ucsc.track <- function( psl, file, name="localBLAT", description="local BLAT alignment", splitchr=F,
                             pack.tracks=c("knownGene", "mgcGenes", "mrna", "multiz17way" ),
                             dense.tracks=c("ruler", "refGene", "xenoRefGene", "ensGene", "intronEst", "snp126", "rmsk", "cytoBand", "stsMapMouseNew", "miRNA"),
                             colour=c(0,0,0),
                             useScore=F,
                             method="psl",
                             header=T ) {

    #
    # reorder the psl to genomic order
    #

    psl <- psl[order(psl$"T name", psl$"T start"), ]

    if( splitchr ) {

        FILE <- file
        chrs <- grep("chr", sort.num( unique(psl$"T name") ), value=T)
        for(i in 1:length(chrs)) {
            file <- p(FILE, ".", chrs[i], ".track" )
            x <- psl[psl$"T name" == chrs[i], ]

            write.psl.track(x, file, paste(name, chrs[i], sep="."), description, splitchr=F, method=method, header=header)
        }

        #
        # are there any alignments to contigs??
        #
        chrs <- setdiff(sort.num( unique(psl$"T name") ), chrs)
        if( length(chrs) > 0 ) {
            file <- p(FILE, ".", "chrUn_random", ".track")
            x <- psl[psl$"T name" %in% chrs, ]

            write.psl.track(x, file, p(name, ".ChrUn_random"), description, splitchr=F,pack.tracks=pack.tracks,
                            dense.tracks=dense.tracks, colour=colour, useScore=useScore, method=method, header=header)
        }
    }
    else {
        sep <- "\t"
        OUT <- file(file, "w")
        chr <- unique(psl$"T name")[1]
        if( header ) {
            write(p("browser position ", psl2coordinates(psl[1,])), OUT)
            write(p("track name=", name, " description=\"", description, "\" color=", paste(colour, collapse=","),
                    " visibility=2 useScore=", as.numeric(useScore)), OUT)

            for(track in pack.tracks)
                write(p("browser pack ", track), OUT)
            for(track in dense.tracks)
                write(p("browser dense ", track), OUT)
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

        write.delim(psl, OUT, sep=sep, col.names=F)

        close(OUT)
    }
}


#
# Write a psl object out to a psl file.
# Useful if you want to use the kent utilities on the
# psl object that you've been modifying in R.
#
# Parameters:
#   psl: a psl data.frame from import.psl or new.psl (with or without score column)
#   out: the outfilename or con (see write, write.delim)
#   header: write the psl header, T/F
#
# Mark Cowley, 13 March 2007
#
write.psl <- function(psl, out, header=T) {
    if( "score" %in% colnames(psl))
        psl[, "score"] <- NULL

    theheader <-
"psLayout version 3

match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T    T       block   blockSizes      qStarts  tStarts
        match   match           count   bases   count   bases           name            size    start   end     name            size    startend     count
---------------------------------------------------------------------------------------------------------------------------------------------------------------"
    if( header )
        write(theheader, out)

    write.delim(psl, out, append=header, col.names=F)
}
