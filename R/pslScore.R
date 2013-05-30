## Code ported from BLAT source code to determine the pslScore for a psl alignment
##
## http://genome.ucsc.edu/FAQ/FAQblat#blat4
##
## Mark Cowley, 12 April 2006
##

#
# colnames and their indices of a psl result table:
#
# [,1]    [,2]         [,3]         [,4]  [,5]         [,6]         [,7]
# "match" "mis-match" "rep-match" "N's" "Q gapcount" "Q gapbases" "T gapcount"
# [,8]         [,9]     [,10]    [,11]    [,12]     [,13]   [,14]    [,15]
# "T gapbases" "strand" "Q name" "Q size" "Q start" "Q end" "T name" "T size"
# [,16]     [,17]   [,18]        [,19]        [,20]     [,21]
# "T start" "T end" "blockcount" "blockSizes" "qStarts" "tStarts"
#

#' calculate the score from a psl alignment
#' 
#' This calculates the score for a PSL alignment, based on C code from
#' Jim Kent, see comment in pslScore.R.
#' This has been optimised, breaking the problem into smaller chunks.
#' 
#' @param psl a PSL object
#' @param isMrna logical: protein scores are 3x size of mrna/nucleotide scores. 
#' 
#' @return vector psl alignment score
#' 
#' @author Mark Cowley, 12 April 2006
#' @export
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)
#' psl$score <- pslScore(psl, FALSE)
#' head(psl)
#' 
#' # or the simpler appraoch
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=TRUE)
#' head(psl)
#' 
pslScore <- function(psl, isMrna=!pslIsProtein(psl[1,])) {

	BATCH_SIZE <- 1000

	# if( nrow(psl) > 1 )
		# init.progress.meter2(nrow(psl))

	if( nrow(psl) > BATCH_SIZE ) {
		#
		# Running pslScore on all 37837 records in psl for the compugen22k 65mers using
		# a loop from 1:37837 takes 30 mins:
		# date(); scores <- pslScore(psl); date()
		# [1] "Wed Apr 12 16:41:55 2006"
		# [1] "Wed Apr 12 17:16:19 2006"
		#
		# ## using this loop:
		# for(i in 1:nrow(psl)) {
		#     scores[i] <- pslScore(psl[i,])
		# }
		#
		# Splitting the 37837 into 38 loops of 1000 rows, and for each set of 1000,
		# compute pslScore via loops from 1:1000 takes 90 seconds:
		# date(); scores <- pslScore(psl); date()
		# [1] "Wed Apr 12 15:48:43 2006"
		# [1] "Wed Apr 12 15:50:26 2006"
		#
		# Note I tried the loop from 1:37837 both with full and empty RAM and still
		# got approx 30 mins both times?
		#
		# What gives?
		#
		# Note, I've tried replacing bitShiftR below with floor(x/2) which didn't
		# speed things up much...
		#
		scores <- rep(0, nrow(psl))
		for(idx in split_into_batches(1:nrow(psl), batch.size=BATCH_SIZE)) {
			scores[idx] <- .pslScore(psl[idx,], isMrna)
		}
		return( scores )
	}
	else {
		return( .pslScore(psl, isMrna) )
	}
}

#' @importFrom bitops bitShiftR
.pslScore <- function(psl, isMrna) {
	# stop("I dropped bitops dependency, as it doesn't have namespace; thus bitops::bitShiftR no longer works.")
	if( nrow(psl) > 1 ) {
		scores <- rep(0, nrow(psl))

		for(i in 1:nrow(psl)) {
			scores[i] <- pslScore(psl[i,], isMrna)
			# update.progress.meter2()
		}
		return( scores )
	}
	else {
		sizeMul <- ifelse(isMrna, 1, 3)

		return( sizeMul * (psl$"match" + ( bitShiftR(psl$"rep-match", 1) )) -
				sizeMul * psl$"mis-match" - psl$"Q gapcount" - psl$"T gapcount" )
	}
}
## C-code from Jim Kent.
## int pslScore(const struct psl *psl) {
##     /* Return score for psl. */
##     int sizeMul = pslIsProtein(psl) ? 3 : 1;
##
##     # x >> 1 right shifts x by 1 bit (ie 0111 -> 0011; 7 -> 3)
##     return sizeMul * (psl->match + ( psl->repMatch>>1)) -
##              sizeMul * psl->misMatch - psl->qNumInsert - psl->tNumInsert;
## }


#' calculate the % identity of a psl alignment
#' 
#' Many years later, this version looks more robust than pslPercentIdentity2 and pslPercentIdentity3
#' 
#' @inheritParams pslScore
#' @export numeric vector
#' @author Mark Cowley, 12 April 2006
#' @export
#' 
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=TRUE)[1:5,]
#' pslPercentIdentity(psl, FALSE)
pslPercentIdentity <- function(psl, isMrna=!pslIsProtein(psl)) {
	BATCH_SIZE <- 1000

	if( nrow(psl) > BATCH_SIZE ) {
		#
		# see comment in pslScore re why this takes so damn long looping
		# from 1:nrow(psl) when if you split into chunks of 1000 rows it
		# runs many orders faster!?
		#
        scores <- rep(0, nrow(psl))

# 		# init.progress.meter2( nrow(psl) )
		for(idx in split_into_batches(1:nrow(psl), batch.size=BATCH_SIZE)) {
			scores[idx] <- pslPercentIdentity(psl[idx,], isMrna)
        }

        return( scores )
	}
	else if( nrow(psl) <= BATCH_SIZE && nrow(psl) > 1 ) {
		scores <- rep(0, nrow(psl))
		for(i in 1:nrow(psl)) {
			scores[i] <- pslPercentIdentity(psl[i,], isMrna)
# 			# update.progress.meter2()
		}
		return( scores )
	}
	else {
    	return( round(100.0 - pslCalcMilliBad(psl, isMrna) * 0.1, 1) )
	}
}

#' calculate the % identity of a psl alignment V2
#' 
#' @inheritParams pslScore
#' @author Mark Cowley, 12 April 2006
#' @importFrom mjcbase split_into_batches
pslPercentIdentity2 <- function(psl, isMrna=!pslIsProtein(psl)) {
	BATCH_SIZE <- 1000

	if( nrow(psl) > BATCH_SIZE ) {
		#
		# see comment in pslScore re why this takes so damn long looping
		# from 1:nrow(psl) when if you split into chunks of 1000 rows it
		# runs many orders faster!?
		#
        scores <- rep(0, nrow(psl))

# 		# init.progress.meter2( nrow(psl) )
		for(idx in split_into_batches(1:nrow(psl), batch.size=BATCH_SIZE)) {
			for(i in 1:length(idx)) {
				scores[idx[i]] <- round(100.0 - pslCalcMilliBad(psl[idx[i],], isMrna) * 0.1, 1)
# 				# update.progress.meter2()
			}
        }
        return( scores )
	}
}

#' calculate the % identity of a psl alignment V3
#' 
#' should be slower than V2
#' 
#' @inheritParams pslScore
#' @author Mark Cowley, 12 April 2006
#' @importFrom mjcbase split_into_batches
pslPercentIdentity3 <- function(psl, isMrna=!pslIsProtein(psl)) {
	scores <- rep(0, nrow(psl))

# 	# init.progress.meter2( nrow(psl) )
	for(i in 1:nrow(psl)) {
		scores[i] <- round(100.0 - pslCalcMilliBad(psl[i,], isMrna) * 0.1, 1)
# 		# update.progress.meter2()
	}
	return( scores )
}




#' internal function used to determing a psl alignment score.
#' 
#' @inheritParams pslScore
#' @return integer. largely undocumented
#' @author Mark Cowley, 12 April 2006
#' 
#' @examples
#' ## Note the order within psl differs if import.psl(score=TRUE)
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)[1:5,]
#' pslCalcMilliBad(psl[3,])
#' psl[3,]
#' 
#' psl <- import.psl(f, score=TRUE)[1:5,]
#' pslCalcMilliBad(psl[3,])
#' psl[3,]
#' 
pslCalcMilliBad <- function(psl, isMrna=!pslIsProtein(psl)) {

	if( nrow(psl) > 1 ) {
		stop( "psl should have only 1 row\n" )
	}
	else {
		sizeMul <- ifelse( isMrna, 1, 3 )
		milliBad <- 0

		qAliSize <- sizeMul * (psl$"Q end" - psl$"Q start") # Qend - Qstart
		tAliSize <- psl$"T end" - psl$"T start" # Tend - Tstart
		aliSize <- min(qAliSize, tAliSize)

		if (aliSize <= 0)
			return( 0 )

		sizeDif <- qAliSize - tAliSize
		if (sizeDif < 0) {
			if (isMrna)
				sizeDif <- 0
			else
				sizeDif <- -sizeDif
		}

		insertFactor <- psl$"Q gapcount" # QgapNum
		if (!isMrna)
			insertFactor <- insertFactor + psl$"T gapcount" # TgapNum

		total <- (sizeMul * (psl$"match" + psl$"rep-match" + psl$"mis-match")) # match + rep.match + mis-match

		if (total != 0)
			milliBad <- (1000 * (psl$"mis-match"*sizeMul + insertFactor + round(3*log(1+sizeDif)))) / total # psl$"mis-match" is mis-match

		return( milliBad )
	}
}
## int pslCalcMilliBad(struct psl *psl, boolean isMrna) {
##     /* Calculate badness in parts per thousand. */
##     int sizeMul = pslIsProtein(psl) ? 3 : 1;
##     int qAliSize, tAliSize, aliSize;
##     int milliBad = 0;
##     int sizeDif;
##     int insertFactor;
##     int total;
##
##     qAliSize = sizeMul * (psl->qEnd - psl->qStart);
##     tAliSize = psl->tEnd - psl->tStart;
##     aliSize = min(qAliSize, tAliSize);
##     if (aliSize <= 0)
##         return 0;
##     sizeDif = qAliSize - tAliSize;
##     if (sizeDif < 0)
##         {
##         if (isMrna)
##             sizeDif = 0;
##         else
##             sizeDif = -sizeDif;
##         }
##     insertFactor = psl->qNumInsert;
##     if (!isMrna)
##         insertFactor += psl->tNumInsert;
##
##     total = (sizeMul * (psl->match + psl->repMatch + psl->misMatch));
##     if (total != 0)
##         milliBad = (1000 * (psl->misMatch*sizeMul + insertFactor +
##     	round(3*log(1+sizeDif)))) / total;
##     return milliBad;
## }


#' is the psl alignment from a protein?
#'
#' fixed the bug on the - strand where psl$"T end" - (_ + 3*_)
#'
#' @inheritParams pslScore
#' 
#' @return logical
#' 
#' @author Mark Cowley, 12 April 2006
#' @export
#' @importFrom mjcbase uncsv
#' 
#' @examples
#' f <- file.path(system.file(package="blat"), "examples", "test.psl")
#' psl <- import.psl(f, score=FALSE)[1:5,]
#' pslIsProtein(psl)
#' 
#' psl <- import.psl(f, score=TRUE)[1:5,]
#' pslIsProtein(psl)
#' 
pslIsProtein <- function(psl) {
	if( nrow(psl) > 1 ) {
		scores <- rep(F, nrow(psl))
		for(i in 1:nrow(psl))
			scores[i] <- pslIsProtein(psl[i,])
		return( scores )
	}
	else {
		## is psl a protein psl (are it's blockSizes and scores in protein space)
		lastBlock <- psl$"blockcount"
		tstarts <- as.numeric( uncsv(psl$"tStarts") )
		blockSizes <- as.numeric( uncsv(psl$"blockSizes") )

		return( ((psl$"strand" == '+' ) && (psl$"T end" == tstarts[lastBlock] + 3 * blockSizes[lastBlock]))
				||
				((psl$"strand" == '-') && (psl$"T start" == (psl$"T end"-(tstarts[lastBlock] + 3*blockSizes[lastBlock]))))
			)
	}
}
## boolean pslIsProtein(const struct psl *psl) {
##     /* is psl a protein psl (are it's blockSizes and scores in protein space)
##     */
##     int lastBlock = psl->blockCount - 1;
##
##     return ( ((psl->strand[1] == '+' ) && (psl->tEnd == psl->tStarts[lastBlock] + 3*psl->blockSizes[lastBlock]))
##               ||
##              ((psl->strand[1] == '-') && (psl->tStart == (psl->tSize-(psl->tStarts[lastBlock] + 3*psl->blockSizes[lastBlock]))))
##            );
## }



