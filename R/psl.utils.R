# sort a psl file such that each query is together, then by score
# then by target name
#
# Mark Cowley, 3 August 2006
#
order.psl <- function(psl) {
	psl <- psl[order(psl$"Q name", -psl$score, psl$"T name"),]
	rownames(psl) <- 1:nrow(psl)
	return( psl )
}

# Like google's "i'm feeling lucky", it will get the
# best hit for each query
#
# Mark Cowley, 3 August 2006
psl.besthit <- function(psl, is.sorted=F) {
	if(!is.sorted)
		psl <- order.psl(psl)

	psl <- psl[match(unique(psl$"Q name"), psl$"Q name"),]
	return( psl )
}

# generate the coordinates that can be searched at UCSC
#
# Mark Cowley, 9 August 2006
#
psl2coordinates <- function(psl) {
	rowapply(psl, function(x) paste(trim(x[[14]]), ":", trim(x[[16]]), "-", trim(x[[17]]), sep=""))
}
