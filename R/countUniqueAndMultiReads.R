# countUniqueAndMultiReads.R


`countUniqueAndMultiReads` <- function( filein, readBufferSize=1000000, verbose=TRUE) {

	if ( ! file.exists( filein)) {
		cat( "\nBAM file not found: ", filein)
		return( list( "Unique"=0, "Multi"=0))
	}

	con <- bamReader( filein)
	
	# read in the alignment file one buffer at a time
	hasMore <- TRUE
	nUniq <- nMult <- 0
	repeat {
		if ( ! hasMore) break
		if (verbose) cat( ".")
		chunk <- getNextChunk( con, n=readBufferSize, alignedOnly=TRUE)
		nNow <- size(chunk)
		if ( nNow < 1) break
		if ( nNow < readBufferSize) hasMore <- FALSE

		whoUniq <- which.unique.align(chunk)
		whoMult <- which.multi.read(chunk)

		nUniq <- nUniq + length( whoUniq)
		nMult <- nMult + length( whoMult)
	} # end of each buffer...

	bamClose( con)
	return( list( "Unique"=nUniq, "Multi"=nMult))
}
