# BAM.fileTools.R -- wrappers to SAMTOOLS calls, for various BAM file manipulations


`BAM.sort` <- function( file, what=c("position", "readID"), index=TRUE, memory=2147000000) {

	if ( ! file.exists( file)) {
		cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
		return()
	}
	bamfile <- sub( "bam$", "sorted", file)

	what <- match.arg( what)
	sortOpt <- if ( what == "readID") " -n " else ""
	memValue <- as.integer(memory)
	if ( is.na( memValue)) memValue <- 500000000
	memOpt <- if ( is.null( memory)) "" else paste( " -m", as.integer(memValue))

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " sort ", sortOpt, memOpt, " ", file, " ", bamfile)
	cat( "\nSorting BAM file:  ", basename(file))
	system( cmdline)
	cat( "\nDone.")

	bamfile <- paste( bamfile, "bam", sep=".")
	if ( index) BAM.index( bamfile)

	return( bamfile)
}


`BAM.index` <- function( file) {

	if ( ! file.exists( file)) {
		cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
		return()
	}
	idxfile <- sub( "bam$", "bam.bai", file)
	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cat( "\nIndexing BAM file:  ", basename(file))
	cmdline <- paste( samtools, " index ", file)
	system( cmdline)
	cat( "\nDone.")
	return( idxfile)
}


`BAM.verifySorted` <- function( file, index=TRUE) {

	# make sure that all the names end up with the sort prefix
	file <- sub( "sorted.bam$", "bam", file)
	bamfile <- sub( "bam$", "sorted.bam", file)

	# we can allow the sorted file to exist without the original!!
	if ( ! file.exists( file)) {
		if ( ! file.exists( bamfile)) {
			cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
			return()
		}
	}

	# see if the 'sorted' one is there yet
	if ( ! file.exists( bamfile)) {
		ans <- BAM.sort( file, index=index)
	}

	# see if the 'index' is there too
	if (index) {
		idxfile <- sub( "bam$", "bam.bai", bamfile)
		if ( ! file.exists( idxfile)) ans <- BAM.index( bamfile)
	}

	return( bamfile)
}


`BAM.merge` <- function( files, newfile, index=TRUE) {

	for ( file in files) {
	    if ( ! file.exists( file)) {
		cat( "\nFile not found:  ", file, "\nFailed to find existing BAM file.")
		return()
	    }
	}

	if ( regexpr( "sorted.bam$", newfile) < 1) {
		newfile <- paste( newfile, "sorted.bam", sep=".")
	}

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " merge -f ", newfile, " ", paste( files, collapse=" "))
	cat( "\nMerging BAM file:  ", basename(files))
	system( cmdline)
	cat( "\nDone.  Created new merged file: ", basename(newfile), "\n")

	if ( index) BAM.index( newfile)
	return( newfile)
}


`BAM.indexFASTA` <- function( file) {

	if ( ! file.exists( file)) {
		cat( "\nFile not found:  ", file, "\nFailed to index FASTA file.")
		return()
	}
	newfile <- paste( file, "fai", sep=".")
	if ( file.exists( newfile)) return()

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " faidx ", file)
	cat( "\nIndexing FASTA file:  ", basename(file))
	system( cmdline)
	cat( "\nDone.")

	return( newfile)
}


`BAM.mpileup` <- function( files, seqID, fastaFile, start=NULL, stop=NULL, min.depth=3, max.depth=1000,
			min.gap.fraction=0.25, mpileupArgs="", summarize.calls=FALSE, verbose=TRUE) {

	N <- length( files)
	for ( i in 1:N) {
		files[i] <- BAM.verifySorted( files[i], index=TRUE)
	}
	fileArg <- files
	if ( N > 1) fileArg <- paste( files, collapse="  ")

	# SAMTOOLS expects these to be '.sorted.bam' files that already have '.sorted.bam.bai' files...  Check!
	allFound <- TRUE
	for ( f in files) {
		ff <- paste( f, "bai", sep=".")
		if ( ! file.exists( ff)) {
			cat( "\nBAM index file not found: ", ff)
			allFound <- FALSE
		}
	}
	if ( ! allFound) {
		cat( "\nCall to 'SAMTOOLS MPILEUP' failed.\n")
		return(NULL)
	}

	tmpFile <- tempfile()

	region <- seqID
	if ( ! is.null( start)) region <- paste( region, as.integer(start), sep=":")
	if ( ! is.null( stop)) region <- paste( region, as.integer(stop), sep="-")

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	cmdline <- paste( samtools, " mpileup -B -r ", region, " -f ", fastaFile, " -m ", min.depth, 
			" -d ", max.depth, " -F ", min.gap.fraction, " -L ", max.depth, 
			" ", mpileupArgs, " ", fileArg, " > ", tmpFile)
	if ( !verbose) cmdline <- paste( cmdline, "  2> /dev/null")

	if (verbose) cat( "\nGenerating Pileups for ", region, " of:  ", basename(files), "\n")
	system( cmdline)
	ans <- try( read.delim( tmpFile, header=FALSE, comment.char="", quote="", as.is=T), silent=TRUE)
	if ( class( ans) == "try-error") return( data.frame())
	file.remove( tmpFile)

	if ( N == 1) {
		colnames( ans)[1:6] <- c( "SEQ_ID", "POSITION", "REF_BASE", "DEPTH", "CALL_BASE", "CALL_SCORE")
	} else {
		# extract sampleIDs from the file names...
		sampleIDs <- sub( "(^.+)(\\.{1}?)(.+)", "\\1", basename(files))
		colnames( ans)[1:3] <- c( "SEQ_ID", "POSITION", "REF_BASE")
		colnames( ans)[4:ncol(ans)] <- paste( rep( c( "DEPTH", "CALL_BASE", "CALL_SCORE"), times=N),
							rep( sampleIDs, each=3), sep="_")
	}

	# we may want the base call summarized
	if ( summarize.calls) {
		if (verbose) cat( "\nSummarizing Pileups into final Base Calls..")
		calledBaseColumns <- grep( "CALL_BASE", colnames(ans))
		for ( baseColumn in calledBaseColumns) {

			rawBases <- ans[[ baseColumn]]
			callAns <- MPU.callBases( rawBases, ans$REF_BASE)
			ans[[ baseColumn]] <- callAns$call
			ans[[ baseColumn + 1]] <- MPU.callTableToString( callAns$depth.table)
			colnames(ans)[ baseColumn + 1] <- sub( "CALL_SCORE", "BASE_TABLE", colnames(ans)[baseColumn+1])
		}
	}

	if (verbose) cat( "\nDone.\nN_Lines: ", nrow(ans))
	return( ans)
}


`BAM.variantCalls` <- function( files, seqID, fastaFile, start=NULL, stop=NULL, 
				prob.variant=0.5, min.depth=3, max.depth=1000, min.gap.fraction=0.25,
				mpileupArgs="", vcfArgs="", geneMap=getCurrentGeneMap(), verbose=TRUE) {

	N <- length( files)
	fileArg <- files
	if ( N > 1) fileArg <- paste( files, collapse="  ")

	tmpFile <- tempfile()

	region <- seqID
	if ( ! is.null( start)) region <- paste( region, as.integer(start), sep=":")
	if ( ! is.null( stop)) region <- paste( region, as.integer(stop), sep="-")

	samtools <- Sys.which( "samtools")
	if ( samtools == "") stop( "Executable not found on search path:  'samtools'")
	bcftools <- Sys.which( "bcftools")
	if ( samtools == "") stop( "Executable not found on search path:  'bcftools'")

	cmdline <- paste( samtools, " mpileup -B -v -u -t DP -r ", region, " -f ", fastaFile, 
			" -d ", max.depth, " -m ", min.depth, " -F", min.gap.fraction,
			" -L", max.depth, mpileupArgs, "  ", fileArg, 
			" | ", bcftools, " call -mv ", " -p ", prob.variant, " ", 
			vcfArgs, " -O v -o ", tmpFile)
	if (verbose) {
		cat( "\nGenerating Variant Calls for ", region, " of:  ", basename(files), "\n")
		cat( "Command Line:   ", cmdline, "\n")
	}
	system( cmdline)
	ans <- try( read.delim( tmpFile, header=FALSE, comment.char="#", as.is=T), silent=TRUE)
	file.remove( tmpFile)

	if ( class( ans) == "try-error") {
		return(data.frame())
	} else {

		colnames( ans)[1:9] <- c( "SEQ_ID", "POSITION", "GENE_ID", "REF_BASE", "ALT_BASE", "QUAL", "FILTER", "INFO", "FORMAT")

		# there could be 'N's in the reference, that are of no use to us...
		isN <- which( ans$REF_BASE == "N")
		if ( length(isN) > 0) ans <- ans[ -isN, ]

		# extract sampleIDs from the file names...
		sampleIDs <- sub( "(^.+)(\\.{1}?)(.+)", "\\1", basename(files))
		colnames( ans)[10:(10+length(sampleIDs)-1)] <- sampleIDs

		# let's push the GeneID in...
		geneMap <- dropAntiSenseGenes( geneMap)
		gmap <- subset( geneMap, SEQ_ID == seqID)
		ptrs <- findInterval( ans$POSITION, gmap$POSITION)
		# tiny chance that no gene at very front edge of this map
		ptrs[ ptrs < 1] <- 1
		ans$GENE_ID <- gmap$GENE_ID[ ptrs]

		# the POSITION is relative to the size of the Indel, so adjust to point at the true culprit
		#nBase <- nchar( ans$REF_BASE)
		#isIndel <- which( nBase > 1)
		#newPos <- ans$POSITION[isIndel] + nBase[isIndel] - 1
		#ans$POSITION[isIndel] <- newPos

		return( ans)
	}
}
