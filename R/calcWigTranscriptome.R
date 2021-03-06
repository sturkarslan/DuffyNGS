# calcWigTranscriptome.R

# turn a Wiggles object into a transcriptome.   Converting from WB to 'wiggles' data strucure...

`calcWigTranscriptome` <- function( WIG, geneMap=NULL, useBothStrands=FALSE, keepIntergenics=FALSE, 
			byExon=FALSE, fileout=NA, verbose=!interactive() ) {

	setCurrentSpecies( WIG$Info$Species)

	if ( is.null( geneMap)) {
		geneMap <- getCurrentGeneMap()
	} else {
		sameSpecies <- all( unique.default( geneMap$SEQ_ID) %in% getCurrentSeqMap()$SEQ_ID)
		if ( ! sameSpecies) stop( "calcTranscriptome:  given 'geneMap' does not go with current Species")
	}

	cat( "\n\nCalculating Transcriptome:");
	cat( "\nSpecies:         \t", getCurrentSpecies());
	cat( "\nAlignment files: \n  ", paste( WIG$Info$FileName, collapse="\n  "), "\n", sep="")

	# allocate vectors for result
	geneSet <- gProdSet <- vector( length=nrow( geneMap))
	totCntSet <- sigmaSet <- rpkmSet <- strandSet <- vector( length=nrow( geneMap))
	totCntSetM <- sigmaSetM <- rpkmSetM <- strandSetM <- vector( length=nrow( geneMap))
	nBaseSet <- vector( length=nrow( geneMap))

	# pre evaluate what we can
	hasNonGenes <- ("REAL_G" %in% colnames( geneMap))
	totalWIGreads <- WIG_getTotalReadsForRPKM( WIG)
	readLength <- WIG$Info$ReadLength
	# catch the very rare case of no reads
	if ( totalWIGreads$Unique < 1) totalWIGreads$Unique <- 1
	if ( totalWIGreads$Multi < 1) totalWIGreads$Multi <- 1

	# visit every gene
	ngenes <- nTranscribed <- 0
	curSeq <- ""

	for( ig in 1:nrow( geneMap)) {

		# bypass this one ?
		if ( !keepIntergenics && hasNonGenes && ( geneMap$REAL_G[ig] == FALSE)) next

		gene <- geneMap$GENE_ID[ig]
		seqID <- geneMap$SEQ_ID[ig]
		if ( seqID != curSeq) {
			wiggleChunk <- WIG_getWigglesOneSeq( WIG, seqID)
			curSeq <- seqID
		}

		ans <- calcWigTranscriptOneGene( wiggleChunk, gene, gmapPtr=ig, geneMap=geneMap, 
				totalReads=totalWIGreads, readLength=readLength,
				useBothStrands=useBothStrands, byExon=byExon)
		if ( is.null(ans)) next

		# load the vectors
		ngenes <- ngenes + 1
		geneSet[ ngenes] <- gene
		gProdSet[ ngenes] <- geneMap$PRODUCT[ig]
		rpkmSet[ ngenes] <- ans$rpkm
		totCntSet[ ngenes] <- ans$rawReads
		sigmaSet[ ngenes] <- ans$sigma
		strandSet[ ngenes] <- ans$strandness
		rpkmSetM[ ngenes] <- ans$rpkm.Multi
		totCntSetM[ ngenes] <- ans$rawReads.Multi
		sigmaSetM[ ngenes] <- ans$sigma.Multi
		strandSetM[ ngenes] <- ans$strandness.Multi
		nBaseSet[ ngenes] <- ans$nBases

		if (ngenes %% 1000 == 0) cat( "\n", ngenes, "\t", gene,"\t", ans$rawReads)

		if ( ans$rpkm.Multi > 0) nTranscribed <- nTranscribed + 1
	}

	out <- data.frame( geneSet, gProdSet, 
			rpkmSetM, totCntSetM, sigmaSetM, strandSetM, 
			rpkmSet, totCntSet, sigmaSet, strandSet, 
			nBaseSet, stringsAsFactors=FALSE)
	colnames( out) <- c("GENE_ID", "PRODUCT", 
			"RPKM_M", "READS_M", "SIGMA_M", "STRAND_M", 
			"RPKM_U", "READS_U", "SIGMA_U", "STRAND_U", 
			"N_EXON_BASES")
			
	# now trim to the true size
	out <- out[ 1:ngenes, ]

	# now sort into intensity order
	ord <- base::order( out$RPKM_M, decreasing=TRUE)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)
	
	if ( ! is.na( fileout)) {

		if ( getCurrentSpecies() == "Hs_grc") out <- addHumanIDterms( out)
		if ( getCurrentSpecies() %in% c( "Pf3D7", "PbANKA", "Py17X", "PyYM", "PvSal1")) out <- addOrigIDterms( out)

		write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)
		cat( "\nWrote Transcriptome file:  \t", fileout)
	}
	
	if (verbose) {
		cat( "\nN_Gene regions processed:         \t", nrow(out))
		cat( "\nN_Genes regions with RPKM > 0:     \t", nTranscribed, "\n")
	}

	return( out)
}

