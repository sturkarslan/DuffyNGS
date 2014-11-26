# pipe.ConsensusBaseCalls.R

`pipe.ConsensusBaseCalls` <- function( sampleID, geneID=NULL, seqID=NULL, position=NULL, end=NULL, 
				annotationFile="Annotation.txt", optionsFile="Options.txt", results.path=NULL,
				aaToo=TRUE) {
				
	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	exonMap <- getCurrentExonMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	genomicFastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	bamfile <- file.path( results.path, "align", paste( sampleID, "genomic.bam", sep="."))

	# specify by either gene or range
	if ( is.null( geneID)) {
		gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & end >= POSITION & position <= END)
		emap <- subset.data.frame( exonMap, GENE_ID %in% gmap$GENE_ID)
		geneID <- gmap$GENE_ID[1]
		geneID <- shortGeneName( geneID, keep=1)
	} else {
		gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
		emap <- subset.data.frame( exonMap, GENE_ID == geneID)
		if (is.null(position)) position <- gmap$POSITION[1]
		if (is.null(end)) end <- gmap$END[1]
		seqID <- gmap$SEQ_ID[1]
	}
	geneStrand <- gmap$STRAND[1]
	gProd <- gmap$PRODUCT[1]

	# load that portion of the PILEUPS
	curMPU <- BAM.mpileup( bamfile, seqID, genomicFastaFile, start=position, stop=end, summarize.calls=TRUE,
					verbose=FALSE)
	saveMPU <<- curMPU

	# add the reference genome
	genomicStr <- getFastaSeqFromFilePath( genomicFastaFile, seqID)
	curGenomeDNA <- strsplit( as.character(genomicStr), split="")[[1]]

	# gather the reference genome bases in this range
	allBases <- position : end
	genomeSNPtext <- genomeBaseText <- curGenomeDNA[ allBases]
	names(genomeSNPtext) <- names(genomeBaseText) <- allBases

	hasBaseCalls <- ( nrow( curMPU) > 0)
	if ( hasBaseCalls) {

		xLocs <- curMPU$POSITION

		# turn the base calls into an explicit matrix in order: genomic,A,C,G,T,N, indel
		flips <- MPU.callStringsToMatrix( curMPU$BASE_TABLE)
		rownames(flips) <- curMPU$POSITION
		saveFlips <<- flips

		# get the majority base at each SNP spot, if not a SNP it will be ',' (matches the reference
		snpTopBase <- apply( flips, MARGIN=1, function(x) colnames(flips)[ which.max(x)])
		names( snpTopBase) <- xLocs

		# the Indels are tougher, we need to get the actual bases for both the reference and the indels and figure if
		# its an insertion or deletion
		indelDetails <- getIndelDetails( flips, curMPU)
		if ( indelDetails$nIndels > 0) {
			isIndel <- indelDetails$who
			snpTopBase[isIndel] <- indelDetails$bases
		}
	
		# see what bases are different from genomic, and draw the different AA
		whoSNP <- base::which( snpTopBase != ",")
		if ( length( whoSNP) > 0) {
			myLocs <- as.integer( names(snpTopBase)[whoSNP])
			where <- base::match( myLocs, allBases)
			genomeSNPtext[ where] <- snpTopBase[ whoSNP]
		}

		# draw the protein amino acid letters. with and without SNP effect
		genomeAminoText <- ""
		if (aaToo) {
			ans <- convertGenomicBasesToCodingAminoAcids( seqID, position=position, end=end, 
					strand=geneStrand, dnaQuery=genomeSNPtext, genomeDNA=curGenomeDNA,
					geneMap=gmap, exonMap=emap)
			genomeAminoText <- ans$genomic
			names( genomeAminoText) <- allBases
		}
	} else {
		genomeAminoText <- ""
	}

	out <- list( "ref"=genomeBaseText, "calls"=flips, "dna.consensus"=genomeSNPtext, "aa.consensus"=genomeAminoText)
	return( out)
}

