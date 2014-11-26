# pipe.PlotSNP.R

`pipe.PlotSNP` <- function( sampleIDs, seqID, position, tailWidth=100, groups=sampleIDs,
				annotationFile="Annotation.txt", optionsFile="Options.txt", results.path=NULL, 
				asPNG=FALSE, pngPath="pngPath", label="", SNPtablePath="~/SNPs/", 
				mf=NULL, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	if ( ! is.null(SNPtablePath)) SNP_curSNPpath <<- SNPtablePath

	if (asPNG) {
		if ( ! file.exists( pngPath)) dir.create( pngPath, recursive=T, showWarnings=F)
	}
	bamfiles <- file.path( results.path, "align", paste( sampleIDs, "genomic.sorted.bam", sep="."))
	vcffiles <- file.path( results.path, "VariantCalls", sampleIDs, 
				paste( sampleIDs, seqID, "VCF.txt", sep="."))

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID)
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples))
	}
	doMulti <- ( length( doSamples) > 1)

	gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & position >= POSITION & position <= END)
	if ( nrow( gmap) > 1) {
		bestOne <- findInterval( position, gmap$POSITION)
		gmap <- gmap[ bestOne, ]
	}

	if ( dev.cur() < 2) X11( bg="white", width=14, height=10)

	if (doMulti) {
		multiSample.plotSNP( position, seqID, sampleSet=sampleIDs, bamfileSet=bamfiles, 
				vcffileSet=vcffiles, fastaFile=fastaFile, groupSet=groups, 
				tailWidth=tailWidth, gmap=gmap, mf=mf, ...)
	} else {
		sid <- sampleIDs[1]
		bamfile <- bamfiles[1]
		plotSNP( position=position, seqID=seqID, sampleID=sid, bamfile=bamfiles, 
				vcffile=vcffiles, fastaFile=fastaFile,
				tailWidth=tailWidth, mode="single", gmap=gmap, label=label, ...)
	}
	if ( asPNG) {
		geneID <- shortGeneName( gmap$GENE_ID[1], keep=1)
		plotfile <- paste( geneID, position, "png", sep=".")
		plotfile <- file.path( pngPath, plotfile)
		dev.print( png, plotfile, width=1000, height=700, "bg"="white")
	}
	return()
}


`pipe.PlotGeneSNPs` <- function( sampleIDs, geneID, tailWidth=100, groups=sampleIDs,
				annotationFile="Annotation.txt", optionsFile="Options.txt", results.path=NULL, 
				asPNG=FALSE, pngPath="pngPath", label="", SNPtablePath="~/SNPs/", 
				mf=NULL, ...) {

	# get needed paths, etc. from the options file
	optT <- readOptionsTable( optionsFile)
	curSpecies <- getCurrentSpecies()
	geneMap <- getCurrentGeneMap()
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	}
	fastaFile <- getOptionValue( optT, "genomicFastaFile", verbose=F)

	if ( ! is.null(SNPtablePath)) SNP_curSNPpath <<- SNPtablePath

	# map from a gene name to the seq, position that the SNP tool wants
	gmap <- subset.data.frame( geneMap, GENE_ID == geneID)
	if ( nrow( gmap) < 1) {
		cat( "\nNo Gene found that matches: ", geneID)
		return()
	}
	seqID <- gmap$SEQ_ID[1]
	midpt <- (gmap$POSITION[1] + gmap$END[1]) / 2


	if (asPNG) {
		if ( ! file.exists( pngPath)) dir.create( pngPath, recursive=T, showWarnings=F)
	}
	if ( dev.cur() < 2) X11( bg="white", width=14, height=10)

	# get the data we need
	bamfiles <- file.path( results.path, "align", paste( sampleIDs, "genomic.sorted.bam", sep="."))
	vcffiles <- file.path( results.path, "VariantCalls", sampleIDs, 
				paste( sampleIDs, seqID, "VCF.txt", sep="."))

	# make sure we know all these sampleIDs
	annT <- readAnnotationTable( annotationFile)
	doSamples <- intersect( sampleIDs, annT$SampleID)
	if ( length( doSamples) < 1) {
		cat( "\nNo Samples found that match annotation: ", annT$SampleID)
		return()
	}
	if ( length( doSamples) < length( sampleIDs)) {
		cat( "\nSome Samples not found: ", setdiff( sampleIDs, doSamples))
	}
	doMulti <- ( length( doSamples) > 1)

	if (doMulti) {
		multiSample.plotSNP( position=midpt, seqID=seqID, sampleSet=sampleIDs, bamfileSet=bamfiles, 
				vcffileSet=vcffiles, fastaFile=fastaFile, groupSet=groups, 
				tailWidth=tailWidth, gmap=gmap, show.legends=FALSE, mf=mf, ...)
	} else {
		sid <- sampleIDs[1]
		bamfile <- bamfiles[1]
		plotSNP( position=midpt, seqID=seqID, sampleID=sid, bamfile=bamfiles, 
				vcffile=vcffiles, fastaFile=fastaFile,
				tailWidth=tailWidth, mode="single", gmap=gmap, label=label, show.legends=FALSE, ...)
	}
	if ( asPNG) {
		geneID <- shortGeneName( gmap$GENE_ID[1], keep=1)
		plotfile <- paste( geneID, "FullLength", "png", sep=".")
		plotfile <- file.path( pngPath, plotfile)
		dev.print( png, plotfile, width=1000, height=700, "bg"="white")
	}
	return()
}

