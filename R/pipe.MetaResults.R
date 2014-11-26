# pipe.MetaResults.R

# run the 3 different DE tools and combine their results

`pipe.MetaResults` <- function( sampleIDset, speciesID="Pf3D7", annotationFile="Annotation.txt",
				optionsFile="Options.txt", 
				useMultiHits=TRUE, results.path=NULL,  folderName="", 
				groupColumn="Group", colorColumn="Color", average.FUN=sqrtmean, 
				tools=c("RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL,
				Ngenes=50, geneColumnHTML=if (speciesID %in% c("Hs_grc","MacMu","Mmu_grc")) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=TRUE, label="", doDE=TRUE, nFDRsimulations=0,
				PLOT.FUN=NULL, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'MetaResults' on Sample Set: \n")
		print(sampleIDset)
		cat("\n", label, "\n\nUsing results from Species:  ", speciesID,"\n")
	}

	# set up for this species...
	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	if ( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)
	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".")
	}
	
	# first call all 4 DE tools
	if ( doDE) {

		toolFuncList <- vector( mode="list")
		if ( "RoundRobin" %in% tools) toolFuncList <- c( toolFuncList, pipe.RoundRobin)
		if ( "EdgeR" %in% tools) toolFuncList <- c( toolFuncList, pipe.EdgeR)
		if ( "DESeq" %in% tools) toolFuncList <- c( toolFuncList, pipe.DESeq)
		if ( "RankProduct" %in% tools) toolFuncList <- c( toolFuncList, pipe.RankProduct)
		if ( "SAM" %in% tools) toolFuncList <- c( toolFuncList, pipe.SAM)

		multicore.lapply( toolFuncList, FUN=function(x) x( sampleIDset, speciesID, annotationFile, optionsFile,
					useMultiHits=useMultiHits, results.path=results.path, 
					groupColumn=groupColumn, colorColumn=colorColumn,
					folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
					Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
					targetID=targetID, verbose=verbose, label=label, PLOT.FUN=PLOT.FUN, ...))
	
		#if ( "RoundRobin" %in% tools) pipe.RoundRobin( sampleIDset, speciesID, annotationFile, optionsFile,
			#useMultiHits=useMultiHits, results.path=results.path, 
			#groupColumn=groupColumn, colorColumn=colorColumn,
			#folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
			#Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
			#targetID=targetID, verbose=verbose, label=label, PLOT.FUN=PLOT.FUN, ...)
	
		#if ( "EdgeR" %in% tools) pipe.EdgeR( sampleIDset, speciesID, annotationFile, optionsFile,
			#useMultiHits=useMultiHits, results.path=results.path, 
			#groupColumn=groupColumn, colorColumn=colorColumn,
			#folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
			#Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
			#targetID=targetID, verbose=verbose, label=label, PLOT.FUN=PLOT.FUN, ...)
	
		#if ( "DESeq" %in% tools) pipe.DESeq( sampleIDset, speciesID, annotationFile, optionsFile,
			#useMultiHits=useMultiHits, results.path=results.path, 
			#groupColumn=groupColumn, colorColumn=colorColumn,
			#folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
			#Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
			#targetID=targetID, verbose=verbose, label=label, PLOT.FUN=PLOT.FUN, ...)
	
		#if ( "RankProduct" %in% tools) pipe.RankProduct( sampleIDset, speciesID, annotationFile, optionsFile,
			#useMultiHits=useMultiHits, results.path=results.path, 
			#groupColumn=groupColumn, colorColumn=colorColumn,
			#folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
			#Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
			#targetID=targetID, verbose=verbose, label=label, PLOT.FUN=PLOT.FUN, ...)
	
		#if ( "SAM" %in% tools) pipe.SAM( sampleIDset, speciesID, annotationFile, optionsFile,
			#useMultiHits=useMultiHits, results.path=results.path, 
			#groupColumn=groupColumn, colorColumn=colorColumn,
			#folderName=folderName, altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel,
			#Ngenes=Ngenes, geneColumnHTML=geneColumnHTML, keepIntergenics=keepIntergenics,
			#targetID=targetID, verbose=verbose, label=label, PLOT.FUN=PLOT.FUN, ...)

	} # end of 'doDE'

	# now we can do meta analysis
	if (verbose) cat( "\nAll DE steps done.\n\nStarting 'MetaResults' phase...\n")

	metaPath <- file.path( results.path, "MetaResults", paste( prefix, folderName, sep="."))
	if ( ! file.exists( metaPath)) dir.create( metaPath, recursive=TRUE)
	pngPath <- "./pngPlots"

	flatSamples <- sort( unique( unlist( sampleIDset)))
	annT2 <- subset( annT, SampleID %in% flatSamples)
	myGrps <- sort( unique( annT2[[ groupColumn]]))
	for ( grp in myGrps) {

		cat( "\n\nDoing MetaResults on:    ", grp)

		out <- metaResults( targetGroup=grp, results.path=results.path, speciesID=speciesID,
				geneColumn="GENE_ID", subfolderName=folderName, tools=tools,
				altGeneMapLabel=altGeneMapLabel,
				rank.average.FUN=average.FUN, value.average.FUN=mean, 
				keepIntergenics=keepIntergenics, topFolder=NULL, 
				other.DE.files=NULL, missingGenes="na", nFDRsimulations=nFDRsimulations)

		# make the text file version
		outDF <- out
		if ( speciesID == "Hs_grc") outDF <- addHumanIDterms( outDF)
		if ( speciesID %in% c( "Pf3D7", "PbANKA", "Py17X", "PvSal1", "PyYM")) outDF <- addOrigIDterms( outDF)

		fileout <- paste( grp, prefix, "Meta.Ratio.txt", sep=".")
		if ( ! is.null( altGeneMapLabel)) {
			fileout <- paste( grp, prefix, altGeneMapLabel, "Meta.Ratio.txt", sep=".")
		}
		fileout <- file.path( metaPath, fileout)
		write.table( outDF, fileout, sep="\t", quote=F, row.names=F)
		rm( outDF)

		# simplify the names?
		fullGname <- out$GENE_ID
		if ( geneColumnHTML != "GENE_ID") {
			gmap <- getCurrentGeneMap()
			wh <- match( fullGname, gmap$GENE_ID, nomatch=0)
			fullGname[ wh > 0] <- gmap[[ geneColumnHTML]][wh]
			out$GENE_ID <- fullGname
		}

		if ( speciesID %in% c( "Pf3D7", "PbANKA", "Py17X", "PvSal1", "PyYM")) out <- addOrigIDterms( out)

		# special mods for altGeneMap...
		Nshow <- Ngenes
		title1 <- paste( "MetaResults: &nbsp;  Genes most up-regulated in group: &nbsp;  ", grp)
		title2 <- paste( "MetaResults: &nbsp;  Genes most down-regulated in group: &nbsp;  ", grp)
		if ( ! is.null( altGeneMapLabel)) {
			title1 <- paste( "MetaResults: &nbsp;  ", altGeneMapLabel, 
					"  most up-regulated in group: &nbsp;  ", grp)
			title2 <- paste( "MetaResults: &nbsp;  ", altGeneMapLabel, 
					"  most down-regulated in group: &nbsp;  ", grp)
			# for plots of varGenes we need to fudge a few items...
			if ( regexpr( "vargene|vardomain", tolower(altGeneMapLabel)) > 0) {
				out <- cbind( "DOMAIN_ID"=fullGname, out)
				out$GENE_ID <- sub( "::.*", "", fullGname)
				fullGname <- out$GENE_ID
				geneColumnHTML <- "GENE_ID"
				if ( "ORIG_ID" %in% colnames(out)) out$ORIG_ID <- gene2OrigID( out$GENE_ID)
			}
		}

		out1 <- out
		out1 <- out1[ out1$LOG2FOLD > 0, ]
		N <- min( nrow( out1), Ngenes)
		htmlout <- sub( "Ratio.txt", "UP.html", fileout)
		metaResultsToHTML( out1, htmlout, title1, maxRows=N, linkColumnName="GENE_ID",
				linkPaths=pngPath)
		rm( out1)

		out2 <- out[ rev( 1:nrow(out)), ]
		out2 <- out2[ out2$LOG2FOLD < 0, ]
		N <- min( nrow( out2), Ngenes)
		htmlout <- sub( "Ratio.txt", "DOWN.html", fileout)
		metaResultsToHTML( out2, htmlout, title2, maxRows=N, linkColumnName="GENE_ID",
				linkPaths=pngPath)
		rm( out2)
	}

	# copy all the gene plots to this new results location
	if ( is.null(PLOT.FUN) || is.function(PLOT.FUN)) {
		if (verbose) cat( "\nCopying plots.. ")
		for ( folder in c( "RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq")) {
			pathFrom <- file.path( results.path, folder, paste( prefix, folderName, sep="."), 
						"pngPlots")
			if ( ! file.exists( pathFrom)) next
	
			cmdLine <- paste( " cp -r ", pathFrom, "  ", metaPath)
			system( cmdLine)
			if (verbose) cat( "  ", folder)
		}
		if (verbose) cat( "  Done.\n")
	}

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'MetaResults' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n", label, "\n")
	}
	return()
}

