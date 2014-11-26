# pipe.RoundRobin.R

# run a N-way round robin comparison of the differential expression between a set of samples

`pipe.RoundRobin` <- function( sampleIDset, speciesID, annotationFile="Annotation.txt",
				optionsFile="Options.txt", 
				useMultiHits=TRUE, results.path=NULL,  folderName="", 
				groupColumn="Group", colorColumn="Color",
				altGeneMap=NULL, altGeneMapLabel=NULL, targetID=NULL, Ngenes=100, 
				geneColumnHTML=if (speciesID %in% c("Hs_grc","MacMu","Mmu_grc")) "NAME" else "GENE_ID", 
				keepIntergenics=FALSE, verbose=TRUE, label="", 
				doDE=TRUE, PLOT.FUN=NULL, ...)
{

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'RoundRobinCompare' on Sample Set: \n")
		print(sampleIDset)
		cat("\n", label, "\n\nUsing results from Species:  ", speciesID,"\n")
	}

	if ( is.null( targetID)) targetID <- getOptionValue( optionsFile, "targetID", notfound="HsPf", verbose=F)
	setCurrentTarget( targetID)

	# build the comparison data
	roundRobinCompare( sampleIDset, speciesID, annotationFile=annotationFile,
		optionsFile=optionsFile, 
		useMultiHits=useMultiHits, results.path=results.path, folderName=folderName, 
		groupColumn=groupColumn, colorColumn=colorColumn, Ngenes=Ngenes,
		altGeneMap=altGeneMap, altGeneMapLabel=altGeneMapLabel, 
		geneColumnHTML=geneColumnHTML, 
		keepIntergenics=keepIntergenics, verbose=verbose, label=label, 
		doDE=doDE, PLOT.FUN=PLOT.FUN, ...)

	geneTableToHTMLandPlotsCleanup()

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'RoundRobinCompare' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n", label, "\n")
	}

	return()
}

