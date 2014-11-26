# pipe.FileCleanup.R

# do the file combining, renaming, deleting temps, etc., after the alignment pipeline has run.

`pipe.FileCleanup` <- function( sampleID,  optionsFile="Options.txt", verbose=FALSE) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nStarting 'File Cleanup' for Sample:     ", sampleID, "\n\n")
	}

	resultsPath <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# delete files we don't need
	if (verbose) cat("\nDeleting temp files...")
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fastq", sep="."))
	file.delete( paste( sampleID, "not", c("ribo","genomic","splice"), "fastq.gz", sep="."))

	# see if any BAM conversions went longer than their genomic parts...
	if (verbose) cat("\nDeleting BAM conversion log files...")
	file.delete( paste( sampleID, "convertRiboBAM.log.txt", sep="."))
	file.delete( paste( sampleID, "convertingRiboBAM.done", sep="."))
	file.delete( paste( sampleID, "convertSpliceBAM.log.txt", sep="."))
	file.delete( paste( sampleID, "convertingSpliceBAM.done", sep="."))

	if ( verbose) cat( "\n...Done.\n")
}

