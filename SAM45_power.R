# setwd('/Users/kameno/Downloads/T100logging')	# Set the directory where data and source codes is
source('SAM45.R')
#--------- Function to read the given logging data and to output power of each array
SAM45Power <- function( fname ){
	tmpdata <- readSAM45( fname )		# Read SAM45 logging file
	A_list  <- levels(tmpdata$scan$IF)	# SAM45 array name (A1, A2, … )
	num_array <- length(A_list)			# Number of arrays
	
	IF_index <- which( tmpdata$scan$IF == levels(tmpdata$scan$IF)[1] & (tmpdata$scan$scan=='OFF' | tmpdata$scan$scan=='ON') )				# Pick up number of time records
	PowerDF <- data.frame(mjd = tmpdata$scan$mjd[IF_index])	# Create Data Frame and Store MJD
	
	for(index in 1:num_array){			# Loop for arrays to extract spectra
		IF_index <- which( tmpdata$scan$IF == levels(tmpdata$scan$IF)[index] & (tmpdata$scan$scan=='OFF' | tmpdata$scan$scan=='ON') )			# Index to pick up the array data
		tmpArray <- tmpdata$spec[, IF_index]	# Spectra for the array
		tmpPower <- colSums(tmpArray) / dim(tmpArray)[1]	# Total power of the array
		PowerDF[A_list[index]] <- tmpPower		# Append the power in the data frame
	}
	return( PowerDF )
}

#--------- Function to Plot the power
PDFpower <- function( fname ){
	tempDF <- SAM45Power(paste(fname, ".dat", sep=""))
	pdf(paste(fname, ".pdf", sep=""))
	for( index in 2:length(names(tempDF)) ){
		plot(tempDF$mjd, tempDF[[index]], type='s', xlab='MJD', ylab='Relative Power', main=paste(fname, names(tempDF)[index], "Power Plot"))
	}
	dev.off()
}

#--------- Execute Power Plot for every file
fileList <- c("T20110104200116", "T20110104214301", "T20110506074419", "T20110520081323", "T20110520094634")
sapply(fileList, PDFpower)


#-------- Power Plot for T20110104200116.dat
# pdf("T20110104200116.pdf")
# T20110104200116 <- SAM45Power('T20110104200116.dat')
# for( index in 2:length(names(T20110104200116)) ){
# 	plot(T20110104200116$mjd, T20110104200116[[index]], type='s', xlab='MJD', main=paste("T20110104200116", names(T20110104200116)[index], "Power Plot"))
# }
# dev.off()

