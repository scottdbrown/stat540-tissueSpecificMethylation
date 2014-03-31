#-------------------------------------------------------------------------------
# Methylation Clean-up
# Jessica Lee
# Date created: March 24, 2014
# Last edit: March 30, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(reshape2)

# Read in data
load("450kMethylationData_probeLevel.RData")
load("450kMethylationData_geneLevelAverage.RData")
load("450kMethylationData_geneLevelPromoterAverage.RData")

# Get me variable names
names <- ls()

# Load helper functions
source("helpers.R")

#-------------------------------
# Functions
#-------------------------------
# Do we have NA's in the dataset?
# If so, how many do we have?
checkNA <- function(frame) {
	hasNA <- is.na(frame)
	if (any(hasNA)) {
		return (sum(hasNA))
	}
	else {
		return (FALSE)
	}
}

# Which row has NA? 
# Return a rows that have NA
whichRowNA <- function(frame) {
	findNA <- function(row){
		any(is.na(row))
	}
	return(frame[apply(frame, 1, findNA), ])
}

# Nuke rows with NA
# We don't want to remove samples. 
# Just nuke probes/genes missing values
nukeNA <- function(frame) {
	return (na.omit(frame))
}

# String match required info so we can parse group info properly
reparse <- function(infoToGrab, frame) {
	# All meta data has colon and a space attached to it
	searchTerm <- paste(infoToGrab, ": ", sep = "")

	# Grep the info I'm looking for in given column
	fetchIndex <- function(col) { 
		grep(searchTerm, col)
	}
	indexMap <- sapply(frame, fetchIndex)

	# Melt away empty columns
	indexMap <- melt(indexMap)
	colnames(indexMap) <- c("index", "colName")
	# return(indexMap)

	# Re-order by index
	indexMap <- indexMap[order(indexMap$index), ]
	row.names(indexMap) <- indexMap$index

	# Grab actual info
	fetchValue <- function(row) {
		# Get index
		index <- as.numeric(row["index"])

		# Get value
		val <- frame[index, row["colName"]]

		# Strip heading
		val <- gsub(searchTerm, "", val)
		return(val)
	}
	values <- apply(indexMap, 1, fetchValue)

	# Fill in empty rows
	mapInfo <- function(index) {
		if(index %in% indexMap$index) {
			return(values[which(index == indexMap$index)])
		} 
		else {
			return(NA)
		}
	}
	info <- as.vector(sapply(1:nrow(frame), mapInfo))
	return(info)
}

# Swap the long name of stem cells based on given hash
swapName <- function(x, swapHash){
	if (x %in% names(swapHash)){
		swapHash[[x]]
	} else {
		x
	}
}

#-------------------------------
# Workspace
#-------------------------------
# Peek
lapply(getData(names), head, n = 5)

# How big
lapply(getData(names), dim)

# Is there NA in the dataset?
lapply(getData(names), checkNA)
# >> Yup

# Find out which row has NA (if interested)
NArows <- lapply(getData(names), whichRowNA)

# Nuke NA
dataNoNA <- lapply(getData(names), nukeNA)

# Recast cleaned methyl data
methylDatClean <- dataNoNA$methylDat
avgMethylByGeneClean <- dataNoNA$avgMethylByGene
avgMethylByGenePromoterClean <- dataNoNA$avgMethylByGenePromoter

# Check dimension
lapply(getData(names), dim) # before nuke
lapply(dataNoNA, dim) # after nuke

# Get cleaned up metadata 
# (there's no NA in metadata so we don't have any omits)
cleanMeta <- dataNoNA$methylMeta

# Are column names decipherable?
colnames(cleanMeta)
# >> Somewhat

# Are columns clear enough to show different groups?
str(cleanMeta)
# >> Errrrr No. 

# For differential methylation analysis, I need group info.
# Peeking into metadata, looks like columns called 'characteristics_ch<xx>' 
# are the ones holding group info
# Peek into the characteristics columns
charMeta <- cleanMeta[ , grep("characteristic", colnames(cleanMeta))]
charMeta

# Reshape this factor hell
methylMetaClean <- with(cleanMeta,
                     		cbind(data.frame(
                     		      sid = as.character(title),
                              geo = as.character(geo_accession),
                              stringsAsFactors= FALSE),
                     			 		data.frame(
                              cellType = reparse("cell type", cleanMeta),
                              gender = reparse("gender", cleanMeta),
                              cellLine = reparse("cell line", cleanMeta),
                              tissue = reparse("tissue type", cleanMeta),
                              fetalAdult = reparse("fetal vs adult tissue", cleanMeta),
                              gestAge = reparse("gestational age", cleanMeta),
                              patID = reparse("patient id", cleanMeta),
                              passNum = reparse("passage number", cleanMeta))))

# Reduce name of stem cell type
swapHash <- list("ES", "iPS", "ES.parthenote")
names(swapHash) <- 
	levels(methylMetaClean$cellType)[grep("stem", 
	                                      levels(methylMetaClean$cellType))]                   
swapHash # check if this is right

# Swap the long name of stem cells
new <- sapply(as.character(methylMetaClean$cellType), swapName, 
              swapHash = swapHash)
new <- unname(new)
methylMetaClean <- cbind(methylMetaClean, cellTypeShort = new)


# Save reshaped table
save(methylDatClean, methylMetaClean, 
     file = "450kMethylationData_probeLevel_clean.RData")
save(avgMethylByGeneClean, 
     file = "450kMethylationData_geneLevelAverage_clean.RData")
save(avgMethylByGenePromoterClean,
     file = "450kMethylationData_geneLevelPromoterAverage_clean.RData")


