######################################
#### Methylation Data Acquisition ####
#### Author(s): Scott Brown
#### Date Created: March 21, 2014
#### Last Edited by: Scott Brown
#### on: March 21, 2014
######################################


# GSE Numbers for the data in the study:
# GSE30652 for Illumina HT12v3 Gene Expression
# GSE30653 for Illumina Infinium 27K DNA Methylation
# GSE31848 for Illumina Infinium 450K DNA Methylation

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"dataAcquisition")
setwd(workingDir)

library(GEOquery)

geoDat <- getGEO("GSE31848")

methylDat <- as.data.frame(exprs(geoDat[[1]]))
methylMeta <- pData(phenoData(geoDat[[1]]))
# the 'title' variable in methylMeta corresponds to the 450K sample name in the supplementary table 1
# the 'geo_accession' variable in methylMeta corresponds to the sample name in the methylDat data.frame


# Just incase I need to go back to the raw geo download, saved in ~/GSAT @ phage.bcgsc.ca
#save(geoDat, file="rawGeoMethyl450k.RData")

meta <- read.delim("../metadata.tsv")
gsmToKeep <- methylMeta$geo_accession[methylMeta$title %in% meta$X450k.dataSample.Name]

methylDat <- methylDat[,gsmToKeep]
methylMeta <- methylMeta[gsmToKeep,]
# Now save the methylation data.frame and the metadata data.frame
save(methylDat, methylMeta, file="450kMethylationData_probeLevel.RData")
