#####################################
#### Expression Data Acquisition ####
#### Author(s): Scott Brown
#### Date Created: March 22, 2014
#### Last Edited by: Scott Brown
#### on: March 22, 2014
#####################################


# GSE Numbers for the data in the study:
# GSE30652 for Illumina HT12v3 Gene Expression
# GSE30653 for Illumina Infinium 27K DNA Methylation
# GSE31848 for Illumina Infinium 450K DNA Methylation

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"expression/")
setwd(workingDir)

library(GEOquery)

geoDat <- getGEO("GSE30652")

expDat <- as.data.frame(exprs(geoDat[[1]]))
expMeta <- pData(phenoData(geoDat[[1]]))
# the 'title' variable in methylMeta corresponds to the 450K sample name in the supplementary table 1
# the 'geo_accession' variable in methylMeta corresponds to the sample name in the methylDat data.frame


# Just incase I need to go back to the raw geo download, saved in ~/GSAT @ phage.bcgsc.ca
#save(geoDat, file="rawGeoMethyl450k.RData")

meta <- read.delim("../metadata.tsv")
gsmToKeep <- expMeta$geo_accession[expMeta$title %in% meta$HT12v3.Sample.Name]

expDat <- expDat[,gsmToKeep]
expMeta <- expMeta[gsmToKeep,]
# Now save the methylation data.frame and the metadata data.frame
save(expDat, expMeta, file="HT12v3_expression.RData")
