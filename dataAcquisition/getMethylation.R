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

# edit working directory path if needed:
workingDir <- "~/stat540-tissueSpecificMethylation/"
setwd(workingDir)

library(GEOquery)
library(IlluminaHumanMethylation450k.db)


#Following command crashed computer after 2.5 hours, not enough RAM.
geoDat <- getGEO("GSE31848")
