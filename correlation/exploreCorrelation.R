####################################################
#### Explore Expression/Methylation Correlation ####
#### Author(s): Scott Brown
#### Date Created: March 31, 2014
#### Last Edited by: Scott Brown
#### on: March 31, 2014
####################################################

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"correlation")
setwd(workingDir)

# load methylation data
load("../methylation/450kMethylationData_geneLevelAverage_clean.RData")

# load expression data
load("../expression/HT12v3_avgExpressionByGene.RData")