#-------------------------------------------------------------------------------
# Differential methylation analysis
# Jessica Lee
# Date created: April 1, 2014
# Last edit: April 1, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(wateRmelon)
library(limma)

# Read in data
load("450kMethylationData_probeLevel_norm.RData")
load("450kMethylationData_geneLevelAverage_norm.RData")
load("450kMethylationData_geneLevelPromoterAverage_norm.RData")
load("450kMethylationData_meta_nuke.RData")

# Get me variable names
names <- ls()
datasets <- names[grep("Meta", names, invert = TRUE)] # data only
meta <- get(names[grep("Meta", names)]) # meta only

# Dataset alias
dataAlias <- c("gene", "promoter", "probe")
names(dataAlias) <- datasets

# Load helper functions
source("helpers.R")

#-------------------------------
# Functions
#-------------------------------
fitAndHit <- function(data, design, ...) {
	dmFit <- lmFit(data, design)
	dmEbFit <- eBayes(dmFit)
	topTable(dmEbfit, ...)
}

fitAndHitAll <- function(data, formula, meta, ...) {
	desMat <- model.matrix(formula, meta)	
	lapply(data, fitAndHit, design = desMat, ...)
}

#-------------------------------
# Workspace
#-------------------------------
# Limma throws an error for some probes that have 
# normalized values of -Inf or Inf
# Nuke probes that have -Inf or Inf
nonInfRows <- apply(methylDatNorm, 1, function(row){ 
	all((row != -Inf) & (row != Inf)) 
})
length(which(nonInfRows == FALSE))
methylDatRmInf <- methylDatNorm[nonInfRows, ]
dim(methylDatRmInf)
dim(methylDatNorm)

# Save RmInf as new data so we can lapply in next step
newData <- gsub("methylDatNorm", "methylDatRmInf", datasets)

# Differential methylation analysis with limma
# Get design matrix
desMat <- model.matrix(~ cellTypeSimple, meta)

# Fit linear model, test that cell types doesn't matter
# Get hits
hit <- lapply(getData(newData), fitAndHit, design = desMat, 
             coef = grep("cellType", colnames(coef(dmFit))), 
             n = Inf)

# One shot method to test various models
hit <- fitAndHitAll(getData(newData), 
                    ~ cellTypeSimple, meta,
                    coef = grep("cellType", colnames(coef(dmFit))), 
                    n = Inf)

# Save top hits
hit <- castGlobal(hit, "Norm|RmInf", "Hit")

# Save data frame
save(methylDatHit,
     file = "450kMethylationData_probeLevel_hit.RData")
save(avgMethylByGeneHit,
     file = "450kMethylationData_geneLevelAverage_hit.RData")
save(avgMethylByGenePromoterHit,
     file = "450kMethylationData_geneLevelPromoterAverage_hit.RData")
