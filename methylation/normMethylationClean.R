#-------------------------------------------------------------------------------
# Normalize methylation data
# With outliers included
# Jessica Lee
# Date created: April 2, 2014
# Last edit: April 3, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(wateRmelon)
library(ggplot2)

# Read in data
load("450kMethylationData_probeLevel_clean.RData")
load("450kMethylationData_geneLevelAverage_clean.RData")
load("450kMethylationData_geneLevelPromoterAverage_clean.RData")

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
# Workspace
#-------------------------------
# Quantile normalization on beta values
betaNorm <- lapply(getData(datasets), function(data){
	mat <- as.matrix(data)
	betaqn(mat)
})

# Convert Beta values to M values
MNorm <- lapply(betaNorm, beta2m)

# Assign normalized data globally
betaNorm <- castGlobal(betaNorm, "Clean", "CleanBetaNorm")
MNorm <- castGlobal(MNorm, "Clean", "CleanNorm")

# Save M values for differential methylation analysis
save(methylDatCleanNorm, methylDatCleanBetaNorm,
     file = "450kMethylationData_probeLevel_clean_norm.RData")
save(avgMethylByGeneCleanNorm, avgMethylByGeneCleanBetaNorm, 
     file = "450kMethylationData_geneLevelAverage_clean_norm.RData")
save(avgMethylByGenePromoterCleanNorm, avgMethylByGenePromoterCleanBetaNorm,
     file = "450kMethylationData_geneLevelPromoterAverage_clean_norm.RData")

#-------------------------------
# Plotting work
#-------------------------------
# Plot Beta Density
densDat <- lapply(betaNorm, getBetaAvg, meta = meta)
densPlot <- lapply(densDat, plotDensity, auto.key = TRUE, plot.points = FALSE, xlab = "Beta Value")
densPlot <- Map(addTitle, densPlot, dataAlias,
    						"Density of Beta Values by Cell Type (alias)")
showMultiPlot(densPlot)
saveMultiPlot(densPlot, dataAlias, "beta-density-alias-after-norm.png", 
    					width = 1000, height = 1000)

# Both before and after
before <- lapply(getData(datasets), getBetaAvg, meta = meta)
after <- lapply(betaNorm, getBetaAvg, meta = meta)
densDat <- Map(function(b, a){
		rbind(cbind(b, norm="before"), cbind(a, norm="after"))
	}, before, after)
densPlot <- lapply(densDat, function(x, ...){
		densityplot(~ beta | norm, x, group = cellType, ...)
	}, auto.key = TRUE, plot.points = FALSE, xlab = "Beta Value")
densPlot <- Map(addTitle, densPlot, dataAlias,
    						"Density of Beta Values by Cell Type (alias)")
showMultiPlot(densPlot)

# Save beta density plot
saveMultiPlot(densPlot, dataAlias, "beta-density-alias-both-norm-clean.png", 
    					width = 1000, height = 1000)

