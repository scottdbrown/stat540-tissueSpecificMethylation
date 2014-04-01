#-------------------------------------------------------------------------------
# Normalize methylation data
# Jessica Lee
# Date created: March 30, 2014
# Last edit: April 1, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(wateRmelon)
library(ggplot2)

# Read in data
load("450kMethylationData_probeLevel_nuke.RData")
load("450kMethylationData_geneLevelAverage_nuke.RData")
load("450kMethylationData_geneLevelPromoterAverage_nuke.RData")

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
saveMultiPlot(densPlot, dataAlias, "beta-density-alias-both-norm.png", 
    					width = 1000, height = 1000)

# Convert Beta values to M values
MNorm <- lapply(betaNorm, beta2m)

# Violin plot MNorm
MVioPlot <- lapply(MNorm, function(data){
	tmp <- melt(data)
	colnames(tmp) <- c("gene", "sample", "exp")
	
	ggplot(a, aes(sample, exp)) + 
		geom_violin(stat = "ydensity", position = position_dodge(0.5), 
	  	          fill = "cadetblue3") +
		stat_summary(fun.y = "median", geom='point') +
		ylab("M Value") +
		theme(axis.text.x  = element_text(angle=90))
})
MVioPlot <- Map(addTitle, MVioPlot, dataAlias, 
                "Distribution of M values (alias)")
showMultiPlot(MVioPlot)

# Save M value plot
saveMultiPlot(MVioPlot, dataAlias, "m-violin-alias.png",
              width = 1000, height = 1000)

# Assign normalized data globally
MNorm <- castGlobal(MNorm, "Nuke", "Norm")

# Save M values for differential methylation analysis
save(methylDatNorm,
     file = "450kMethylationData_probeLevel_norm.RData")
save(avgMethylByGeneNorm,
     file = "450kMethylationData_geneLevelAverage_norm.RData")
save(avgMethylByGenePromoterNorm,
     file = "450kMethylationData_geneLevelPromoterAverage_norm.RData")

