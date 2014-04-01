#-------------------------------------------------------------------------------
# Explore methylation data
# Jessica Lee
# Date created: March 24, 2014
# Last edit: March 31, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(RColorBrewer)
library(gplots)
library(hexbin)
library(lattice)
library(reshape2)

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

# Load color scheme
rdBu <- colorRampPalette(brewer.pal(n = 11, "RdBu"))

#-------------------------------
# Functions
#-------------------------------
# Get beta values
getBeta <- function(data, meta) {
	dens <- sapply(levels(meta$cellTypeShort), 
              	 function(type){
               	 rowMeans(data[ , 
               	 	        which(meta$cellTypeShort == type)])})
	dens <- melt(dens)
	colnames(dens) <- c("gene", "cellType", "beta")
	return(dens)
}

# Draw densityplot
plotDensity <- function(frame, ...) {
	return(densityplot(~ beta, frame, group = cellType, ...))
}

# Return heatmap data
getHeatMat <- function(data, meta) {
	# Matrix
	res <- as.matrix(data)

	# Check if order is right
	check <- all(order(meta$geo) == order(colnames(res)))
	
	# Fix order if not
	if (!check) {
		res <- res[ , order(colnames(res))]
	}

	# Rename columns for convenient reference
	colnames(res) <- with(meta, paste(geo, cellTypeShort, tissue, 
	                            	    cellLine, sep = "_"))
	return(res)
}

# Plot pearson correlation on heatmap
plotCor <- function(frame, ...) {
	frameCor <- cor(frame)
	heatmap(frameCor, scale = "none", col = rdBu(256), ...)
}

# Plot hex scatter
plotHex <- function(frame, ...) {
	return(hexplom(frame, ...))
}

# Examine outlier in scatter plot against other samples
examOutL <- function(frame, outL, nonOutL, ...) {
	# If no outlier group is set, get non-outlier groups
	if (missing(nonOutL)) {
		nonOutL <- which(colnames(frame) != outLName)
	} 

	# We don't want too many samples on scatter plot (SLOW)
	# So restrict to 5 samples
	if(length(nonOutL) > 5) {
		nonOutL <- sample(nonOutL, 5)
	}
	
	# Plot scatter
	cols <- which(colnames(frame) %in% c(outL, nonOutL))
	return(plotHex(frame[ , cols], ...))
}

#-------------------------------
# Workspace
#-------------------------------
# Density of Beta values
densDat <- lapply(getData(datasets), getBeta, meta = meta)
densPlot <- lapply(densDat, plotDensity, auto.key = TRUE, 
                   plot.points = FALSE, xlab = "Beta")
densPlot <- Map(addTitle, densPlot, dataAlias,
    						"Density of Beta Values by Cell Type (alias)")
showMultiPlot(densPlot)

# Save beta density plot
saveMultiPlot(densPlot, dataAlias, "beta-density-alias-before-norm.png", 
    					width = 1000, height = 1000)

# Gene-level pearson correlation
hDat <- lapply(getData(datasets), getHeatMat, meta = meta)
hPlot <- lapply(hDat, function(dat, ...){dev.new(); plotCor(dat, ...)}, 
                Rowv = NA, Colv = NA, margin = c(5, 15), 
                labCol = meta$geo)

hPlot <- lapply(hDat, function(dat, ...){dev.new(); plotCor(dat, ...)}, 
                margin = c(5, 15), 
                labCol = meta$geo) # With clustering
# >> looks like GSM867986, GSM867947, GSM867948, GSM867949 are scary batches

# Save Pearson correlation plot
hPlot <- Map(function(dat, alias){
							png(gsub("alias", alias, "cor-alias-before-norm.png"),
							    width = 1000, height = 1000)
							plotCor(dat, margin = c(c(5, 15)), labCol = methylMetaClean$geo)
							dev.off()
						 }, 
             hDat, dataAlias) # With clustering


# Examine possible outlier GSM867986 by plotting with its replicate GSM867987
outL <- c("GSM867986")

# Examine the outlier in a cellType context
# iPS outlier
iPSOutLPlot <- lapply(getData(datasets), examOutL, outL,
                      unlist(subset(meta, cellTypeShort == "iPS" & 
            											  !(geo %in% outL), select = geo)),
                      varname.cex = .8, axis.cex = 0.8)
iPSOutLPlot <- Map(addTitle, iPSOutLPlot, dataAlias, 
                   paste("Outlier", paste(outL, collapse = ", "), "at alias-level"))
showMultiPlot(iPSOutLPlot)

# Save plot
saveMultiPlot(iPSOutLPlot, dataAlias, "outlier-ips-alias.png", 
    					width = 1000, height = 1000)

# ES outliers
outL <- c("GSM867947", "GSM867948", "GSM867949")
ESOutLPlot <- lapply(getData(datasets), examOutL, outL,
                      unlist(subset(meta, cellTypeShort == "ES" & 
            											  !(geo %in% outL), select = geo)),
                      varname.cex = .8, axis.cex = 0.8)
ESOutLPlot <- Map(addTitle, ESOutLPlot, dataAlias, 
                  paste("Outlier", paste(outL, collapse = ", "), "at alias-level"))
showMultiPlot(ESOutLPlot)

# Save plot
saveMultiPlot(ESOutLPlot, dataAlias, "outlier-es-alias.png", 
    					width = 1000, height = 1000)


# Remove outliers
outL <- c("GSM867986", "GSM867947", "GSM867948", "GSM867949")
nuke <- lapply(getData(datasets), function(data){
								data[ , !(colnames(data) %in% outL)]
							})
lapply(nuke, dim)

# New name for nuked sets
names(nuke) <- gsub("Clean", "Nuke", datasets)
nuke <- Map(function(data, name){
					assign(name, data, .GlobalEnv)
				}, nuke, names(nuke))
datasetsNuke <- names(nuke)

methylMetaNuke <-
	methylMetaClean[!(methylMetaClean$geo %in% outL), ]
dim(methylMetaNuke)
metaNuke <- methylMetaNuke

# What does it look like after outliers removed
# Heatmap
hDat <- lapply(getData(datasetsNuke), getHeatMat, meta = metaNuke)
hPlot <- lapply(hDat, function(dat, ...){dev.new(); plotCor(dat, ...)}, 
                Rowv = NA, Colv = NA, margin = c(5, 15), 
                labCol = metaNuke$geo)
hPlot <- lapply(hDat, function(dat, ...){dev.new(); plotCor(dat, ...)}, 
                margin = c(5, 15), 
                labCol = metaNuke$geo) # With clustering

# Save Pearson correlation plot
hDat <- lapply(getData(datasetsNuke), getHeatMat, meta = metaNuke)
hPlot <- Map(function(dat, alias){
							png(gsub("alias", alias, "cor-alias-before-norm-no-out.png"),
							    width = 1000, height = 1000)
							plotCor(dat, margin = c(c(5, 15)), labCol = methylMetaClean$geo)
							dev.off()
						 }, 
             hDat, dataAlias) # With clustering

# Beta density
densDat <- lapply(getData(datasetsNuke), getBeta, meta = metaNuke)
densPlot <- lapply(densDat, plotDensity, auto.key = TRUE, 
                   plot.points = FALSE, xlab = "Beta")
densPlot <- Map(addTitle, densPlot, dataAlias,
    						"Density of Beta Values by Cell Type (alias)")
showMultiPlot(densPlot)

# Save beta density plot
saveMultiPlot(densPlot, dataAlias, "beta-density-alias-before-norm-no-out.png", 
    					width = 1000, height = 1000)


# Save no-outlier tables
save(outL, methylMetaNuke, methylDatNuke, 
     file = "450kMethylationData_probeLevel_nuke.RData")
save(avgMethylByGeneNuke, 
     file = "450kMethylationData_geneLevelAverage_nuke.RData")
save(avgMethylByGenePromoterNuke,
     file = "450kMethylationData_geneLevelPromoterAverage_nuke.RData")


