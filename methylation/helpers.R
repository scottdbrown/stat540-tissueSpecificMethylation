#-------------------------------------------------------------------------------
# Some repeatedly used functions
# Some often used libraries
# Saved here separately for easier reference
#-------------------------------------------------------------------------------
# Load libraries
library(RColorBrewer)
library(reshape2)
library(gplots)
library(hexbin)
library(lattice)

# Load color scheme
rdBu <- colorRampPalette(brewer.pal(n = 11, "RdBu"))

# Get datasets specified in names vector
getData <- function(names) {
	getDataByName <- function(name) {
		return (get(name))
	}
	data <- lapply(names, getDataByName)
	names(data) <- names
	return (data)
}

#-------------------------------
# PLOT HELPERS
#-------------------------------
# Add title to a plot where place holder "alias" is (the word "alias" will be 
# replaced by actual alias); then return the plot
addTitle <- function(plot, alias, title) {
	plot$main <- gsub("alias", alias, title)
	return(plot)
}

# Show multiple plots on screen
showMultiPlot <- function(plotList, ...) {
	# Plot displaying function
	show <- function(plot){
		dev.new()
		plot(plot, ...)
	}
	lapply(plotList, show)
}

# Show multiple plots on screen
showMultiHeatmap <- function(data, ...) {
	show <- function(plot){
		dev.new()
		heatmap(plot, ...)
	}
	lapply(data, show)
}

# Save a given plot to png format
savePlot <- function(plot, ...) {
	png(...)
	plot(plot)
	dev.off()
}

# Save multiple plot
saveMultiPlot <- function(plots, aliases, fileName, ...) {
	callSave <- function(plot, alias, fileName, ...) {
		args <- c(list(plot, gsub("alias", alias, fileName)), ...)
		do.call(savePlot, args)
	}

	Map(callSave, plots, aliases, fileName, ...)
}

#-------------------------------
# EXPLORE HELPERS
#-------------------------------
# Get beta values
getBetaAvg <- function(data, meta) {
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

# Rename dataset and assign it at global scope
castGlobal <- function(datasets, oldVar, newVar){
	names(datasets) <- gsub(oldVar, newVar, names(datasets))
	Map(function(data, name){
		assign(name, data, .GlobalEnv)
	}, datasets, names(datasets))
}
