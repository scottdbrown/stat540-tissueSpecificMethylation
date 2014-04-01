#-------------------------------------------------------------------------------
# Some repeatedly used functions
# Saved separately for easier reference
#-------------------------------------------------------------------------------
# Get datasets specified in names vector
getData <- function(names) {
	getDataByName <- function(name) {
		return (get(name))
	}
	data <- lapply(names, getDataByName)
	names(data) <- names
	return (data)
}

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

# Save a heatmap to png format
saveHeatmap <- function(data, ...) {
	png(...)
	heatmap(data, ...)
	dev.off()
}

# Save multiple heatmap
saveMultiHeatmap <- function(data, alias, fileName, ...) {
	args <- c(list(data, gsub("alias", alias, fileName)), ...)
	do.call(savePlot, args)
}

