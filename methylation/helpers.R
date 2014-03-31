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

# Save a given plot to png format
savePlot <- function(plot, ...) {
	png(...)
	plot(plot)
	dev.off()
}
