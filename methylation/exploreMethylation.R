#-------------------------------------------------------------------------------
# Explore methylation data
# Jessica Lee
# Date created: March 24, 2014
# Last edit: March 30, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(RColorBrewer)
library(gplots)
library(hexbin)
library(lattice)

# Read in data
load("450kMethylationData_probeLevel_clean.RData")
load("450kMethylationData_geneLevelAverage_clean.RData")

# Get me variable names
names <- ls()

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

# Save a plot
savePlot <- function(plot, ...) {
	png(...)
	plot(plot)
	dev.off()
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
# Density of Beta values - gene level
dGene <- getBeta(avgMethylByGeneClean, methylMetaClean)
dGenePlot <- plotdity(dGene, auto.key = TRUE, 
         							plot.points = FALSE,
          						main = "Density of Beta Values by Cell Type (Gene Level)", 
          						xlab = "Beta")
dGenePlot

# Save beta density plot
savePlot(dGenePlot, "beta-density-before-norm-gene.png", 
         width = 1000, height = 1000)

# Density of Beta values - probe level
dProbe <- getBeta(methylDatClean, methylMetaClean)
dProbePlot <- plotdity(dProbe, auto.key = TRUE, 
         							 plot.points = FALSE,
          						 main = "Density of Beta Values by Cell Type (Probe Level)", 
          						 xlab = "Beta")
dProbePlot

# Save beta density plot
savePlot(dProbePlot, "beta-density-before-norm-probe.png", 
         width = 1000, height = 1000)


# Gene-level pearson correlation
hGene <- getHeatMat(avgMethylByGeneClean, methylMetaClean)
hGenePlot <- plotCor(hGene, Rowv = NA, Colv = NA, margin = c(4, 10), 
                     labCol = methylMetaClean$geo)
hGenePlot <- plotCor(hGene, margin = c(4, 10), 
                     labCol = methylMetaClean$geo) # With clustering

# Probe-level pearson correlation
hProbe <- getHeatMat(methylDatClean, methylMetaClean)
hProbePlot <- plotCor(hProbe, Rowv = NA, Colv = NA, margin = c(4, 10), 
                      labCol = methylMetaClean$geo)
hProbePlot <- plotCor(hProbe, margin = c(4, 10), 
                      labCol = methylMetaClean$geo) # With clustering
# >> looks like GSM867986, GSM867947, GSM867948, GSM867949 are scary batches


# Examine possible outlier GSM867986 by plotting with its replicate GSM867987
outL <- c("GSM867986")

# Examine the outlier in a cellType context
# iPS outlier
outLGenePlot <- examOutL(avgMethylByGeneClean, outL, 
         								 unlist(subset(methylMetaClean, 
         								        			 cellTypeShort == "iPS" & 
                											 !(geo %in% outL), select = geo)),
         								 main = "Outlier at gene-level")
outLGenePlot

outLProbePlot <- examOutL(methylDatClean, outL, 
         									unlist(subset(methylMetaClean, 
         									       				cellTypeShort == "iPS" & 
                								 				!(geo %in% outL), select = geo)),
         									main = "Outlier at probe-level")
outLProbePlot

# ES outliers
outL <- c("GSM867947", "GSM867948", "GSM867949")
outLGenePlot <- examOutL(avgMethylByGeneClean, outL, 
         								 unlist(subset(methylMetaClean, 
         								        			 cellTypeShort == "ES" & 
                											 !(geo %in% outL), select = geo)),
         								 main = "Outlier at gene-level")
outLGenePlot

outLProbePlot <- examOutL(methylDatClean, outL, 
         								 unlist(subset(methylMetaClean, 
         								        			 cellTypeShort == "ES" & 
                											 !(geo %in% outL), select = geo)),
         								 main = "Outlier at probe-level")
outLProbePlot

bwplot(geneExp ~ devStage, oDat,
       panel = panel.violin)


# Remove outliers
outL <- c("GSM867986", "GSM867947", "GSM867948", "GSM867949")
avgMethylByGeneNuke <- 
	avgMethylByGeneClean[ , !(colnames(avgMethylByGeneClean) %in% outL)]
dim(avgMethylByGeneNuke)

methylDatNuke <- 
	methylDatClean[ , !(colnames(methylDatClean) %in% outL)]
dim(methylDatNuke)

methylMetaNuke <-
	methylMetaClean[!(methylMetaClean$geo %in% outL), ]
dim(methylMetaNuke)

# What does it look like after outliers removed
# Heat map
hGene <- getHeatMat(avgMethylByGeneNuke, methylMetaNuke)
hGenePlot <- plotCor(hGene, margin = c(4, 10), 
                     labCol = methylMetaNuke$geo) # With clustering
hProbe <- getHeatMat(methylDatNuke, methylMetaNuke)
hProbePlot <- plotCor(hProbe, margin = c(4, 10), 
                      labCol = methylMetaNuke$geo) # With clustering

# Beta density
dGene <- getBeta(avgMethylByGeneNuke, methylMetaNuke)
dGenePlot <- plotdity(dGene, auto.key = TRUE, 
         							plot.points = FALSE,
          						main = "Density of Beta Values by Cell Type (Gene Level)", 
          						xlab = "Beta")
dGenePlot

dProbe <- getBeta(methylDatNuke, methylMetaNuke)
dProbePlot <- plotdity(dProbe, auto.key = TRUE, 
         							 plot.points = FALSE,
          						 main = "Density of Beta Values by Cell Type (Probe Level)", 
          						 xlab = "Beta")
dProbePlot


# Save no-outlier tables
save(outL, methylMetaNuke, methylDatNuke, 
     file = "450kMethylationData_probeLevel_clean_nuke.RData")
save(avgMethylByGeneNuke, 
     file = "450kMethylationData_geneLevelAverage_nuke.RData")

