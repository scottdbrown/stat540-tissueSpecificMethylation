#-------------------------------------------------------------------------------
# Differential methylation analysis
# With outliers included
# Jessica Lee
# Date created: April 2, 2014
# Last edit: April 3, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

# Load libraries
library(wateRmelon)
library(limma)

# Read in data
load("450kMethylationData_probeLevel_clean_norm.RData")
load("450kMethylationData_geneLevelAverage_clean_norm.RData")
load("450kMethylationData_geneLevelPromoterAverage_clean_norm.RData")
load("450kMethylationData_meta_clean.RData")

# Get me variable names
names <- ls()
datasets <- names[grep("Meta", names, invert = TRUE)] # data only
betasets <- names[grep("Beta", names)]
datasets <- names[grep("Beta", datasets, invert = TRUE)]
meta <- get(names[grep("Meta", names)]) # meta only

# Dataset alias
dataAlias <- c("gene", "promoter", "probe")
names(dataAlias) <- datasets

# Load helper functions
source("helpers.R")

#-------------------------------
# Functions
#-------------------------------
fitAndHit <- function(data, design, coefName, ...) {
	dmFit <- lmFit(data, design)
	dmEbFit <- eBayes(dmFit)
	topTable(dmEbFit, 
	         coef = grep(coefName, colnames(coef(dmEbFit))), 
	         ...)
}

fitAndHitAll <- function(data, formula, meta, coefName, ...) {
	desMat <- model.matrix(formula, meta)	
	lapply(data, fitAndHit, design = desMat, 
	       coefName = coefName, ...)
}

#-------------------------------
# Workspace
#-------------------------------
# Limma throws an error for some probes that have 
# normalized values of -Inf or Inf
# Nuke probes that have -Inf or Inf
nonInfRows <- apply(methylDatCleanNorm, 1, function(row){ 
	all((row != -Inf) & (row != Inf)) 
})
length(which(nonInfRows == FALSE))
methylDatCleanRmInf <- methylDatCleanNorm[nonInfRows, ]
dim(methylDatCleanRmInf)
dim(methylDatCleanNorm)

# Save RmInf as new data so we can lapply in next step
newData <- gsub("methylDatCleanNorm", "methylDatCleanRmInf", datasets)

# Differential methylation analysis with limma
# Get design matrix
desMat <- model.matrix(~ cellTypeSimple, meta)

# Fit linear model, test that cell types doesn't matter
# Get hits
hit <- lapply(getData(newData), fitAndHit, design = desMat, 
             coefName = "cellType", n = Inf)

# One shot method to test various models
hit <- fitAndHitAll(getData(newData), 
                    ~ cellTypeSimple, meta,
                    coefName = "cellType", n = Inf)

# Save top hits
hit <- castGlobal(hit, "Norm|RmInf", "Hit")

# Save data frame
save(methylDatCleanHit,
     file = "450kMethylationData_probeLevel_hit_clean.RData")
save(avgMethylByGeneCleanHit,
     file = "450kMethylationData_geneLevelAverage_hit_clean.RData")
save(avgMethylByGenePromoterCleanHit,
     file = "450kMethylationData_geneLevelPromoterAverage_hit_clean.RData")


#-------------------------------
# Plotting work
#-------------------------------
# Visualize top 100 hits in heatmap
top100 <- lapply(hit, function(x){head(x, n = 100)})
hitBeta <- Map(function(top, beta){
	subset(beta, rownames(beta) %in% rownames(top))
}, top100, getData(betasets))

getHeatmap <- function(data, meta) {
	dev.new()
	col <- c("darkgoldenrod1", "forestgreen")
	heatmap.2(data, col = rdBu, ColSideColors = col[meta$cellTypeSimple], density.info = "none", 
	    trace = "none", Rowv = TRUE, Colv = TRUE, labCol = FALSE, 
	    dendrogram = "row", margins = c(1, 15))
	legend("topright", c("stem cell", "somatic cell"), 
	       col = c("darkgoldenrod1", "forestgreen"), 
	    	 pch = 15)	
}
lapply(hitBeta, getHeatmap, meta)

# Visualize as Stripplot
top5 <- lapply(hit, function(x){head(x, n = 5)})
hitBeta <- Map(function(top, beta){
	subset(beta, rownames(beta) %in% rownames(top))
}, top5, getData(betasets))

getStrip <- function(data, meta) {
	tmp <- melt(data)
	colnames(tmp) <- c("gene", "sample", "beta")
	tmp <- cbind(tmp, cellType = meta$cellTypeSimple[tmp$sample])
	
	return (ggplot(tmp, aes(cellType, beta, color = cellType)) +
		geom_point(position = position_jitter(width = 0.05)) +
		stat_summary(fun.y = mean, aes(group = 1), geom = "line", color = "black") + 
		facet_grid(~ gene) + xlab("Cell Type") + ylab("Beta Value") +
		labs(title = "Top 5 DMR") + 
		theme(legend.title=element_blank()))
}
top5Strip <- lapply(hitBeta, getStrip, meta)
top5Strip <- showMultiPlot(top5Strip)

