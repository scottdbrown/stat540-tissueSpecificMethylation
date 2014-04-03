#-------------------------------------------------------------------------------
# Explore methylation data
# Jessica Lee
# Date created: March 24, 2014
# Last edit: April 1, 2014
#-------------------------------------------------------------------------------
# Set working directories
setwd("~/workspace/stat540.proj/methylation")

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
# Remove outliers
outL <- c("GSM867986", "GSM867947", "GSM867948", "GSM867949")
nuke <- lapply(getData(datasets), function(data){
                data[ , !(colnames(data) %in% outL)]
              })
lapply(nuke, dim)

# New name for nuked sets
nuke <- castGlobal(nuke, "Clean", "Nuke")
datasetsNuke <- names(nuke)

methylMetaNuke <-
  methylMetaClean[!(methylMetaClean$geo %in% outL), ]
dim(methylMetaNuke)
metaNuke <- methylMetaNuke


# Save no-outlier tables
save(outL,
     file = "450kMethylationData_outlier.RData")
save(methylMetaNuke, 
     file = "450kMethylationData_meta_nuke.RData")
save(methylDatNuke,  
     file = "450kMethylationData_probeLevel_nuke.RData")
save(avgMethylByGeneNuke, 
     file = "450kMethylationData_geneLevelAverage_nuke.RData")
save(avgMethylByGenePromoterNuke,
     file = "450kMethylationData_geneLevelPromoterAverage_nuke.RData")


#-------------------------------
# Plotting work
#-------------------------------
# Density of Beta values
densDat <- lapply(getData(datasets), getBetaAvg, meta = meta)
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
densDat <- lapply(getData(datasetsNuke), getBetaAvg, meta = metaNuke)
densPlot <- lapply(densDat, plotDensity, auto.key = TRUE, 
                   plot.points = FALSE, xlab = "Beta")
densPlot <- Map(addTitle, densPlot, dataAlias,
    						"Density of Beta Values by Cell Type (alias)")
showMultiPlot(densPlot)

# Save beta density plot
saveMultiPlot(densPlot, dataAlias, "beta-density-alias-before-norm-no-out.png", 
    					width = 1000, height = 1000)


