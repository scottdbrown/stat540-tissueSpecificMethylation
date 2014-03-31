####################################
#### Assign Methylation to Gene ####
#### Author(s): Scott Brown
#### Date Created: March 21, 2014
#### Last Edited by: Scott Brown
#### on: March 31, 2014
####################################

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"methylation")
setwd(workingDir)

library(IlluminaHumanMethylation450k.db)

# Load in the probe level methylation data
load("450kMethylationData_probeLevel.RData")

# get probe chromosome, coordinates, and genes
coor <- IlluminaHumanMethylation450kCPGCOORDINATE
mapped_probes <- mappedkeys(coor)
coor <- as.data.frame(coor[mapped_probes])
CHR <- as.data.frame(IlluminaHumanMethylation450kCHR37)
CG_gene <- as.data.frame(IlluminaHumanMethylation450kSYMBOL)

GenoCoor<-merge(merge(CHR, coor,by="Probe_ID"), CG_gene, by.x="Probe_ID", by.y="probe_id")

prom <- as.data.frame(IlluminaHumanMethylation450kREGULATORYGROUP)
GenoCoor<-merge(GenoCoor, prom, by.x="Probe_ID", by.y="regview.Probe_ID")
names(GenoCoor)[5] <- "UCSC_group"

# Merge probe info with data
methylDat <- merge(GenoCoor, cbind(methylDat, probe=row.names(methylDat)), by.x="Probe_ID", by.y="probe", sort=F, all=T)

# Remove unneccesary variables
rm(mapped_probes, CHR, coor, CG_gene, prom)
gc()

# Remove probes not associated with a gene
methylDat <- methylDat[!is.na(methylDat$symbol),]

#### Method 1 ####
# Find average methylation value for each gene
# don't do the mean on the first 5 columns (probe, chromosome, coordinate, gene, promoter)
# Note: takes quite a while to run (~1 hour?)
avgMethylByGene <- apply(methylDat[,-(1:4)],2,by, methylDat$symbol, mean, na.rm=T)
save(avgMethylByGene, file="450kMethylationData_geneLevelAverage.RData")


#### Method 2 ####
# Find average methylation value for each gene's promoter
# don't do the mean on the first 5 columns (probe, chromosome, coordinate, gene, promoter)
# Note: takes quite a while to run (~1 hour?)
methylDatProm <- methylDat[grep("Promoter",methylDat$UCSC_group),]
avgMethylByGenePromoter <- apply(methylDatProm[,-(1:5)],2,by, methylDatProm$symbol, mean, na.rm=T)
save(avgMethylByGenePromoter, file="450kMethylationData_geneLevelPromoterAverage.RData")
