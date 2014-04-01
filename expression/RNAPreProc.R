# =====================================
# RNA-Array Preprocessed Code (Nathaniel)
# 31st March 2014
# Last edit: Scott Brown - March 31, 2014
# =====================================

# Session Cleanup
# rm(list = ls())

# Loading Libraries

# Loading RData Container
load('HT12v3_expression.RData')

# Obtaining Illumina HT12v3 Annotation from File
load('IlluminaHTMappingTable.HUGO')

# Obtaining Illumina Mappings for Illumina HT12v3 Platform (HUGO SYMBOL)
annot.expDat <- data.frame(HUGO = Illumina.HUGO[rownames(expDat),]$HUGO, expDat)

# Get average expression of each gene (collapse data to one value per gene)
avgExpByGene <- apply(annot.expDat[,-1],2,by, annot.expDat$HUGO, mean, na.rm=T)
save(avgExpByGene, expMeta, file="HT12v3_avgExpressionByGene.RData")
