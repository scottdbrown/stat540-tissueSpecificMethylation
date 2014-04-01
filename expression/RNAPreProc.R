# =====================================
# RNA-Array Preprocessed Code (Nathaniel)
# 31st March 2014
# =====================================

# Session Cleanup
# rm(list = ls())

# Loading Libraries

# Loading RData Container
load('HT12v3_expression.RData')

# Obtaining Illumina HT12v3 Annotation from File
load('IlluminaHTMappingTable.HUGO')

# Obtaining Illumina Mappings for Illumina HT12v3 Platform (HUGO SYMBOL)
annot.expDat <- data.frame(expDat, HUGO = Illumina.HUGO[rownames(expDat),]$HUGO)
