# =====================================
# RNA-Seq Preprocessed Code (Nathaniel)
# =====================================

# Session Cleanup
# rm(list = ls())

# Loading Libraries

# Loading RData Container
load('HT12v3_expression.RData')

# Obtaining Illumina HT12v3 Annotation from File
Illumina.HUGO <- dget('IlluminaHTMappingTable.HUGO')
rownames(Illumina.HUGO) <- Illumina.HUGO$ID

# Obtaining Illumina Mappings for Illumina HT12v3 Platform (HUGO SYMBOL)
annot.expDat <- data.frame(expDat, HUGO = Illumina.HUGO[rownames(expDat),]$HUGO)
