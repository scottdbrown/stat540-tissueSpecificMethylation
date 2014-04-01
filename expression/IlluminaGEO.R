# ===================================================
# Illumina HT12v3 GEO-querying Subroutine (Nathaniel)
# 31st March 2014
# ===================================================

# Loading Library
library(GEOquery)

# Obtaining Illumina HT12v3 Annotation
# Illumina.GEO <- getGEO('GPL6947')

# Alternative Calling
# Obtaining Illumina HT12v3 Annotation from File
Illumina.GEO <- getGEO(filename = 'GPL6947.soft')

# Annotation Dataframe Restructuring (for HUGO SYMBOLs only)
Illumina.GEO.TRIM <- Illumina.GEO@dataTable@table
Illumina.HUGO <- data.frame(ID = Illumina.GEO.TRIM$ID, HUGO = Illumina.GEO.TRIM$ILMN_Gene)
rownames(Illumina.HUGO) <- Illumina.HUGO$ID

# RData Repacking
save(Illumina.HUGO, file = 'IlluminaHTMappingTable.HUGO')