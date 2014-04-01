# ===================================================
# Illumina HT12v3 GEO-querying Subroutine (Nathaniel)
# 1st April 2014
# ===================================================

# Loading Library
library(GEOquery)

# Obtaining Illumina HT12v3 Annotation
# Illumina.GEO <- getGEO('GPL6947')

# Alternative Calling
# Obtaining Illumina HT12v3 Annotation from File
Illumina.GEO <- getGEO(filename = 'GPL6947.soft')

# Extracting DataTable
Illumina.GEO.TRIM <- Illumina.GEO@dataTable@table

# Annotation Dataframe Restructuring (for HUGO SYMBOLs only)
Illumina.HUGO <- data.frame(ID = Illumina.GEO.TRIM$ID, HUGO = Illumina.GEO.TRIM$ILMN_Gene)
rownames(Illumina.HUGO) <- Illumina.HUGO$ID

# Annotation Dataframe Restructuring (for HUGO + ILMN_Trascript ID only) *Experimental*
Illumina.EXPR <- data.frame(ID = Illumina.GEO.TRIM$ID, HUGO = Illumina.GEO.TRIM$ILMN_Gene, TranscriptID = Illumina.GEO.TRIM$Transcript)
rownames(Illumina.EXPR) <- Illumina.EXPR$ID

# RData Repacking
save(Illumina.HUGO, file = 'IlluminaHTMappingTable.HUGO')
save(Illumina.EXPR, file = 'IlluminaHTMappingTable.EXPR')