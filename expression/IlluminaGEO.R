# ===================================================
# Illumina HT12v3 GEO-querying Subroutine (Nathaniel)
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
Illumina.Table.Repack <- data.frame(ID = Illumina.GEO.TRIM$ID, HUGO = Illumina.GEO.TRIM$ILMN_Gene)

# RData Repacking
dput(Illumina.Table.Repack, file = 'IlluminaHTMappingTable.HUGO')
