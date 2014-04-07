# Deep Analysis Script
# by Nathaniel (7 April)

# Loads Final Result Package
load('Analysis.RDAT')

# Create subset Tables using FDR < 1e-5
tinyTable <- topTableType[topTableType$adj.P.Val < 1e-5, ]
tinyPromoter <- avgMethylByGenePromoterCleanHit[avgMethylByGenePromoterCleanHit$adj.P.Val<1e-5,]
tinyGene <- avgMethylByGeneCleanHit[avgMethylByGeneCleanHit$adj.P.Val<1e-5,]

# Differential Methylation
methylation.union <- union(rownames(tinyGene), rownames(tinyPromoter))
write.table(methylation.union, file = 'methylation.all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

triple.intersect <- intersect(intersect(rownames(tinyGene), rownames(tinyPromoter)), rownames(tinyTable))
triple.intersect.table <- data.frame(ID = triple.intersect, gMethyl = tinyGene[triple.intersect,]$adj.P.Val, pMethyl = tinyPromoter[triple.intersect,]$adj.P.Val, gExp = tinyTable[triple.intersect,]$adj.P.Val)

# Writing HUGO files for downstream GO-Analysis
write.table(triple.intersect.table.names, file = '3IHUGO.txt', quote = FALSE,col.names = FALSE, row.names = FALSE)
