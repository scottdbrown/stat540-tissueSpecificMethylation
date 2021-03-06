# Deep Analysis Script
# by Nathaniel (7 April)

# Loads Final Result Package
load('Analysis.RDAT')

# Create subset Tables using FDR < 1e-5
tinyTable <- topTableType[topTableType$adj.P.Val < 1e-5, ]
tinyPromoter <- avgMethylByGenePromoterCleanHit[avgMethylByGenePromoterCleanHit$adj.P.Val < 1e-5,]
tinyGene <- avgMethylByGeneCleanHit[avgMethylByGeneCleanHit$adj.P.Val < 1e-5,]

# Differential Methylation
# Union
methylation.union <- union(rownames(tinyGene), rownames(tinyPromoter))
write.table(methylation.union, file = 'methylation.all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Individual
write.table(rownames(tinyGene), file = 'methylation.gene.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(rownames(tinyPromoter), file = 'methylation.promoter.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Differential Expression
write.table(rownames(tinyTable), file = 'expression.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Methylation Union against Expression
all.methylation.intersect.expression <- intersect(rownames(tinyTable), union(rownames(tinyGene), rownames(tinyPromoter)))
write.table(all.methylation.intersect.expression, file = 'all.methylation.intersect.expression.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Triple Intersect
triple.intersect <- intersect(intersect(rownames(tinyGene), rownames(tinyPromoter)), rownames(tinyTable))
write.table(triple.intersect, file = 'triple.intersect.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Promoter Methylation Intersect Expression
prom.intersect <- intersect(rownames(tinyPromoter), rownames(tinyTable))
write.table(prom.intersect, file = 'promoter.intersect.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Gene Methylation Intersect Expression
gene.intersect <- intersect(rownames(tinyGene), rownames(tinyTable))
write.table(gene.intersect, file = 'gene.intersect.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)