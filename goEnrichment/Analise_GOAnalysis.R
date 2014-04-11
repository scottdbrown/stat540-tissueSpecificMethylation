Exploration of Embryonic Stem Cell specific genes for GO Term enrichment
========================================================

topTableType  <- read.table(file = "expression//expTypeTable.tsv")
load("methylation//450kMethylationData_geneLevelPromoterAverage_hit_clean.RData")

#make FDR cut off subset
expSubset  <- subset(topTableType, adj.P.Val < 1e-5)
methPromoterSubset  <- subset(avgMethylByGenePromoterHit, adj.P.Val < 1e-5)

#get intersect of Methylation based on Methylation by promoter with expression
MethPromoterExp.Intersect  <- intersect(row.names(methPromoterSubset), row.names(expSubset))
write.table(MethPromoterExp.Intersect, file = "goEnrichment/methpro.exp.Intersect.txt", row.names = F, col.names = F)

#try to find some biological meaning by looking at some Stem cell specific genes and looking where they fall in the intersects
hESCGenes  <- read.table(file = "goEnrichment/hESC_expressedGenes.txt")
hESCGenes.Hits.Intersect  <- intersect(MethPromoterExp.Intersect, hESCGenes[,"V1"])
#hESCGenes.Hits.Intersect  <- as.character(hESCGenes.Hits.Intersect)
write.table(hESCGenes.Hits.Intersect, file = "goEnrichment/methpro.exp.hESC.Intersect.txt", row.names = F, col.names = F)

#make a heatmap of the hESC specific genes that overlap with MethPromoter + Expression
heatmap(as.matrix(log.avgExpCleanNorm[c(hESCGenes.Hits.Intersect),]), scale = "none", Rowv = NA, main = "hESC Specific Genes in Intersecting Genes Promoter")

#Need to get a copy of the cleaned and normalized promoter/avg gene meth data
load(file = "methylation//450kMethylationData_geneLevelAverage.RData")
load(file = "methylation//450kMethylationData_geneLevelPromoterAverage.RData")
load(file = "450kMethylationData_geneLevelAverage_clean.RData") #loaded top table of Methylation by gene average (avgMethylByGeneCleanHit)
load(file = "methylation//450kMethylationData_geneLevelAverage.RData")

methAvgGene  <- subset(avgMethylByGeneCleanHit, adj.P.Val <1e-5)

intersect.methAvg.Exp  <- intersect(row.names(methAvgGene), row.names(expSubset))
write.table(intersect.methAvg.Exp, file = "goEnrichment/methAvg.Exp.Intersect.txt", row.names = F, col.names = F)
intersect.methAvg.hESC  <- intersect(intersect.methAvg.Exp, hESCGenes[,"V1"])
write.table(intersect.methAvg.hESC, file = "goEnrichment/methAvg.Exp.hESC.Intersect.txt", row.names = F, col.names = F)
heatmap(as.matrix(log.avgExpCleanNorm[c(intersect.methAvg.hESC),]), scale = "none", Rowv = NA, main = "hESC Specific Genes & Methylation Gene Average")

intersect.exp.hESC  <- intersect(row.names(expSubset), hESCGenes[,"V1"]) 
write.table(intersect.exp.hESC, file = "goEnrichment/exp.hESC.Intersect.txt", row.names = F, col.names = F)
heatmap(as.matrix(log.avgExpCleanNorm[c(intersect.exp.hESC),]), scale = "none", main = "Expression of hESC and Expression Gene Overlap", Rowv = NA)

intersect.methAvg.hESC  <- intersect(row.names(methAvgGene), hESCGenes[, "V1"])
write.table(intersect.methAvg.hESC, file = "goEnrichment/methPromoter.hESC.Intersect.txt", row.names = F, col.names = F)

intersect.methPromoter.hESC  <- intersect(row.names(methPromoterSubset), hESCGenes[,"V1"])
write.table(intersect.methPromoter.hESC, file = "goEnrichment/methPromoter.hESC.Intersect.txt", row.names = F, col.names = F)

intersect.all.meth.exp  <- intersect(intersect.exp.hESC, intersect.methPromoter.hESC)
intersect.all.meth.exp  <- intersect(intersect.all.meth.exp, intersect.methAvg.hESC)
heatmap(as.matrix(log.avgExpCleanNorm[c(intersect.all.meth.exp),]), scale = "none")

write.table(intersect.all.meth.exp, file = "goEnrichment/methAvg.methPro.Exp.hESC.Intersect.txt", row.names = F, col.names = F)

#want the hESC genes that don't overlap with meth promoter + exp
diff.methprot.exp.hESC  <- setdiff(hESCGenes[,"V1"], hESCGenes.Hits.Intersect)
write.table(diff.methprot.exp.hESC, file = "goEnrichment/diff.hESC.overlap.methpro.exp.txt", col.names= F, row.names = F)

