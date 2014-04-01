#### Make file for ID conversion ####
#### Scott Brown
#### April 1, 2014

projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- projRoot
setwd(workingDir)

## Expression:

load("/expression/HT12v3_avgExpressionByGene.RData")

#Read in metadata title IDs
IDs <- read.table("IDs.tsv", header=T, sep="\t")
expMeta$title <- as.character(expMeta$title)
expMeta$title <- as.factor(expMeta$title)
expMeta <- merge(expMeta, IDs[,c("gsatid","HT12v3_ID")], by.x="title", by.y="HT12v3_ID")


expGSM_gsatid <- expMeta[order(expMeta$gsatid),c("geo_accession","gsatid")]
rownames(expGSM_gsatid) <- expGSM_gsatid$geo_accession
write.table(expGSM_gsatid, "expression/expGSM_gsatid.tsv", sep="\t")


## Methylation:

load("methylation/450kMethylationData_probeLevel.RData")  #to get methylMeta

#Read in metadata title IDs
methylMeta$title <- as.character(methylMeta$title)
methylMeta$title <- as.factor(methylMeta$title)
methylMeta <- merge(methylMeta, IDs[,c("gsatid","X450k_ID")], by.x="title", by.y="X450k_ID")


methylGSM_gsatid <- methylMeta[order(methylMeta$gsatid),c("geo_accession","gsatid")]
rownames(methylGSM_gsatid) <- methylGSM_gsatid$geo_accession
write.table(methylGSM_gsatid, "methylation/methylGSM_gsatid.tsv", sep="\t")
