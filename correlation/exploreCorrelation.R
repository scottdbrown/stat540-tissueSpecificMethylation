####################################################
#### Explore Expression/Methylation Correlation ####
#### Author(s): Scott Brown
#### Date Created: March 31, 2014
#### Last Edited by: Scott Brown
#### on: March 31, 2014
####################################################

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"correlation")
setwd(workingDir)

# load methylation data
load("../methylation/450kMethylationData_geneLevelPromoterAverage_clean.RData") #change to promoter?

# load expression data
#load("../expression/HT12v3_avgExpressionByGene.RData")
avgExpByGeneClean <- read.table("../expression/expGeneCleanedData.tsv", header=T, sep=" ")

shared_genes <- intersect(rownames(avgExpByGeneClean), rownames(avgMethylByGenePromoterClean))


expconv <- read.table("../expression/expGSM_gsatid.tsv", header=T, sep="\t")
methylconv <- read.table("../methylation/methylGSM_gsatid.tsv", header=T, sep="\t")

colnames(avgExpByGeneClean) <- expconv[colnames(avgExpByGeneClean),"gsatid"]
colnames(avgMethylByGenePromoterClean) <- methylconv[colnames(avgMethylByGenePromoterClean),"gsatid"]

samples <- as.character(expconv$gsatid)

avgExpByGeneClean.shared <- avgExpByGeneClean[shared_genes,samples]
avgMethylByGenePromoterClean.shared <- avgMethylByGenePromoterClean[shared_genes,samples]

# get correlation (genome wide) for each sample between methylation values and expression values.
cor(log(avgExpByGeneClean.shared[,3]), avgMethylByGenePromoterClean.shared[,3], use="complete.obs")

genome.cor <- sapply(samples, function(x){
  cor(log(avgExpByGeneClean.shared[,x]), avgMethylByGenePromoterClean.shared[,x], use="complete.obs")
})

ids <- read.table("../IDs.tsv", header=T, sep="\t")

cordat <- data.frame(gsatid=ids$gsatid, group=ids$Group)
cordat$gsatid <- as.character(paste0("g",cordat$gsatid))
cordat$gsatid <- factor(cordat$gsatid, levels=cordat$gsatid, ordered=T)
cordat$cor <- genome.cor[cordat$gsatid]

cordat$group2 <- "Stem"
cordat$group2[grepl("Somatic",cordat$group)] <- "Somatic"

library(ggplot2)
ggplot(cordat, aes(group2, cor)) + geom_violin(fill="lightblue") + 
  geom_boxplot(width=.1) + stat_summary(fun.y="mean", geom="point", color="#FF6600", size=3) + 
  ylab("Correlation") + xlab("Cell type") + ggtitle("Correlation between Expression and Promoter Methylation")


t.test(cordat$cor~cordat$group2)

ggplot(cordat, aes(gsatid, cor, color=group2)) + geom_point(stat="identity", size=3)

ggplot(cordat, aes(gsatid, cor, color=group)) + geom_point(stat="identity", size=3)

