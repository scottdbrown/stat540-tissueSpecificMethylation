####################################################
#### Explore Expression/Methylation Correlation ####
#### Author(s): Scott Brown
#### Date Created: March 31, 2014
#### Last Edited by: Scott Brown
#### on: April 2, 2014
####################################################

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"correlation")
setwd(workingDir)

# load methylation data
#load("../methylation/450kMethylationData_geneLevelPromoterAverage_clean.RData")
load("../methylation/450kMethylationData_geneLevelAverage_clean.RData")

# load expression data
#load("../expression/HT12v3_avgExpressionByGene.RData")
avgExpByGeneClean <- read.table("../expression/expGeneCleanedData.tsv", header=T, sep=" ")

#set variables
methyl <- get(ls()[grepl("avgMethyl",ls())])
exp <- avgExpByGeneClean

shared_genes <- intersect(rownames(exp), rownames(methyl))


expconv <- read.table("../expression/expGSM_gsatid.tsv", header=T, sep="\t")
methylconv <- read.table("../methylation/methylGSM_gsatid.tsv", header=T, sep="\t")

colnames(exp) <- expconv[colnames(exp),"gsatid"]
colnames(methyl) <- methylconv[colnames(methyl),"gsatid"]

samples <- as.character(expconv$gsatid)

exp.shared <- exp[shared_genes,samples]
methyl.shared <- methyl[shared_genes,samples]

# get correlation (genome wide) for each sample between methylation values and expression values.
# cor(log(exp.shared[,3]), methyl.shared[,3], use="complete.obs")

genome.cor <- sapply(samples, function(x){
  cor(log(exp.shared[,x]), methyl.shared[,x], use="complete.obs")
})

ids <- read.table("../IDs.tsv", header=T, sep="\t")

cordat <- data.frame(gsatid=ids$gsatid, group=ids$Group)
cordat$gsatid <- as.character(paste0("g",cordat$gsatid))
cordat$gsatid <- factor(cordat$gsatid, levels=cordat$gsatid, ordered=T)
cordat$cor <- genome.cor[cordat$gsatid]

cordat$group2 <- "Stem"
cordat$group2[grepl("Somatic",cordat$group)] <- "Somatic"

library(ggplot2)
ggplot(cordat, aes(group, cor)) + geom_violin(fill="lightblue") + 
  geom_boxplot(width=.1) + stat_summary(fun.y="mean", geom="point", color="#FF6600", size=3) + 
  geom_point(position=position_jitter(width=.1), alpha=1/2, pch=4) +
  ylab("Correlation") + xlab("Cell type") + ggtitle("Correlation between Expression and Promoter Methylation")


t.test(cordat$cor~cordat$group2)

ggplot(cordat, aes(gsatid, cor, color=group2)) + geom_point(stat="identity", size=3)

ggplot(cordat, aes(gsatid, cor, color=group)) + geom_point(stat="identity", size=3)


summary(lm(cor~group, data=cordat))


## TODO:
## - same analysis for differentially expressed genes.
## - xyplot showing what the data looks like for some samples (genomewide and differential ones)


topExp <- read.table("../expression/expTypeTable.tsv", header=T, sep=" ")
#load("../methylation/450kMethylationData_geneLevelPromoterAverage_hit_clean.RData")
load("../methylation/450kMethylationData_geneLevelAverage_hit_clean.RData")
topMethyl <- avgMethylByGeneCleanHit

#length(diff_shared <- intersect(rownames(topExp[1:1000,]), rownames(topMethyl[1:1000,])))
length(diff_shared <- intersect(rownames(topExp[topExp$adj.P.Val<1e-5,]), rownames(topMethyl[topMethyl$adj.P.Val<1e-5,])))

topExp.shared <- exp[diff_shared,samples]
topMethyl.shared <- methyl[diff_shared,samples]

# get correlation (genome wide) for each sample between methylation values and expression values.
# cor(log(exp.shared[,3]), methyl.shared[,3], use="complete.obs")

top.genome.cor <- sapply(samples, function(x){
  cor(log(topExp.shared[,x]), topMethyl.shared[,x], use="complete.obs")
})

cordat$topcor <- top.genome.cor[cordat$gsatid]


library(ggplot2)
ggplot(cordat, aes(group, cor)) + geom_violin(fill="lightblue") + 
  geom_boxplot(width=.1) + stat_summary(fun.y="mean", geom="point", color="#FF6600", size=3) + 
  geom_point(position=position_jitter(width=.1), alpha=1/2, pch=4) +
  ylab("Correlation") + xlab("Cell type") + ggtitle("Correlation between Expression and Promoter Methylation")


t.test(cordat$topcor~cordat$group2)
t.test(cordat$cor~cordat$group2)

ggplot(cordat, aes(gsatid, topcor, color=group2)) + geom_point(stat="identity", size=3)


ggplot(cordat, aes(group2, topcor)) + geom_boxplot(width=.2) + geom_jitter() + ylim(-.4,0) + geom_hline(aes(yintercept=0))
ggplot(cordat, aes(group2, cor)) + geom_boxplot(width=.2) + geom_jitter() + ylim(-.4,0) + geom_hline(aes(yintercept=0))


ggplot(cordat, aes(gsatid, cor, color=group)) + geom_point(stat="identity", size=3)


summary(lm(cor~group, data=cordat))
