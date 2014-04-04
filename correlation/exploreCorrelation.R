####################################################
#### Explore Expression/Methylation Correlation ####
#### Author(s): Scott Brown
#### Date Created: March 31, 2014
#### Last Edited by: Scott Brown
#### on: April 2, 2014
####################################################

library(ggplot2)

# edit project root folder if needed:
projRoot <- "~/stat540-tissueSpecificMethylation/"

workingDir <- paste0(projRoot,"correlation")
setwd(workingDir)

# load methylation data
load("../methylation/450kMethylationData_geneLevelPromoterAverage_clean.RData")
load("../methylation/450kMethylationData_geneLevelAverage_clean.RData")

# load expression data
avgExpByGeneClean <- read.table("../expression/expGeneCleanedData.tsv", header=T, sep=" ")

#set variables
methylPromoter <- avgMethylByGenePromoterClean
methylGene <- avgMethylByGeneClean
exp <- avgExpByGeneClean

shared_genesG <- intersect(rownames(exp), rownames(methylGene))
shared_genesP <- intersect(rownames(exp), rownames(methylPromoter))


expconv <- read.table("../expression/expGSM_gsatid.tsv", header=T, sep="\t")
methylconv <- read.table("../methylation/methylGSM_gsatid.tsv", header=T, sep="\t")

colnames(exp) <- expconv[colnames(exp),"gsatid"]
colnames(methylPromoter) <- methylconv[colnames(methylPromoter),"gsatid"]
colnames(methylGene) <- methylconv[colnames(methylGene),"gsatid"]

samples <- as.character(expconv$gsatid)

expP.shared <- exp[shared_genesP,samples]
expG.shared <- exp[shared_genesG,samples]
methylP.shared <- methylPromoter[shared_genesP,samples]
methylG.shared <- methylGene[shared_genesG,samples]

# get correlation (genome wide) for each sample between methylation values and expression values.
# cor(log(exp.shared[,3]), methyl.shared[,3], use="complete.obs")

genome.corP <- sapply(samples, function(x){
  cor(log(expP.shared[,x]), methylP.shared[,x], use="complete.obs")
})

genome.corG <- sapply(samples, function(x){
  cor(log(expG.shared[,x]), methylG.shared[,x], use="complete.obs")
})

ids <- read.table("../IDs.tsv", header=T, sep="\t")

cordat <- data.frame(gsatid=ids$gsatid, group=ids$Group)
cordat$gsatid <- as.character(paste0("g",cordat$gsatid))
cordat$gsatid <- factor(cordat$gsatid, levels=cordat$gsatid, ordered=T)
cordat$corP <- genome.corP[cordat$gsatid]
cordat$corG <- genome.corG[cordat$gsatid]


cordat$group2 <- "Stem"
cordat$group2[grepl("Somatic",cordat$group)] <- "Somatic"



topExp <- read.table("../expression/expTypeTable.tsv", header=T, sep=" ")
load("../methylation/450kMethylationData_geneLevelPromoterAverage_hit_clean.RData")
load("../methylation/450kMethylationData_geneLevelAverage_hit_clean.RData")

topMethylG <- avgMethylByGeneCleanHit
topMethylP <- avgMethylByGenePromoterCleanHit

#length(diff_shared <- intersect(rownames(topExp[1:1000,]), rownames(topMethyl[1:1000,])))
diff_sharedP <- intersect(rownames(topExp[topExp$adj.P.Val<1e-5,]), 
                          rownames(topMethylP[topMethylP$adj.P.Val<1e-5,]))
diff_sharedG <- intersect(rownames(topExp[topExp$adj.P.Val<1e-5,]), 
                          rownames(topMethylG[topMethylG$adj.P.Val<1e-5,]))


topExpP.shared <- exp[diff_sharedP,samples]
topExpG.shared <- exp[diff_sharedG,samples]
topMethylP.shared <- methylPromoter[diff_sharedP,samples]
topMethylG.shared <- methylGene[diff_sharedG,samples]

# get correlation (genome wide) for each sample between methylation values and expression values.
# cor(log(exp.shared[,3]), methyl.shared[,3], use="complete.obs")

top.genome.corP <- sapply(samples, function(x){
  cor(log(topExpP.shared[,x]), topMethylP.shared[,x], use="complete.obs")
})
top.genome.corG <- sapply(samples, function(x){
  cor(log(topExpG.shared[,x]), topMethylG.shared[,x], use="complete.obs")
})

cordat$topcorP <- top.genome.corP[cordat$gsatid]
cordat$topcorG <- top.genome.corG[cordat$gsatid]

allcordat <- data.frame(cor=c(cordat$corG, cordat$corP, cordat$topcorG, cordat$topcorP), 
                        cell=rep(cordat$group2, 4), 
                        group=c(rep("All",2*nrow(cordat)),rep("Differential",2*nrow(cordat))),
                        methyl=rep(c(rep("GeneAverage",nrow(cordat)),rep("PromoterAverage",nrow(cordat))),2))


pdf("../plots/correlations_by_cell_group_CpG.pdf", width=10, height=8)
ggplot(allcordat, aes(cell, cor, fill=cell)) + geom_boxplot(width=.2) + geom_jitter()+ facet_wrap(methyl~group)
dev.off()


pdf("../plots/correlations_by_cell_group_CpG2.pdf", width=10, height=8)
ggplot(allcordat, aes(cell, cor, fill=cell)) + geom_boxplot(width=.2) + facet_wrap(methyl~group) +
  theme(panel.background=element_rect(fill="white"))
dev.off()

#TODO: Statistical analysis of result.

#ggplot(cordat, aes(group, cor)) + geom_violin(fill="lightblue") + 
#  geom_boxplot(width=.1) + stat_summary(fun.y="mean", geom="point", color="#FF6600", size=3) + 
#  geom_point(position=position_jitter(width=.1), alpha=1/2, pch=4) +
#  ylab("Correlation") + xlab("Cell type") + ggtitle("Correlation between Expression and Promoter Methylation")


#t.test(cordat$topcor~cordat$group2)
#t.test(cordat$cor~cordat$group2)

#ggplot(cordat, aes(gsatid, topcor, color=group2)) + geom_point(stat="identity", size=3)


#ggplot(cordat, aes(group2, topcor)) + geom_boxplot(width=.2) + geom_jitter() + ylim(-.35,-.05) #+ geom_hline(aes(yintercept=0))
#ggplot(cordat, aes(group2, cor)) + geom_boxplot(width=.2) + geom_jitter() + ylim(-.35,-.05) #+ geom_hline(aes(yintercept=0))


#ggplot(cordat, aes(gsatid, cor, color=group)) + geom_point(stat="identity", size=3)


#summary(lm(cor~group, data=cordat))

#ggplot(allcordat, aes(cell, cor)) + geom_boxplot(width=.2) + geom_jitter() + facet_wrap(~group)


## trying to get contrasts to test differences i am interested in.
#contrasts:
#c(1,-1,0,0)
#c(0,0,1,-1)
#c(1,0,-1,0)
#c(0,1,0,-1)
#summary(lm(cor~0+cell+group+cell*group, data=allcordat))
