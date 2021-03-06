####################################################
#### Explore Expression/Methylation Correlation ####
#### Author(s): Scott Brown
#### Date Created: March 31, 2014
#### Last Edited by: Scott
#### on: April 6, 2014
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

#testing m-values
library(wateRmelon)
methylPromoter <- beta2m(betaqn(avgMethylByGenePromoterClean))
methylGene <- beta2m(betaqn(avgMethylByGeneClean))

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
  cor(expP.shared[,x], methylP.shared[,x], use="complete.obs")
})

genome.corG <- sapply(samples, function(x){
  cor(expG.shared[,x], methylG.shared[,x], use="complete.obs")
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
  cor(topExpP.shared[,x], topMethylP.shared[,x], use="complete.obs")
})
top.genome.corG <- sapply(samples, function(x){
  cor(topExpG.shared[,x], topMethylG.shared[,x], use="complete.obs")
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

a <- with(subset(allcordat, group=="All" & methyl=="GeneAverage"),
     t.test(cor[cell=="Somatic"], cor[cell=="Stem"]))$p.val

b <- with(subset(allcordat, group=="Differential" & methyl=="GeneAverage"),
     t.test(cor[cell=="Somatic"], cor[cell=="Stem"]))$p.val

c <- with(subset(allcordat, group=="All" & methyl=="PromoterAverage"),
     t.test(cor[cell=="Somatic"], cor[cell=="Stem"]))$p.val

d <- with(subset(allcordat, group=="Differential" & methyl=="PromoterAverage"),
     t.test(cor[cell=="Somatic"], cor[cell=="Stem"]))$p.val

#----------------

e <- with(subset(allcordat, methyl=="GeneAverage" & cell=="Somatic"),
     t.test(cor[group=="All"], cor[group=="Differential"]))$p.val

f <- with(subset(allcordat, methyl=="GeneAverage" & cell=="Stem"),
     t.test(cor[group=="All"], cor[group=="Differential"]))$p.val

g <- with(subset(allcordat, methyl=="PromoterAverage" & cell=="Somatic"),
     t.test(cor[group=="All"], cor[group=="Differential"]))$p.val

h <- with(subset(allcordat, methyl=="PromoterAverage" & cell=="Stem"),
     t.test(cor[group=="All"], cor[group=="Differential"]))$p.val

#----------------

i <- with(subset(allcordat, group=="All" & cell=="Somatic"),
     t.test(cor[methyl=="GeneAverage"], cor[methyl=="PromoterAverage"]))$p.val

j <- with(subset(allcordat, group=="All" & cell=="Stem"),
     t.test(cor[methyl=="GeneAverage"], cor[methyl=="PromoterAverage"]))$p.val

k <- with(subset(allcordat, group=="Differential" & cell=="Somatic"),
     t.test(cor[methyl=="GeneAverage"], cor[methyl=="PromoterAverage"]))$p.val

l <- with(subset(allcordat, group=="Differential" & cell=="Stem"),
     t.test(cor[methyl=="GeneAverage"], cor[methyl=="PromoterAverage"]))$p.val


t.test.res <- data.frame(Comparison=c("Cell","Cell","Cell","Cell",
                                      "Group","Group","Group","Group",
                                      "Methyl_strategy","Methyl_strategy","Methyl_strategy","Methyl_strategy"),
                         Context=c("All_GeneAverage","Differential_GeneAverage","All_PromoterAverage","Differential_PromoterAverage",
                                   "GeneAverage_Somatic","GeneAverage_Stem","PromoterAverage_Somatic","PromoterAverage_Stem",
                                   "All_Somatic","All_Stem","Differential_Somatic","Differential_Stem"),
                         p.value=c(a,b,c,d,e,f,g,h,i,j,k,l))

library(knitr)
kable(t.test.res)

# Venn Diagram Code (Nat)

# Plotting Venn Diagram to compare (raw) expressed vs methylation
library(VennDiagram)
venn.combine.1 <- list(Promoter_methylation = rownames(methylPromoter), 
                       Gene_methylation = rownames(methylGene), 
                       Gene_expression = rownames(exp))
pdf("../plots/gene_overlap_venn.pdf")
plot.new()
venn.plot.1 <- venn.diagram(venn.combine.1, filename = NULL, fill = c('red', 'blue', 'green'))
grid.draw(venn.plot.1)
dev.off()

# Plotting Venn Diagram to compare differentially expressed vs methylation
venn.combine.2 <- list(Differential_Promoter_methylation = rownames(topMethylP[topMethylP$adj.P.Val<1e-5,]), 
                       Differential_Gene_methylation = rownames(topMethylG[topMethylG$adj.P.Val<1e-5,]), 
                       Differential_expression = rownames(topExp[topExp$adj.P.Val<1e-5,]))
pdf("../plots/differential_gene_overlap_venn.pdf")
plot.new()
venn.plot.2 <- venn.diagram(venn.combine.2, filename = NULL, fill = c('red', 'blue', 'green'))
grid.draw(venn.plot.2)
dev.off()


#### Example correlation xyplot
library(lattice)
pdf("../plots/correlation_example.pdf")
xyplot(topExpG.shared[,1]~topMethylG.shared[,1], type=c("p","r"), xlab="Mehtylation M Values",
       ylab="Log(expression)", main="Relationship between methylation and expression")
dev.off()
xyplot(topExpP.shared[,2]~topMethylP.shared[,2], type=c("p","r"))


library(ggplot2)
ggdat <- data.frame(exp=topExpP.shared[,2], meth=topMethylP.shared[,2])

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}


pdf("../plots/example_correlation.pdf")
ggplot(ggdat, aes(meth, exp)) + geom_point(alpha=0.5) + geom_smooth(method='lm') +
  annotate("text", x = 3, y = 14, label = lm_eqn(lm(exp ~ meth, ggdat)), colour="black", size = 5, parse=TRUE) +
  xlab("Methylation") + ylab("log(Expression)") + ggtitle("Correlation between Methylation and Expression")
dev.off()
# fit <- lm(cor~0+cell+group+methyl+cell*group+cell*methyl+group*methyl+cell*group*methyl, allcordat)
# summary(fit)
# 
# cont.matrix <- makeContrasts()
# 
# ###############
# (desmat <- model.matrix(~cell+group+methyl+cell*group+cell*methyl+group*methyl+cell*group*methyl, allcordat))
# colnames(desmat)
# 
# cont.matrix <- makeContrasts(a=cellStem-Intercept,
#                              b=cellStem:groupDifferential-groupDifferential,
#                              c=methylPromoterAverage*cellStem-methylPromoterAverage,
#                              d=methylPromoterAverage*cellStem*groupDifferential-methylPromoterAverage*groupDifferential,
#                              levels=desmat)

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

