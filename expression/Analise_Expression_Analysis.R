##############
#Created by: Analise Hofmann
#Created: March 31, 2014
#Last updated: March 31, 2014
################


#Get Data: See Scott's .R for downloading the data.

#ExpDat : expression data downloaded
#expMeta : the meta data required

####Libraries:
library(preprocessCore)
library(limma)

#------------------------------------------
#log2 transform the data:
expDatlog  <- log2(expDat)

#get dimensions, did visual check of column names in expDat, and row names in expMeta 
dim(expDatlog) 
dim(expMeta) 

#check for NAs
checkExpNAs  <-  is.na(expDatlog)
omitExpNA  <- na.exclude(expDatlog)

#check dimensions against expDat
dim(omitExpNA)
#the dimensions are the same so no genes have been omitted due to NAs

#CLEANING DONE#

#-----------------------------------------------
##############
##NORMALIZATION#######

heatmap(cor(expDatlog), scale = "none", main="prenormalization Expression Correlation")
heatmap(cor(expDatlog), scale = "none", Rowv = NA, Colv = NA, main="prenormalization Expression Correlation") #saved

#Check for outliers: before normalizing quantiles
corMatrix <- cor(expDatlog)
meanCor <- data.frame(apply(corMatrix, 2, FUN = mean))
meanCor <- t(meanCor)

#colnames(meanCor) = colnames(transCombineData)
# apply(meanCor, 1, FUN = min) #prints the minimum

orderMeanCor <- order(meanCor)
orderedMeanCor <- data.frame(meanCor[, orderMeanCor])
plot(orderedMeanCor[1:100, ], main = "prenormalized Expression Correlation") #saved

#Now that the data is log2 transformed it is harder to visually remove outliers and make a judgement where the cutoff is.. 
#
#Check for outliers: after normalizing quantiles
NormExpDat  <- data.frame(normalize.quantiles(as.matrix(expDatlog))) #need to give NormExpDat column names. 
colnames(NormExpDat)  <- colnames(expDat)
row.names(NormExpDat)  <- row.names(expDat)
corMatrix <- cor(NormExpDat)
meanCor <- data.frame(apply(corMatrix, 2, FUN = mean))
meanCor <- t(meanCor)
# 
# #colnames(meanCor) = colnames(transCombineData)
# # apply(meanCor, 1, FUN = min) #prints the minimum
# 
orderMeanCor <- order(meanCor)
orderedMeanCor2 <- data.frame(meanCor[, orderMeanCor])
plot(orderedMeanCor2[1:100, ], main = "post normalization Expression Correlation") #saved

# heatmap(cor(NormExpDat), scale = "none", Rowv = NA, Colv = NA, main="post normalization Expression Correlation") #saved
# heatmap(cor(NormExpDat), scale = "none", main="post normalization Expression Correlation") #saved
# 
# #The three most apparent outliers ar GSM760199, GSM760185, and GSM60176. (The rest of the probes are within 70% correlation with the other samples)
# #Look at outliers in context of meta data:
# 
# expMeta["GSM760176",] #Somatic Diaphragm Adult Female
# expMeta["GSM760199",] #Somatic Skeletal Muscle Adult Male
# expMeta["GSM760185",] #Somatic Skeletal Muscle Adult Female
# 
# #remove the outliers, and renormalize the data:
# 
# newExpDat  <- expDatlog
# row.names(newExpDat)  <- row.names(expDat)
# newExpDat$GSM760199  <- NULL
# newExpDat$GSM760185  <- NULL
# newExpDat$GSM760176  <- NULL
# 
# NormNew.ExpDat  <- normalize.quantiles(as.matrix(newExpDat))
# row.names(NormNew.ExpDat)  <- row.names(newExpDat)
# colnames(NormNew.ExpDat)  <- colnames(newExpDat)
# 
# heatmap(cor(NormNew.ExpDat), main = "PostOutlier, Norm Exp Data", scale = "none")

##NORMALIZATION DONE##


#--------------------------------------------------
###################################
#Differential Expression Analysis##
###################################

#make model matrix

#remove rows from Meta data that correspond to our outliers that were removed. 
##newMeta <- expMeta[-c(235,221,212),]  #don't need this anymore since the outliers are not taken out yet... 

#no outliers removed at this point. 
# model.ExpDat  <- model.matrix(~characteristics_ch1.2,expMeta)
# 
# ExpFit <- lmFit(NormExpDat, model.ExpDat)
# ExpEbFit <- eBayes(ExpFit)
# ExpTissue.Table <- topTable(ExpEbFit, number = Inf, coef = "characteristics_ch1.2cell type: Somatic.Tissue")


#-------------------------------
####################################
#Try to make a more usable Meta Tabe 
#(Sample.id, cell.type, cell.line.tissue.type, gender)
####################################

#shortMeta is the condensed and decluttered version of expMeta
shortMeta  <- expMeta[,1:4]
colnames(shortMeta)  <- c("Sample.id", "cell.type", "cell.line.tissue.type", "gender")
colnames(shortMeta)
shortMeta$Sample.id  <- row.names(shortMeta)
shortMeta$cell.type  <- expMeta$characteristics_ch1.2
shortMeta$cell.line.tissue.type  <- expMeta$characteristics_ch1.3

gender  <- read.table(file = "expression/genderForMeta.txt") #I made this manually in excel, and used ch1.4-ch1.7 columns to figure out the genders
shortMeta$gender  <- gender$V1
write.table(shortMeta, file = "expression/Modified.ExpMeta.tsv", row.names = F, col.names = T)

#--------------------------------
###############################
#Retry Linear Regression Analysis with new Meta table
##############################

#Make Model Matrix
#adding the (-1) to the mode.matrix argument makes it so that there is no Intercept (Is this correct? or should the ESC be the Intercept? or ESC + PESC together should be the Intercept)
model.ExpDat2  <- model.matrix(~cell.type ,shortMeta)


ExpFit2 <- lmFit(NormExpDat, model.ExpDat2)
ExpEbFit2 <- eBayes(ExpFit2)

#I am not sure how to proceed with the TopTable here in regards to my new Meta table
ExpTissue.Table <- topTable(ExpEbFit2, number = Inf, coef = grep("cell.type", colnames(coef(ExpEbFit2))))


##############RETRY EVERYTHING WITH THE avgExpByGene#############
#log2 transform
log.avgExpByGene  <- log2(avgExpByGene)

#get dimensions, did visual check of column names in expDat, and row names in expMeta 
dim(log.avgExpByGene) #need to cut down to the 87 we are using.

dim(expMeta) #need to cut down to the 87 we are using

#check for NAs
checkExpNAs  <-  is.na(log.avgExpByGene)
log.avgExpGeneClean <- na.exclude(log.avgExpByGene)

#check dimensions against expDat
dim(log.avgExpGeneClean) #visual check confirms no more NAs


heatmap(cor(log.avgExpGeneClean), scale = "none", main="prenormalization Expression Correlation") #saved
#I don't see any crazy outliers

#shorten the Meta file to the 87 Samples.
expGDM_gsatid  <- read.table("expression/expGSM_gsatid.tsv")
geo_accession  <- as.vector(expGDM_gsatid[,"geo_accession"])
row.names(newMeta)  <- newMeta$Sample.id
expMetaClean  <- newMeta[geo_accession,] 
intersect(expMetaClean$Sample.id, expGDM_gsatid$geo_accession) #check intersect, good!
#add column to Meta file ES/Somatic for all samples.##


#Normalize
log.avgExpCleanNorm  <- data.frame(normalize.quantiles(as.matrix(log.avgExpGeneClean)))
colnames(log.avgExpCleanNorm)   <- colnames(log.avgExpGeneClean)
row.names(log.avgExpCleanNorm)  <- row.names(log.avgExpGeneClean)
heatmap(cor(log.avgExpCleanNorm), scale = "none", main=" post normalization Expression Correlation") #saved

#first need to reorder the expGSM file numerically to get the gsatid's
#add make a column for the meta data that is ES/Somatic
#load in the table made in xcel 
#newMeta  <- read.table(file = "expression/expMetaNewColumns.txt")
type  <- read.table(file = "expression/ESvsSomatic.txt")


expMetaClean$Type  <- type[,"V1"]

write.table(log.avgExpCleanNorm, file="expression/expGeneCleanedData.tsv", row.names = T, col.names = T)

#limma: differential expression analysis
model.ExpDat3  <- model.matrix(~Type,expMetaClean)


ExpFit3 <- lmFit(log.avgExpCleanNorm, model.ExpDat3)
ExpEbFit3 <- eBayes(ExpFit3)
topTableType  <- topTable(ExpEbFit3, number = Inf, coef = grep("Type", colnames(coef(ExpEbFit3))))
write.table(topTableType, file = "expression/expTypeTable.tsv", col.names = T, row.names = T)

#Make a heatmap of the top 100 genes.
heatmap(as.matrix(log.avgExpCleanNorm[c(row.names(topTableType[1:100,])),]),scale = "none", Rowv = NA, Colv = NA, main = "no clustering top 100 genes")
heatmap(as.matrix(log.avgExpCleanNorm[c(row.names(topTableType[1:100,])),]),scale = "none", Rowv = NA, main = "clustering top 100 genes")


#reorder genes in expression data based on top table.
