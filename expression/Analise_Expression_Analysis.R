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


#get dimensions, did visual check of column names in expDat, and row names in expMeta 
dim(expDat) 
dim(expMeta) 

#check for NAs
checkExpNAs  <-  is.na(expDat)
omitExpNA  <- na.exclude(expDat)

#check dimensions against expDat
dim(omitExpNA)
#the dimensions are the same so no genes have been omitted due to NAs

#CLEANING DONE#


##############
##NORMALIZATION#######

heatmap(cor(expDat), scale = "none", main="prenormalization Expression Correlation")
heatmap(cor(expDat), scale = "none", Rowv = NA, Colv = NA, main="prenormalization Expression Correlation") #saved

#Check for outliers: before normalizing quantiles
corMatrix <- cor(expDat)
meanCor <- data.frame(apply(corMatrix, 2, FUN = mean))
meanCor <- t(meanCor)

#colnames(meanCor) = colnames(transCombineData)
# apply(meanCor, 1, FUN = min) #prints the minimum

orderMeanCor <- order(meanCor)
orderedMeanCor <- data.frame(meanCor[, orderMeanCor])
plot(orderedMeanCor[1:100, ], main = "prenormalized Expression Correlation") #saved


#Check for outliers: after normalizing quantiles
NormExpDat  <- data.frame(normalize.quantiles(as.matrix(expDat))) #need to give NormExpDat column names. 
colnames(NormExpDat)  <- colnames(expDat)
row.names(NormExpDat)  <- row.names(expDat)
corMatrix <- cor(NormExpDat)
meanCor <- data.frame(apply(corMatrix, 2, FUN = mean))
meanCor <- t(meanCor)

#colnames(meanCor) = colnames(transCombineData)
# apply(meanCor, 1, FUN = min) #prints the minimum

orderMeanCor <- order(meanCor)
orderedMeanCor2 <- data.frame(meanCor[, orderMeanCor])
plot(orderedMeanCor2[1:100, ], main = "post normalization Expression Correlation") #saved

heatmap(cor(NormExpDat), scale = "none", Rowv = NA, Colv = NA, main="post normalization Expression Correlation") #saved
heatmap(cor(NormExpDat), scale = "none", main="post normalization Expression Correlation") #saved

#The three most apparent outliers ar GSM760199, GSM760185, and GSM60176. (The rest of the probes are within 70% correlation with the other samples)
#Look at outliers in context of meta data:

expMeta["GSM760176",] #Somatic Diaphragm Adult Female
expMeta["GSM760199",] #Somatic Skeletal Muscle Adult Male
expMeta["GSM760185",] #Somatic Skeletal Muscle Adult Female

#remove the outliers, and renormalize the data:

newExpDat  <- expDat
row.names(newExpDat)  <- row.names(expDat)
newExpDat$GSM760199  <- NULL
newExpDat$GSM760185  <- NULL
newExpDat$GSM760176  <- NULL

NormNew.ExpDat  <- normalize.quantiles(as.matrix(newExpDat))
row.names(NormNew.ExpDat)  <- row.names(newExpDat)
colnames(NormNew.ExpDat)  <- colnames(newExpDat)

heatmap(cor(NormNew.ExpDat), main = "PostOutlier, Norm Exp Data", scale = "none")

