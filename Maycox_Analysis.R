#Maycox Analysis Code
#Liam

# TO-DO:
# 1) Diagnosis is not included in characteristics;I see it on GSE tho, so I need to either manually input it or write up some code to extract it

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612")

library(GEOquery)
gse <- getGEO("GSE17612", GSEMatrix = FALSE)

head(Meta(gse))
str(Meta(gse))
Meta(GSMList(gse)$GSM439778)$characteristics_ch1
# [1] "gender: Male"            "age: 74"                 "post-mortem delay: 4.5h"
# [4] "ph: 6"  

# Diagnosis is not included in characteristics, so I extract it later on using the samples' titles
Meta(GSMList(gse)$GSM439778)$title
# [1] "S014_Scz_M_74"

sub("gender: ","", Meta(GSMList(gse)$GSM439778)$characteristics_ch1[1])
as.numeric(sub("age: ","", Meta(GSMList(gse)$GSM439778)$characteristics_ch1[2]))
PMIlabelremoved<-sub("post-mortem delay: ","", Meta(GSMList(gse)$GSM439778)$characteristics_ch1[3], fixed = T)
as.numeric(substr(PMIlabelremoved, 1, nchar(PMIlabelremoved)-1))
rm(PMIlabelremoved)
as.numeric(sub("ph: ","", Meta(GSMList(gse)$GSM439778)$characteristics_ch1[4]))

SampleID<-as.matrix(names(GSMList(gse)))
Gender<-matrix("a", nrow=51, ncol=1)
Age<-matrix(0, nrow=51, ncol=1)
BrainpH<-matrix(0, nrow=51, ncol=1)
PMI<-matrix(0, nrow=51, ncol=1)
Diagnosis2<-matrix("a", nrow=51, ncol=1)


i<-1
for(i in c(1:51)){
  Gender[i]<-sub("gender: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[1])
  Age[i]<-as.numeric(sub("age: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[2]))
  PMIlabelremoved<-sub("post-mortem delay: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[3], fixed = T)
  PMI[i]<-as.numeric(substr(PMIlabelremoved, 1, nchar(PMIlabelremoved)-1))
  rm(PMIlabelremoved)
  BrainpH[i]<-as.numeric(sub("ph: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[4]))
  if(all(substr(Meta(GSMList(gse)[[i]])$title, 6,8)=="Scz")) Diagnosis2[i]<-"Scz" else Diagnosis2[i]<-"Control"
}

# This entry does not have a gender listed, so the data needs to be assigned outside of the for loop
Gender[11]<-NA
Age[11]<-as.numeric(sub("age: ","", Meta(GSMList(gse)$GSM439788)$characteristics_ch1[1]))
PMIlabelremoved<-sub("post-mortem delay: ","", Meta(GSMList(gse)$GSM439788)$characteristics_ch1[2], fixed = T)
PMI[11]<-as.numeric(substr(PMIlabelremoved, 1, nchar(PMIlabelremoved)-1))
rm(PMIlabelremoved)
BrainpH[11]<-as.numeric(sub("ph: ","", Meta(GSMList(gse)$GSM439788)$characteristics_ch1[3]))

Diagnosis<-Diagnosis2
Diagnosis<-relevel(as.factor(Diagnosis), ref="Control")



library(org.Hs.eg.db)
library(plyr)
library(affy)

# CDF and the chip.
install.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/pd.hgu133plus2.hs.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/hgu133plus2hsentrezgprobe_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/hgu133plus2hsentrezgcdf_25.0.0.tar.gz", repos = NULL, type = "source")
# Reads .CEL files into an Affybatch
data2<-ReadAffy(cdfname ="hgu133plus2hsentrezg")
str(data2)
data2
# AffyBatch object
# size of arrays=1164x1164 features (40 kb)
# cdf=hgu133plus2hsentrez (??? affyids)
# number of samples=51

ScanDate2<-protocolData(data2)$ScanDate
ScanDateDayOnly2<-matrix("a", nrow=50, ncol=1)

i<-1
for(i in c(1:51)){
  ScanDateDayOnly2[i]<-substr(ScanDate2[i], 1, 8)
}

ScanDate<-ScanDate2
ScanDate<-relevel(as.factor(ScanDate2), ref= "01/21/04 10:44:01")
ScanDateDayOnly<-ScanDateDayOnly2
ScanDateDayOnly<-relevel(as.factor(ScanDateDayOnly), ref="01/22/04")

Maycox_SampleCharacteristics<-data.frame(SampleID, Gender, Age, PMI, BrainpH, Diagnosis, ScanDateDayOnly, stringsAsFactors=F)

head(Maycox_SampleCharacteristics)

write.csv(Maycox_SampleCharacteristics, "Maycox_SampleCharacteristics.csv")

#Converts the data2 AffyBatch into an ExpressionSet object using the robust multi-array average (RMA) expression measure. 
#The expression measure is given in log base 2 scale.
eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDFplus2.txt")
RMAExpression_customCDFplus2<-read.delim("data_customCDFplus2.txt", sep="\t")
str(RMAExpression_customCDFplus2)
#'data.frame':	8551 obs. of  51 variables
write.csv(RMAExpression_customCDFplus2, "RMAExpression_customCDFplus2.csv")

library(AffyRNADegradation)

#Generate and visualize tongs plots
tongs <- GetTongs(data2, chip.idx = 4)
PlotTongs(tongs)
tongs <- GetTongs(data2, chip.idx = 5)
PlotTongs(tongs)

rna.deg<- RNADegradation(data2, location.type = "index")
RNADegradPerSample<-d(rna.deg)
str(RNADegradPerSample)

head(RMAExpression_customCDFplus2)
RMAExpression_EntrezIDplus2<-sub("_at", "", RMAExpression_customCDFplus2[,1])
head(RMAExpression_EntrezIDplus2)
RMAExpression_customCDFAnnotationplus2<-data.frame(RMAExpression_customCDFplus2[,1], RMAExpression_EntrezIDplus2, stringsAsFactors = F )
colnames(RMAExpression_customCDFAnnotationplus2)<-c("ProbesetID", "EntrezGeneID")

x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1
# 61217

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

RMAExpression_customCDFAnnotation2plus2<-join(RMAExpression_customCDFAnnotationplus2, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")

write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

SignalSortedNoNA3<-as.matrix(RMAExpression_customCDFplus2[,-1])

cbind(SampleID, colnames(SignalSortedNoNA3))

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Graphs")

#Quality Control

RMAExpression_customCDFAnnotation2plus2[which(RMAExpression_customCDFAnnotation2plus2[,3]=="XIST"),]

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]

#checking for gender switches
png("XIST_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="XIST"),]~Gender, col=2)
dev.off()

png("RPS4Y1_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="RPS4Y1"),]~Gender, col=2)
dev.off()


png("DDX3Y_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="DDX3Y"),]~Gender, col=2)
dev.off()


png("RPS4Y1vsXIST_GenderCheck.png")
plot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="RPS4Y1"),]~SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="XIST"),], col=as.numeric(as.factor(Gender)))
dev.off()

png("DDX3YvsXIST_GenderCheck.png")
plot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="DDX3Y"),]~SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]=="XIST"),], col=as.numeric(as.factor(Gender)))
dev.off()



##Removing bad samples
KnownBadSamples<-c("GSM439786","GSM439795")

SampleCharacteristics_NoOutliers<-Maycox_SampleCharacteristics[(SampleID%in%KnownBadSamples)==F,]
dim(SampleCharacteristics_NoOutliers)

SignalSortedNoNA3NoOutliers<-SignalSortedNoNA3[,(SampleID%in%KnownBadSamples)==F]
dim(SignalSortedNoNA3NoOutliers)

RNADegradPerSampleNoOutliers<-RNADegradPerSample[(SampleID%in%KnownBadSamples)==F]


#Redefining the variables to only include good data:
SignalSortedNoNA3<-SignalSortedNoNA3NoOutliers

Gender<-SampleCharacteristics_NoOutliers$Gender
Age<-SampleCharacteristics_NoOutliers$Age
PMI<-SampleCharacteristics_NoOutliers$PMI
BrainpH<-SampleCharacteristics_NoOutliers$BrainpH
Diagnosis<-SampleCharacteristics_NoOutliers$Diagnosis
RNADegradPerSample<-RNADegradPerSampleNoOutliers
ScanDateDayOnly<-SampleCharacteristics_NoOutliers$ScanDateDayOnly

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/RemovedOutliers/Graphs")


#Variable Analysis
SubjectFactorVariables<-cbind(Diagnosis, Gender, ScanDateDayOnly)
colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender", "ScanDateDayOnly")

SubjectContinuousVariables<-cbind(BrainpH, PMI, Age, RNADegradPerSample)
colnames(SubjectContinuousVariables)<-c("BrainpH", "PMI", "Age", "RNADegradPerSample")

for (i in 1:length(SubjectContinuousVariables[1,])){
  png(paste(paste("Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep="."))	
  hist(SubjectContinuousVariables[, i], col=i+1)
  dev.off()		
}

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}

#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}

#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("Cross Tabs Between Subject Factors.txt")
out<-c(
  capture.output(
    
    summary(Diagnosis),
    summary(Gender),
    summary(ScanDateDayOnly),
    
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(ContingencyTable)
        print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
      }		
    }
  )
)
cat(out, file="Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)

library(car)

StatisticalRelationshipsIV<-file("Statistical Relationships between Subject Variables.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. 
    vif(lm(SignalSortedNoNA3[1,]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample+ScanDateDayOnly))
    
  ),
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(chisq.test(ContingencyTable))		
      }		
    }
  )
  
)
cat(out, file="Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)

#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("Flagged Relationships Between Subject Variables.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
        if(chisq.test(ContingencyTable)$p.value<0.05){
          print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))
        }else{}
      }
    }
  )
)
cat(out, file="Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)

#PCA
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/RemovedOutliers/PCA")



pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")


PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, RMAExpression_customCDFAnnotation2plus2)
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

#Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()


#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off()

#After outlier removal:
png("10 PC1 vs PC2_byDiagnosis.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Diagnosis))
dev.off()


# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off()

#After outlier removal:
png("10 PC3 vs PC4_byGender.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Gender))
dev.off()

png("10 PC3 vs PC4_byDiagnosis.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Diagnosis))
dev.off()

SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

PCAoutput<-cbind(SubjectFactorVariables, SubjectContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#Visualize the sample-sample correlations using a heatmap:
png("09 Sample Sample Correlations Heatmap.png")
image(cor(SignalSortedNoNA3), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=1000, height=600)
boxplot(data.frame(cor(SignalSortedNoNA3)), cex=0.25, las=3, par(cex.axis=0.75, mar=c(15,4,4,4)), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()


##Plotting PCA vs Subject Variables
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/RemovedOutliers/Graphs/PCA_vs_Subject_Variables")

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectPCA[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
    RegressionLine<-lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}

#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectPCA[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectPCA)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}

#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("15 Statistical Relationships between Subject Variables and PCA.txt")
out<-c(
 
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
    summary.lm(lm(SubjectPCA[,1]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample+ScanDateDayOnly))
  ),
  
  capture.output(
    summary.lm(lm(SubjectPCA[,2]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample+ScanDateDayOnly))
  ),
  
  capture.output(
    summary.lm(lm(SubjectPCA[,3]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample+ScanDateDayOnly))
  ),
  
  capture.output(
    summary.lm(lm(SubjectPCA[,4]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample+ScanDateDayOnly))
  ),
  
  
  #Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between PCA and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  )
  
)
cat(out, file="15 Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)