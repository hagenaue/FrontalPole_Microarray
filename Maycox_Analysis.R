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
PMDlabelremoved<-sub("post-mortem delay: ","", Meta(GSMList(gse)$GSM439778)$characteristics_ch1[3], fixed = T)
as.numeric(substr(PMDlabelremoved, 1, nchar(PMDlabelremoved)-1))
rm(PMDlabelremoved)
as.numeric(sub("ph: ","", Meta(GSMList(gse)$GSM439778)$characteristics_ch1[4]))

SampleID<-as.matrix(names(GSMList(gse)))
Gender<-matrix("a", nrow=51, ncol=1)
Age<-matrix(0, nrow=51, ncol=1)
BrainpH<-matrix(0, nrow=51, ncol=1)
PMD<-matrix(0, nrow=51, ncol=1)
Diagnosis2<-matrix("a", nrow=51, ncol=1)


i<-1
for(i in c(1:51)){
  Gender[i]<-sub("gender: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[1])
  Age[i]<-as.numeric(sub("age: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[2]))
  PMDlabelremoved<-sub("post-mortem delay: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[3], fixed = T)
  PMD[i]<-as.numeric(substr(PMDlabelremoved, 1, nchar(PMDlabelremoved)-1))
  rm(PMDlabelremoved)
  BrainpH[i]<-as.numeric(sub("ph: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[4]))
  if(all(substr(Meta(GSMList(gse)[[i]])$title, 6,8)=="Scz")) Diagnosis2[i]<-"Scz" else Diagnosis2[i]<-"Control"
}

# This entry does not have a gender listed, so the data needs to be assigned outside of the for loop
Gender[11]<-NA
Age[11]<-as.numeric(sub("age: ","", Meta(GSMList(gse)$GSM439788)$characteristics_ch1[1]))
PMDlabelremoved<-sub("post-mortem delay: ","", Meta(GSMList(gse)$GSM439788)$characteristics_ch1[2], fixed = T)
PMD[11]<-as.numeric(substr(PMDlabelremoved, 1, nchar(PMDlabelremoved)-1))
rm(PMDlabelremoved)
BrainpH[11]<-as.numeric(sub("ph: ","", Meta(GSMList(gse)$GSM439788)$characteristics_ch1[3]))

Diagnosis<-Diagnosis2
Diagnosis<-relevel(as.factor(Diagnosis), ref="Control")

Maycox_SampleCharacteristics<-data.frame(SampleID, Gender, Age, PMD, BrainpH, Diagnosis, stringsAsFactors=F)

head(Maycox_SampleCharacteristics)

write.csv(Maycox_SampleCharacteristics, "Maycox_SampleCharacteristics.csv")

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
ScanDateDayOnly<-relevel(as.factor(ScanDateDayOnly), ref="01/21/04")

#Converts the data2 AffyBatch into an ExpressionSet object using the robust multi-array average (RMA) expression measure. 
#The expression measure is given in log base 2 scale.
eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDFplus2.txt")
RMAExpression_customCDFplus2<-read.delim("data_customCDFplus2.txt", sep="\t")
str(RMAExpression_customCDFplus2)
#'data.frame':	8551 obs. of  51 variables
write.csv(RMAExpression_customCDFplus2, "RMAExpression_customCDFplus2.csv")

BiocManager::install("AffyRNADegradation")
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

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Graphs")


#Variable Analysis
SubjectFactorVariables<-cbind(Diagnosis, Gender, ScanDateDayOnly)
colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender", "ScanDateDayOnly")

SubjectContinuousVariables<-cbind(BrainpH, PMD, Age, RNADegradPerSample)
colnames(SubjectContinuousVariables)<-c("BrainpH", "PMD", "Age", "RNADegradPerSample")

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
    vif(lm(SignalSortedNoNA3[1,]~BrainpH+PMD+Diagnosis+Gender+Age+RNADegradPerSample+ScanDateDayOnly))
    
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