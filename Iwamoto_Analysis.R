#Iwamoto Analysis Code
#Liam

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654")

library(GEOquery)
gse <- getGEO("GSE12654", GSEMatrix = FALSE)

head(Meta(gse))
str(Meta(gse))
Meta(GSMList(gse)$GSM317649)$characteristics_ch1
# [1] "prefrontal cortex (BA10)" "control"


New_Iwamoto_MetaData<-read.csv("demo_kato_short.csv", header=T, stringsAsFactors=FALSE)

SampleID<-as.matrix(New_Iwamoto_MetaData$Trimmed_Sample_id)
Gender<-as.matrix(New_Iwamoto_MetaData$Trimmed_Sex)
Age<-as.matrix(New_Iwamoto_MetaData$Trimmed_Age)
BrainpH<-as.matrix(New_Iwamoto_MetaData$Trimmed_Brain.PH)
PMI<-as.matrix(New_Iwamoto_MetaData$Trimmed_PMI..h.)
Diagnosis2<-as.matrix(New_Iwamoto_MetaData$Trimmed_Profile)
LeftBrainStatus<-as.matrix(New_Iwamoto_MetaData$Trimmed_Left.Brain)
SuicideStatus<-as.matrix(New_Iwamoto_MetaData$Trimmed_Suicide.Status)
PsychoticFeature<-as.matrix(New_Iwamoto_MetaData$Trimmed_Psychotic.Feature)
RateofDeath<-as.matrix(New_Iwamoto_MetaData$Trimmed_Rate.of.Death)
SmokingatTimeofDeath<-as.matrix(New_Iwamoto_MetaData$Trimmed_Smoking.at.Time.of.Death)
LifetimeAlcohol<-as.matrix(New_Iwamoto_MetaData$Trimmed_Lifetime.Alcohol)
LifetimeDrugs<-as.matrix(New_Iwamoto_MetaData$Trimmed_Lifetime.Drugs)

Diagnosis<-Diagnosis2
Diagnosis<-relevel(as.factor(Diagnosis), ref="Control")

Iwamoto_SampleCharacteristics<-data.frame(SampleID, Gender, Age, BrainpH, PMI, Diagnosis, LeftBrainStatus, SuicideStatus, PsychoticFeature, RateofDeath, SmokingatTimeofDeath, LifetimeAlcohol, LifetimeDrugs)

head(Iwamoto_SampleCharacteristics)

library(org.Hs.eg.db)
library(plyr)
library(affy)

#CDF and the chip.
install.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/pd.hgu95av2.hs.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/hgu95av2hsentrezgprobe_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/hgu95av2hsentrezgcdf_25.0.0.tar.gz", repos = NULL, type = "source")

#Reads .CEL files into an Affybatch
data2<-ReadAffy(cdfname ="hgu95av2hsentrezg")
str(data2)
data2

ScanDate2<-protocolData(data2)$ScanDate
ScanDateDayOnly2<-matrix("a", nrow=50, ncol=1)
ScanDate<-matrix("a", nrow=50, ncol=1)

i<-1
for(i in c(1:50)){
  ScanDateDayOnly2[i]<-substr(ScanDate2[i], 1, 8)
  ScanDate[i]<-ScanDate2[i]
}

#ScanDate<-ScanDate2
#ScanDate<-relevel(as.factor(ScanDate2), ref= "08/24/01 15:28:22")
#ScanDateDayOnly<-ScanDateDayOnly2
#ScanDateDayOnly<-relevel(as.factor(ScanDateDayOnly), ref="08/24/01")

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

sum(is.na(RMAExpression_customCDFAnnotation2plus2[,3])==F)
#[1] 8484
dim(RMAExpression_customCDFAnnotation2plus2)
#[1] 8551     3
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

SignalSortedNoNA3<-as.matrix(RMAExpression_customCDFplus2[,-1])

cbind(SampleID, colnames(SignalSortedNoNA3))

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Graphs")

#Quality Control

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",]

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]

#checking for gender switches
#not working for some reason
png("XIST_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]~Gender, col=2)
dev.off()
#

#Variable Analysis
SubjectFactorVariables<-cbind(Diagnosis, Gender, LeftBrainStatus, SuicideStatus, PsychoticFeature, RateofDeath, SmokingatTimeofDeath, LifetimeAlcohol, LifetimeDrugs, ScanDate)
colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender", "LeftBrainStatus", "SuicideStatus", "PsychoticFeature", "RateofDeath", "SmokingatTimeofDeath", "LifetimeAlcohol", "LifetimeDrugs", "ScanDate")

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
    summary(ScanDate),
    summary(LeftBrainStatus),
    summary(SuicideStatus),
    summary(PsychoticFeature),
    summary(RateofDeath),
    summary(SmokingatTimeofDeath),
    summary(LifetimeAlcohol),
    summary(LifetimeDrugs),
    
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
    vif(lm(SignalSortedNoNA3[1,]~BrainpH+PMI+Diagnosis+Gender+Age+LeftBrainStatus+SuicideStatus+PsychoticFeature+RateofDeath+SmokingatTimeofDeath+LifetimeAlcohol+LifetimeDrugs+RNADegradPerSample+ScanDate))
    
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