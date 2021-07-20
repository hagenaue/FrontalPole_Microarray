#fGSEA analysis using Iwamoto and Maycox limma results for models 1 and 2

library(fgsea)
library(plyr)

#reading in data and adding annotation
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")

Iwamoto_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)

Iwamoto_Model1_NoAnnotation<-read.csv("Limma_results_Model_onlyDiagnosis.csv", header=TRUE, stringsAsFactors = FALSE)
Iwamoto_Model1<-cbind.data.frame(Iwamoto_Model1_NoAnnotation, Iwamoto_Annotation)
dim(Iwamoto_Model1)
colnames(Iwamoto_Model1)
str(Iwamoto_Model1)

Iwamoto_Model2_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath.csv", header=TRUE, stringsAsFactors = FALSE)
Iwamoto_Model2<-cbind.data.frame(Iwamoto_Model2_NoAnnotation, Iwamoto_Annotation)
dim(Iwamoto_Model2)
colnames(Iwamoto_Model2)
str(Iwamoto_Model2)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")

Maycox_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)

Maycox_Model1_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH.csv", header=TRUE, stringsAsFactors = FALSE)
Maycox_Model1<-cbind.data.frame(Maycox_Model1_NoAnnotation, Maycox_Annotation)
dim(Maycox_Model1)
colnames(Maycox_Model1)
str(Maycox_Model1)

Maycox_Model2_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI.csv", header=TRUE, stringsAsFactors = FALSE)
Maycox_Model2<-cbind.data.frame(Maycox_Model2_NoAnnotation, Maycox_Annotation)
dim(Maycox_Model2)
colnames(Maycox_Model2)
str(Maycox_Model2)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/fGSEA")

#removing NAs
sum(is.na(Iwamoto_Model1$GeneSymbol))
# [1] 67
Iwamoto_Model1_noNA<-Iwamoto_Model1[is.na(Iwamoto_Model1$GeneSymbol)==FALSE,]
dim(Iwamoto_Model1_noNA)

sum(is.na(Iwamoto_Model2$GeneSymbol))
# [1] 67
Iwamoto_Model2_noNA<-Iwamoto_Model2[is.na(Iwamoto_Model2$GeneSymbol)==FALSE,]
dim(Iwamoto_Model2_noNA)


sum(is.na(Maycox_Model1$GeneSymbol))
# [1] 68
Maycox_Model1_noNA<-Maycox_Model1[is.na(Maycox_Model1$GeneSymbol)==FALSE,]
dim(Maycox_Model1_noNA)

sum(is.na(Maycox_Model2$GeneSymbol))
# [1] 68
Maycox_Model2_noNA<-Maycox_Model2[is.na(Maycox_Model2$GeneSymbol)==FALSE,]
dim(Maycox_Model2_noNA)


#determining if there are duplicated symbols
sum(duplicated(Iwamoto_Model1_noNA$GeneSymbol))
# [1] 0
sum(duplicated(Iwamoto_Model2_noNA$GeneSymbol))
# [1] 0
sum(duplicated(Maycox_Model1_noNA$GeneSymbol))
# [1] 0
sum(duplicated(Maycox_Model2_noNA$GeneSymbol))
# [1] 0

#Bipolar Disorder Preparation
Iwamoto_Model1_BP_log2FC_forGSEA<-tapply(X=Iwamoto_Model1_noNA$Coef.DiagnosisFactorBipolar, INDEX=Iwamoto_Model1_noNA$GeneSymbol, FUN=mean)
Iwamoto_Model2_BP_log2FC_forGSEA<-tapply(X=Iwamoto_Model2_noNA$Coef.DiagnosisFactorBipolar, INDEX=Iwamoto_Model2_noNA$GeneSymbol, FUN=mean)

names(Iwamoto_Model1_BP_log2FC_forGSEA)<-names(table(Iwamoto_Model1_noNA$GeneSymbol))
names(Iwamoto_Model2_BP_log2FC_forGSEA)<-names(table(Iwamoto_Model2_noNA$GeneSymbol))

length(Iwamoto_Model1_BP_log2FC_forGSEA)
# [1] 8484
length(Iwamoto_Model2_BP_log2FC_forGSEA)
# [1] 8484

Iwamoto_Model1_BP_log2FC_forGSEA_Ranked<-Iwamoto_Model1_BP_log2FC_forGSEA[order(Iwamoto_Model1_BP_log2FC_forGSEA)]
head(Iwamoto_Model1_BP_log2FC_forGSEA_Ranked)
Iwamoto_Model2_BP_log2FC_forGSEA_Ranked<-Iwamoto_Model2_BP_log2FC_forGSEA[order(Iwamoto_Model2_BP_log2FC_forGSEA)]
head(Iwamoto_Model2_BP_log2FC_forGSEA_Ranked)

#Depression Preparation
Iwamoto_Model1_Depression_log2FC_forGSEA<-tapply(X=Iwamoto_Model1_noNA$Coef.DiagnosisFactorDepression, INDEX=Iwamoto_Model1_noNA$GeneSymbol, FUN=mean)
Iwamoto_Model2_Depression_log2FC_forGSEA<-tapply(X=Iwamoto_Model2_noNA$Coef.DiagnosisFactorDepression, INDEX=Iwamoto_Model2_noNA$GeneSymbol, FUN=mean)

names(Iwamoto_Model1_Depression_log2FC_forGSEA)<-names(table(Iwamoto_Model1_noNA$GeneSymbol))
names(Iwamoto_Model2_Depression_log2FC_forGSEA)<-names(table(Iwamoto_Model2_noNA$GeneSymbol))

length(Iwamoto_Model1_Depression_log2FC_forGSEA)
# [1] 8484
length(Iwamoto_Model2_Depression_log2FC_forGSEA)
# [1] 8484

Iwamoto_Model1_Depression_log2FC_forGSEA_Ranked<-Iwamoto_Model1_Depression_log2FC_forGSEA[order(Iwamoto_Model1_Depression_log2FC_forGSEA)]
head(Iwamoto_Model1_Depression_log2FC_forGSEA_Ranked)
Iwamoto_Model2_Depression_log2FC_forGSEA_Ranked<-Iwamoto_Model2_Depression_log2FC_forGSEA[order(Iwamoto_Model2_Depression_log2FC_forGSEA)]
head(Iwamoto_Model2_Depression_log2FC_forGSEA_Ranked)

#Schizophrenia Preparation
Iwamoto_Model1_Schiz_log2FC_forGSEA<-tapply(X=Iwamoto_Model1_noNA$Coef.DiagnosisFactorSchizophrenia, INDEX=Iwamoto_Model1_noNA$GeneSymbol, FUN=mean)
Iwamoto_Model2_Schiz_log2FC_forGSEA<-tapply(X=Iwamoto_Model2_noNA$Coef.DiagnosisFactorSchizophrenia, INDEX=Iwamoto_Model2_noNA$GeneSymbol, FUN=mean)

names(Iwamoto_Model1_Schiz_log2FC_forGSEA)<-names(table(Iwamoto_Model1_noNA$GeneSymbol))
names(Iwamoto_Model2_Schiz_log2FC_forGSEA)<-names(table(Iwamoto_Model2_noNA$GeneSymbol))

length(Iwamoto_Model1_Schiz_log2FC_forGSEA)
# [1] 8484
length(Iwamoto_Model2_Schiz_log2FC_forGSEA)
# [1] 8484

Iwamoto_Model1_Schiz_log2FC_forGSEA_Ranked<-Iwamoto_Model1_Schiz_log2FC_forGSEA[order(Iwamoto_Model1_Schiz_log2FC_forGSEA)]
head(Iwamoto_Model1_Schiz_log2FC_forGSEA_Ranked)
Iwamoto_Model2_Schiz_log2FC_forGSEA_Ranked<-Iwamoto_Model2_Schiz_log2FC_forGSEA[order(Iwamoto_Model2_Schiz_log2FC_forGSEA)]
head(Iwamoto_Model2_Schiz_log2FC_forGSEA_Ranked)


Maycox_Model1_Schiz_log2FC_forGSEA<-tapply(X=Maycox_Model1_noNA$Coef.DiagnosisFactorScz, INDEX=Maycox_Model1_noNA$GeneSymbol, FUN=mean)
Maycox_Model2_Schiz_log2FC_forGSEA<-tapply(X=Maycox_Model2_noNA$Coef.DiagnosisFactorScz, INDEX=Maycox_Model2_noNA$GeneSymbol, FUN=mean)

names(Maycox_Model1_Schiz_log2FC_forGSEA)<-names(table(Maycox_Model1_noNA$GeneSymbol))
names(Maycox_Model2_Schiz_log2FC_forGSEA)<-names(table(Maycox_Model2_noNA$GeneSymbol))

length(Maycox_Model1_Schiz_log2FC_forGSEA)
# [1] 8484
length(Maycox_Model2_Schiz_log2FC_forGSEA)
# [1] 8484

Maycox_Model1_Schiz_log2FC_forGSEA_Ranked<-Maycox_Model1_Schiz_log2FC_forGSEA[order(Maycox_Model1_Schiz_log2FC_forGSEA)]
head(Maycox_Model1_Schiz_log2FC_forGSEA_Ranked)
Maycox_Model2_Schiz_log2FC_forGSEA_Ranked<-Maycox_Model2_Schiz_log2FC_forGSEA[order(Maycox_Model2_Schiz_log2FC_forGSEA)]
head(Maycox_Model2_Schiz_log2FC_forGSEA_Ranked)


#Reading in updated GMT
GMT<-gmtPathways("c5.go.v7.3.symbols_BrainInABlender_Psych.gmt.txt")
str(GMT)


#bipolar disorder results
Iwamoto_Model1_BP_FGSEAResults<-fgsea(GMT, Iwamoto_Model1_BP_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Iwamoto_Model1_BP_FGSEAResults)
Iwamoto_Model1_BP_FGSEAResults$leadingEdge<-vapply(Iwamoto_Model1_BP_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Iwamoto_Model1_BP_FGSEAResults, "Iwamoto_Model1_BipolarDisorder_FGSEAResults.csv")

Iwamoto_Model2_BP_FGSEAResults<-fgsea(GMT, Iwamoto_Model2_BP_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Iwamoto_Model2_BP_FGSEAResults)
Iwamoto_Model2_BP_FGSEAResults$leadingEdge<-vapply(Iwamoto_Model2_BP_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Iwamoto_Model2_BP_FGSEAResults, "Iwamoto_Model2_BipolarDisorder_FGSEAResults.csv")

#depression results
Iwamoto_Model1_Depression_FGSEAResults<-fgsea(GMT, Iwamoto_Model1_Depression_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Iwamoto_Model1_Depression_FGSEAResults)
Iwamoto_Model1_Depression_FGSEAResults$leadingEdge<-vapply(Iwamoto_Model1_Depression_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Iwamoto_Model1_Depression_FGSEAResults, "Iwamoto_Model1_Depression_FGSEAResults.csv")

Iwamoto_Model2_Depression_FGSEAResults<-fgsea(GMT, Iwamoto_Model2_Depression_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Iwamoto_Model2_Depression_FGSEAResults)
Iwamoto_Model2_Depression_FGSEAResults$leadingEdge<-vapply(Iwamoto_Model2_Depression_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Iwamoto_Model2_Depression_FGSEAResults, "Iwamoto_Model2_Depression_FGSEAResults.csv")

#schizophrenia results
Iwamoto_Model1_Schiz_FGSEAResults<-fgsea(GMT, Iwamoto_Model1_Schiz_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Iwamoto_Model1_Schiz_FGSEAResults)
Iwamoto_Model1_Schiz_FGSEAResults$leadingEdge<-vapply(Iwamoto_Model1_Schiz_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Iwamoto_Model1_Schiz_FGSEAResults, "Iwamoto_Model1_Schizophrenia_FGSEAResults.csv")

Iwamoto_Model2_Schiz_FGSEAResults<-fgsea(GMT, Iwamoto_Model2_Schiz_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Iwamoto_Model2_Schiz_FGSEAResults)
Iwamoto_Model2_Schiz_FGSEAResults$leadingEdge<-vapply(Iwamoto_Model2_Schiz_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Iwamoto_Model2_Schiz_FGSEAResults, "Iwamoto_Model2_Schizophrenia_FGSEAResults.csv")

Maycox_Model1_Schiz_FGSEAResults<-fgsea(GMT, Maycox_Model1_Schiz_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Maycox_Model1_Schiz_FGSEAResults)
Maycox_Model1_Schiz_FGSEAResults$leadingEdge<-vapply(Maycox_Model1_Schiz_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Maycox_Model1_Schiz_FGSEAResults, "Maycox_Model1_Schizophrenia_FGSEAResults.csv")

Maycox_Model2_Schiz_FGSEAResults<-fgsea(GMT, Maycox_Model2_Schiz_log2FC_forGSEA_Ranked, nperm=50000, minSize = 10, maxSize = 1000)
str(Maycox_Model2_Schiz_FGSEAResults)
Maycox_Model2_Schiz_FGSEAResults$leadingEdge<-vapply(Maycox_Model2_Schiz_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))
write.csv(Maycox_Model2_Schiz_FGSEAResults, "Maycox_Model2_Schizophrenia_FGSEAResults.csv")
