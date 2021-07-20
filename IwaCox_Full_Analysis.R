##Iwamoto vs. Maycox: how similar are the results?


#reading in data
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
Iwamoto_Results<-read.csv("Iwamoto_All_Models_Annotated.csv", stringsAsFactors=FALSE)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")
Maycox_Results<-read.csv("Maycox_All_Models_Annotated.csv", stringsAsFactors=FALSE)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox")

#removing NAs
Iwamoto_Results_NoNA<-Iwamoto_Results[is.na(Iwamoto_Results$EntrezGeneID)==FALSE,]
Maycox_Results_NoNA<-Maycox_Results[is.na(Maycox_Results$EntrezGeneID)==FALSE,]

#joining the data
library(plyr)
IwamotoVsMaycox_AllResults<-join(Iwamoto_Results_NoNA, Maycox_Results_NoNA, by="EntrezGeneID", type="full", match="all")
write.csv(IwamotoVsMaycox_AllResults, "IwamotoVsMaycox_AllResults.csv")

#correlation matrix
cor(cbind(IwamotoVsMaycox_AllResults$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia, IwamotoVsMaycox_AllResults$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia, IwamotoVsMaycox_AllResults$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, IwamotoVsMaycox_AllResults$Maycox_Model1_Coef.DiagnosisFactorScz, IwamotoVsMaycox_AllResults$Maycox_Model2_Coef.DiagnosisFactorScz, IwamotoVsMaycox_AllResults$Maycox_Model3_Coef.DiagnosisFactorScz), use="pairwise.complete.obs")

#Schiz Scatterplots using Log2FC

pdf("Scatterplot_IwamotoModel1vsMaycoxModel1_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model1_Coef.DiagnosisFactorScz, xlab="Model 1 Maycox Schiz Log2FC", ylab="Model 1 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel1vsMaycoxModel2_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model2_Coef.DiagnosisFactorScz, xlab="Model 2 Maycox Schiz Log2FC", ylab="Model 1 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel1vsMaycoxModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model3_Coef.DiagnosisFactorScz, xlab="Model 3 Maycox Schiz Log2FC", ylab="Model 1 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel2vsMaycoxModel1_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model1_Coef.DiagnosisFactorScz, xlab="Model 1 Maycox Schiz Log2FC", ylab="Model 2 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel2vsMaycoxModel2_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model2_Coef.DiagnosisFactorScz, xlab="Model 2 Maycox Schiz Log2FC", ylab="Model 2 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel2vsMaycoxModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model3_Coef.DiagnosisFactorScz, xlab="Model 3 Maycox Schiz Log2FC", ylab="Model 2 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel3vsMaycoxModel1_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model1_Coef.DiagnosisFactorScz, xlab="Model 1 Maycox Schiz Log2FC", ylab="Model 3 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel3vsMaycoxModel2_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model2_Coef.DiagnosisFactorScz, xlab="Model 2 Maycox Schiz Log2FC", ylab="Model 3 Iwamoto Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel3vsMaycoxModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model3_Coef.DiagnosisFactorScz, xlab="Model 3 Maycox Schiz Log2FC", ylab="Model 3 Iwamoto Schiz Log2FC")
dev.off()

#Schiz Scatterplots using t-statistics

pdf("Scatterplot_IwamotoModel1vsMaycoxModel1_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model1_t.DiagnosisFactorScz, xlab="Model 1 Maycox Schiz T-Stat", ylab="Model 1 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel1vsMaycoxModel2_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model2_t.DiagnosisFactorScz, xlab="Model 2 Maycox Schiz T-Stat", ylab="Model 1 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel1vsMaycoxModel3_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model3_t.DiagnosisFactorScz, xlab="Model 3 Maycox Schiz T-Stat", ylab="Model 1 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel2vsMaycoxModel1_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model1_t.DiagnosisFactorScz, xlab="Model 1 Maycox Schiz T-Stat", ylab="Model 2 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel2vsMaycoxModel2_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model2_t.DiagnosisFactorScz, xlab="Model 2 Maycox Schiz T-Stat", ylab="Model 2 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel2vsMaycoxModel3_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model3_t.DiagnosisFactorScz, xlab="Model 3 Maycox Schiz T-Stat", ylab="Model 2 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel3vsMaycoxModel1_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model1_t.DiagnosisFactorScz, xlab="Model 1 Maycox Schiz T-Stat", ylab="Model 3 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel3vsMaycoxModel2_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model2_t.DiagnosisFactorScz, xlab="Model 2 Maycox Schiz T-Stat", ylab="Model 3 Iwamoto Schiz T-Stat")
dev.off()

pdf("Scatterplot_IwamotoModel3vsMaycoxModel3_Schiz_tstat.pdf", width=5, height=5)
plot(IwamotoVsMaycox_AllResults$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia~IwamotoVsMaycox_AllResults$Maycox_Model3_t.DiagnosisFactorScz, xlab="Model 3 Maycox Schiz T-Stat", ylab="Model 3 Iwamoto Schiz T-Stat")
dev.off()