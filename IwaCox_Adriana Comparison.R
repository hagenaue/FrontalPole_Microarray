#Comparing IwaCox results with Adriana's results

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox")

Adriana_Results<-read.csv("TableS4_BothDatasets_ResultSummary_ForAllTargetGenes_ForLiam.csv", stringsAsFactors = FALSE)

#naming the columns
colnames(Adriana_Results)
colnames(Adriana_Results)<-paste("Adriana_", colnames(Adriana_Results), sep="")
colnames(Adriana_Results)[colnames(Adriana_Results) == "Adriana_ï..Gene.Symbol"] <- "GeneSymbol"

#joining the datasets and writing out the resulting combined file
library(plyr)
AdrianaResults_vs_IwamotoAndMaycox<-join(Adriana_Results, IwamotoVsMaycox_AllResults, by="GeneSymbol", type="left", match="all")
write.csv(AdrianaResults_vs_IwamotoAndMaycox, "AdrianaResults_vs_IwamotoAndMaycox.csv")

#Schiz correlation matrix
cor(cbind(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia, AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia, AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model1_Coef.DiagnosisFactorScz, AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model2_Coef.DiagnosisFactorScz, AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model3_Coef.DiagnosisFactorScz, AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta), use="pairwise.complete.obs")


#Schiz Scatterplots using Log2FC
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox/Adriana_Results/Schiz")

pdf("Scatterplot_IwamotoModel1vsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="Iwamoto Model 1 Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel2vsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="Iwamoto Model 2 Schiz Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel3vsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="Iwamoto Model 3 Schiz Log2FC")
dev.off()

pdf("Scatterplot_MaycoxModel1vsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model1_Coef.DiagnosisFactorScz~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="Maycox Model 1 Schiz Log2FC")
dev.off()

pdf("Scatterplot_MaycoxModel2vsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model2_Coef.DiagnosisFactorScz~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="Maycox Model 2 Schiz Log2FC")
dev.off()

pdf("Scatterplot_MaycoxModel3vsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model3_Coef.DiagnosisFactorScz~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="Maycox Model 3 Schiz Log2FC")
dev.off()

#Schiz Scatterplots using t-statistics
pdf("Scatterplot_IwamotoModel1vsAdriana_Schiz_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Tstat, xlab="Adriana Schiz t-stat", ylab="Iwamoto Model 1 Schiz t-stat")
dev.off()

pdf("Scatterplot_IwamotoModel2vsAdriana_Schiz_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Tstat, xlab="Adriana Schiz t-stat", ylab="Iwamoto Model 2 Schiz t-stat")
dev.off()

pdf("Scatterplot_IwamotoModel3vsAdriana_Schiz_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Tstat, xlab="Adriana Schiz t-stat", ylab="Iwamoto Model 3 Schiz t-stat")
dev.off()

pdf("Scatterplot_MaycoxModel1vsAdriana_Schiz_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model1_t.DiagnosisFactorScz~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Tstat, xlab="Adriana Schiz t-stat", ylab="Maycox Model 1 Schiz t-stat")
dev.off()

pdf("Scatterplot_MaycoxModel2vsAdriana_Schiz_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model2_t.DiagnosisFactorScz~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Tstat, xlab="Adriana Schiz t-stat", ylab="Maycox Model 2 Schiz t-stat")
dev.off()

pdf("Scatterplot_MaycoxModel3vsAdriana_Schiz_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Maycox_Model3_t.DiagnosisFactorScz~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Tstat, xlab="Adriana Schiz t-stat", ylab="Maycox Model 3 Schiz t-stat")
dev.off()

#BP Scatterplots using Log2FC
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox/Adriana_Results/BP")

pdf("Scatterplot_IwamotoModel1vsAdriana_BP_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model1_Coef.DiagnosisFactorBipolar~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_BP_PostHocSummary_MLM_Beta, xlab="Adriana BP Log2FC", ylab="Iwamoto Model 1 BP Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel2vsAdriana_BP_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model2_Coef.DiagnosisFactorBipolar~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_BP_PostHocSummary_MLM_Beta, xlab="Adriana BP Log2FC", ylab="Iwamoto Model 2 BP Log2FC")
dev.off()

pdf("Scatterplot_IwamotoModel3vsAdriana_BP_LogFC.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model3_Coef.DiagnosisFactorBipolar~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_BP_PostHocSummary_MLM_Beta, xlab="Adriana BP Log2FC", ylab="Iwamoto Model 3 BP Log2FC")
dev.off()

#BP Scatterplots using t-statistics
pdf("Scatterplot_IwamotoModel1vsAdriana_BP_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model1_t.DiagnosisFactorBipolar~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_BP_PostHocSummary_MLM_Tstat, xlab="Adriana BP t-stat", ylab="Iwamoto Model 1 BP t-stat")
dev.off()

pdf("Scatterplot_IwamotoModel2vsAdriana_BP_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model2_t.DiagnosisFactorBipolar~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_BP_PostHocSummary_MLM_Tstat, xlab="Adriana BP t-stat", ylab="Iwamoto Model 2 BP t-stat")
dev.off()

pdf("Scatterplot_IwamotoModel3vsAdriana_BP_t-stat.pdf", width=5, height=5)
plot(AdrianaResults_vs_IwamotoAndMaycox$Iwamoto_Model3_t.DiagnosisFactorBipolar~AdrianaResults_vs_IwamotoAndMaycox$Adriana_Diagnosis_BP_PostHocSummary_MLM_Tstat, xlab="Adriana BP t-stat", ylab="Iwamoto Model 3 BP t-stat")
dev.off()

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox")
