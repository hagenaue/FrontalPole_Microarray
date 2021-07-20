##Sensitivity Analysis


##############################
#   Iwamoto Section Begins   #
##############################
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")

Iwamoto_Model_1<-read.csv("Limma_results_Model_onlyDiagnosis.csv", stringsAsFactors=FALSE)
Iwamoto_Model_2<-read.csv("Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath.csv", stringsAsFactors=FALSE)
Iwamoto_Model_3<-read.csv("Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath_Age_PMI_Gender.csv", stringsAsFactors=FALSE)

Iwamoto_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)

#setting up appropriate column names
colnames(Iwamoto_Model_1)<-paste("Iwamoto_Model1_", colnames(Iwamoto_Model_1), sep="")
colnames(Iwamoto_Model_2)<-paste("Iwamoto_Model2_", colnames(Iwamoto_Model_2), sep="")
colnames(Iwamoto_Model_3)<-paste("Iwamoto_Model3_", colnames(Iwamoto_Model_3), sep="")

#combining models and annotation into one dataframe and saving that dataframe as a csv file
Iwamoto_All_Models_Annotated<-cbind.data.frame(Iwamoto_Model_1, Iwamoto_Model_2, Iwamoto_Model_3, Iwamoto_Annotation)

write.csv(Iwamoto_All_Models_Annotated, "Iwamoto_All_Models_Annotated.csv")

##Visualizing similarities

#Correlation matrix using Log2 fold change (Log2FC)
cor(cbind(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia, Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia, Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia))

#Correlation matrix using t-statistics
cor(cbind(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia, Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia, Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia))

#Scatterplots using Log2FC

#schiz
pdf("Scatterplot_Model1vsModel2_Schiz_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia~Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia, xlab="Model 2 Schiz Log2FC", ylab="Model 1 Schiz Log2FC")
dev.off()

pdf("Scatterplot_Model1vsModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorSchizophrenia~Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, xlab="Model 3 Schiz Log2FC", ylab="Model 1 Schiz Log2FC")
dev.off()

pdf("Scatterplot_Model2vsModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorSchizophrenia~Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, xlab="Model 3 Schiz Log2FC", ylab="Model 2 Schiz Log2FC")
dev.off()

#bp
pdf("Scatterplot_Model1vsModel2_Bipolar_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorBipolar~Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorBipolar, xlab="Model 2 Bipolar Log2FC", ylab="Model 1 Bipolar Log2FC")
dev.off()

pdf("Scatterplot_Model1vsModel3_Bipolar_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorBipolar~Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorBipolar, xlab="Model 3 Bipolar Log2FC", ylab="Model 1 Bipolar Log2FC")
dev.off()

pdf("Scatterplot_Model2vsModel3_Bipolar_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorBipolar~Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorBipolar, xlab="Model 3 Bipolar Log2FC", ylab="Model 2 Bipolar Log2FC")
dev.off()

#depression
pdf("Scatterplot_Model1vsModel2_Depression_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorDepression~Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorDepression, xlab="Model 2 Depression Log2FC", ylab="Model 1 Depression Log2FC")
dev.off()

pdf("Scatterplot_Model1vsModel3_Depression_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_Coef.DiagnosisFactorDepression~Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorDepression, xlab="Model 3 Depression Log2FC", ylab="Model 1 Depression Log2FC")
dev.off()

pdf("Scatterplot_Model2vsModel3_Depression_LogFC.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model2_Coef.DiagnosisFactorDepression~Iwamoto_All_Models_Annotated$Iwamoto_Model3_Coef.DiagnosisFactorDepression, xlab="Model 3 Depression Log2FC", ylab="Model 2 Depression Log2FC")
dev.off()


#Scatterplots using t-statistics

#schiz
pdf("Scatterplot_Model1vsModel2_Schiz_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia~Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia, xlab="Model 2 Schiz t-stat", ylab="Model 1 Schiz t-stat")
dev.off()

pdf("Scatterplot_Model1vsModel3_Schiz_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorSchizophrenia~Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia, xlab="Model 3 Schiz t-stat", ylab="Model 1 Schiz t-stat")
dev.off()

pdf("Scatterplot_Model2vsModel3_Schiz_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorSchizophrenia~Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorSchizophrenia, xlab="Model 3 Schiz t-stat", ylab="Model 2 Schiz t-stat")
dev.off()

#bipolar
pdf("Scatterplot_Model1vsModel2_Bipolar_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorBipolar~Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorBipolar, xlab="Model 2 Bipolar t-stat", ylab="Model 1 Bipolar t-stat")
dev.off()

pdf("Scatterplot_Model1vsModel3_Bipolar_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorBipolar~Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorBipolar, xlab="Model 3 Bipolar t-stat", ylab="Model 1 Bipolar t-stat")
dev.off()

pdf("Scatterplot_Model2vsModel3_Bipolar_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorBipolar~Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorBipolar, xlab="Model 3 Bipolar t-stat", ylab="Model 2 Bipolar t-stat")
dev.off()

#depression
pdf("Scatterplot_Model1vsModel2_Depression_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorDepression~Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorDepression, xlab="Model 2 Depression t-stat", ylab="Model 1 Depression t-stat")
dev.off()

pdf("Scatterplot_Model1vsModel3_Depression_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model1_t.DiagnosisFactorDepression~Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorDepression, xlab="Model 3 Depression t-stat", ylab="Model 1 Depression t-stat")
dev.off()

pdf("Scatterplot_Model2vsModel3_Depression_tstat.pdf", width=5, height=5)
plot(Iwamoto_All_Models_Annotated$Iwamoto_Model2_t.DiagnosisFactorDepression~Iwamoto_All_Models_Annotated$Iwamoto_Model3_t.DiagnosisFactorDepression, xlab="Model 3 Depression t-stat", ylab="Model 2 Depression t-stat")
dev.off()


##############################
#    Maycox Section Begins   #
##############################
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")

Maycox_Model_1<-read.csv("Limma_results_Model_Diagnosis_pH.csv", stringsAsFactors=FALSE)
Maycox_Model_2<-read.csv("Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI.csv", stringsAsFactors=FALSE)
Maycox_Model_3<-read.csv("Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI_Gender_ScanDateDayOnly.csv", stringsAsFactors=FALSE)

Maycox_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)

#setting up appropriate column names
colnames(Maycox_Model_1)<-paste("Maycox_Model1_", colnames(Maycox_Model_1), sep="")
colnames(Maycox_Model_2)<-paste("Maycox_Model2_", colnames(Maycox_Model_2), sep="")
colnames(Maycox_Model_3)<-paste("Maycox_Model3_", colnames(Maycox_Model_3), sep="")

#combining models and annotation into one dataframe and saving that dataframe as a csv file
Maycox_All_Models_Annotated<-cbind.data.frame(Maycox_Model_1, Maycox_Model_2, Maycox_Model_3, Maycox_Annotation)

write.csv(Maycox_All_Models_Annotated, "Maycox_All_Models_Annotated.csv")

##Visualizing similarities

#Correlation matrix using Log2 fold change (Log2FC)
cor(cbind(Maycox_All_Models_Annotated$Maycox_Model1_Coef.DiagnosisFactorScz, Maycox_All_Models_Annotated$Maycox_Model2_Coef.DiagnosisFactorScz, Maycox_All_Models_Annotated$Maycox_Model3_Coef.DiagnosisFactorScz))

#Correlation matrix using t-statistics
cor(cbind(Maycox_All_Models_Annotated$Maycox_Model1_t.DiagnosisFactorScz, Maycox_All_Models_Annotated$Maycox_Model2_t.DiagnosisFactorScz, Maycox_All_Models_Annotated$Maycox_Model3_t.DiagnosisFactorScz))

#Scatterplots using Log2FC
pdf("Scatterplot_Model1vsModel2_Schiz_LogFC.pdf", width=5, height=5)
plot(Maycox_All_Models_Annotated$Maycox_Model1_Coef.DiagnosisFactorScz~Maycox_All_Models_Annotated$Maycox_Model2_Coef.DiagnosisFactorScz, xlab="Model 2 Schiz Log2FC", ylab="Model 1 Schiz Log2FC")
dev.off()

pdf("Scatterplot_Model1vsModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(Maycox_All_Models_Annotated$Maycox_Model1_Coef.DiagnosisFactorScz~Maycox_All_Models_Annotated$Maycox_Model3_Coef.DiagnosisFactorScz, xlab="Model 3 Schiz Log2FC", ylab="Model 1 Schiz Log2FC")
dev.off()

pdf("Scatterplot_Model2vsModel3_Schiz_LogFC.pdf", width=5, height=5)
plot(Maycox_All_Models_Annotated$Maycox_Model2_Coef.DiagnosisFactorScz~Maycox_All_Models_Annotated$Maycox_Model3_Coef.DiagnosisFactorScz, xlab="Model 3 Schiz Log2FC", ylab="Model 2 Schiz Log2FC")
dev.off()

#Scatterplots using t-statistics
pdf("Scatterplot_Model1vsModel2_Schiz_tstat.pdf", width=5, height=5)
plot(Maycox_All_Models_Annotated$Maycox_Model1_t.DiagnosisFactorScz~Maycox_All_Models_Annotated$Maycox_Model2_t.DiagnosisFactorScz, xlab="Model 2 Schiz t-stat", ylab="Model 1 Schiz t-stat")
dev.off()

pdf("Scatterplot_Model1vsModel3_Schiz_tstat.pdf", width=5, height=5)
plot(Maycox_All_Models_Annotated$Maycox_Model1_t.DiagnosisFactorScz~Maycox_All_Models_Annotated$Maycox_Model3_t.DiagnosisFactorScz, xlab="Model 3 Schiz t-stat", ylab="Model 1 Schiz t-stat")
dev.off()

pdf("Scatterplot_Model2vsModel3_Schiz_tstat.pdf", width=5, height=5)
plot(Maycox_All_Models_Annotated$Maycox_Model2_t.DiagnosisFactorScz~Maycox_All_Models_Annotated$Maycox_Model3_t.DiagnosisFactorScz, xlab="Model 3 Schiz t-stat", ylab="Model 2 Schiz t-stat")
dev.off()

