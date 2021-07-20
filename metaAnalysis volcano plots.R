#Meta-analysis Volcano Plots


#reading in metadata dataframe and removing NA values
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis/Output")
metaAnalysisOutputFDR<-read.csv("RMAOutputForIwamotoVsMaycoxSchizResults.csv", stringsAsFactors = FALSE)
metaAnalysisOutputFDR_noNA<-metaAnalysisOutputFDR[is.na(metaAnalysisOutputFDR$b)==FALSE,]
temp<-metaAnalysisOutputFDR_noNA
metaAnalysisOutputFDR_noNA<-temp[is.na(temp$se)==FALSE,]
temp<-metaAnalysisOutputFDR_noNA
metaAnalysisOutputFDR_noNA<-temp[is.na(temp$pval)==FALSE,]
temp<-metaAnalysisOutputFDR_noNA
metaAnalysisOutputFDR_noNA<-temp[is.na(temp$BH_adjPval)==FALSE,]
rm(temp)

#reading in Iwamoto limma results
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
Iwamoto_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)
Iwamoto_Annotation_2<-cbind.data.frame(Iwamoto_Annotation$EntrezGeneID, Iwamoto_Annotation$GeneSymbol)
rm(Iwamoto_Annotation)
Iwamoto_Annotation<-Iwamoto_Annotation_2
rm(Iwamoto_Annotation_2)
colnames(Iwamoto_Annotation)[colnames(Iwamoto_Annotation) == "Iwamoto_Annotation$EntrezGeneID"]<-"EntrezGeneID"
colnames(Iwamoto_Annotation)[colnames(Iwamoto_Annotation) == "Iwamoto_Annotation$GeneSymbol"]<-"GeneSymbol"

Iwamoto_Model2_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath.csv", header=TRUE, stringsAsFactors = FALSE)
Iwamoto_Model2<-cbind.data.frame(Iwamoto_Model2_NoAnnotation, Iwamoto_Annotation)
rm(Iwamoto_Annotation)
rm(Iwamoto_Model2_NoAnnotation)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Volcano_Plots")
outDir<-"C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Volcano_Plots/Plots"

#setting up the plot variables for Schiz
estimate<-Iwamoto_Model2$Coef.DiagnosisFactorSchizophrenia
pval<-Iwamoto_Model2$p.value.DiagnosisFactorSchizophrenia
padj<-Iwamoto_Model2$p.value.adj.DiagnosisFactorSchizophrenia

#making the plot for Schiz
tiff(paste0(outDir, "VolcanoPlotIwamotoSchiz.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))


#am I supposed to be using the meta-analysis "metaAnalysisOutputFDR" for the individual study plots as well?
par(mai=c(1.02, 1,0.9,0.40))
with(Iwamoto_Model2, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                 xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1) for Schiz
with(subset(Iwamoto_Model2, abs(estimate)>1), points(estimate, -log10(pval), 
                                                            pch=19, col="red", cex=0.6))
with(subset(Iwamoto_Model2, padj< .05 ), points(estimate, -log10(pval), 
                                                     pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()

##BP

#setting up the plot variables for BP
estimate<-Iwamoto_Model2$Coef.DiagnosisFactorBipolar
pval<-Iwamoto_Model2$p.value.DiagnosisFactorBipolar
padj<-Iwamoto_Model2$p.value.adj.DiagnosisFactorBipolar

#making the plot for BP
tiff(paste0(outDir, "VolcanoPlotIwamotoBipolarDisorder.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))


#am I supposed to be using the meta-analysis "metaAnalysisOutputFDR" for the individual study plots as well?
par(mai=c(1.02, 1,0.9,0.40))
with(Iwamoto_Model2, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                          xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1) for BP
with(subset(Iwamoto_Model2, abs(estimate)>1), points(estimate, -log10(pval), 
                                                     pch=19, col="red", cex=0.6))
with(subset(Iwamoto_Model2, padj< .1 ), points(estimate, -log10(pval), 
                                                pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()


#Maycox
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")
Maycox_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)
Maycox_Annotation_2<-cbind.data.frame(Maycox_Annotation$EntrezGeneID, Maycox_Annotation$GeneSymbol)
rm(Maycox_Annotation)
Maycox_Annotation<-Maycox_Annotation_2
rm(Maycox_Annotation_2)
colnames(Maycox_Annotation)[colnames(Maycox_Annotation) == "Maycox_Annotation$EntrezGeneID"]<-"EntrezGeneID"
colnames(Maycox_Annotation)[colnames(Maycox_Annotation) == "Maycox_Annotation$GeneSymbol"]<-"GeneSymbol"

Maycox_Model2_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI.csv", header=TRUE, stringsAsFactors = FALSE)
Maycox_Model2<-cbind.data.frame(Maycox_Model2_NoAnnotation, Maycox_Annotation)
rm(Maycox_Annotation)
rm(Maycox_Model2_NoAnnotation)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Volcano_Plots")

#setting up the plot variables
estimate<-Maycox_Model2$Coef.DiagnosisFactorScz
pval<-Maycox_Model2$p.value.DiagnosisFactorScz
padj<-Maycox_Model2$p.value.adj.DiagnosisFactorScz

#making the plot
tiff(paste0(outDir, "VolcanoPlotMaycoxSchiz.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))


#am I supposed to be using the meta-analysis "metaAnalysisOutputFDR" for the individual study plots as well?
par(mai=c(1.02, 1,0.9,0.40))
with(Maycox_Model2, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                          xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(Maycox_Model2, abs(estimate)>1), points(estimate, -log10(pval), 
                                                     pch=19, col="red", cex=0.6))
with(subset(Maycox_Model2, padj< .1 ), points(estimate, -log10(pval), 
                                                pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()


#IwamotoVsMaycox Meta Analysis

#setting up the plot variables
estimate<-metaAnalysisOutputFDR_noNA$b
pval<-metaAnalysisOutputFDR_noNA$pval
padj<-metaAnalysisOutputFDR_noNA$BH_adjPval

#making the plot
tiff(paste0(outDir, "VolcanoPlotIwamotoVsMaycoxSchiz.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))


par(mai=c(1.02, 1,0.9,0.40))
with(metaAnalysisOutputFDR_noNA, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                         xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1)
# NONE OF THE VALUES FIT WITHIN THE abs(estimate)>1) or the padj< .05.
with(subset(metaAnalysisOutputFDR_noNA, abs(estimate)>1), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
with(subset(metaAnalysisOutputFDR_noNA, padj< .1 ), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))

legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()
