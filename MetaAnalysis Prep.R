#MetaAnalysis Prep


library(plyr)
##############################
#   Iwamoto Section Begins   #
##############################
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
Iwamoto_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)
Iwamoto_Annotation_2<-cbind.data.frame(Iwamoto_Annotation$EntrezGeneID, Iwamoto_Annotation$GeneSymbol)
rm(Iwamoto_Annotation)
Iwamoto_Annotation<-Iwamoto_Annotation_2
rm(Iwamoto_Annotation_2)
colnames(Iwamoto_Annotation)[colnames(Iwamoto_Annotation) == "Iwamoto_Annotation$EntrezGeneID"]<-"EntrezGeneID"
colnames(Iwamoto_Annotation)[colnames(Iwamoto_Annotation) == "Iwamoto_Annotation$GeneSymbol"]<-"GeneSymbol"

Iwamoto_Model2_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath.csv", header=TRUE, stringsAsFactors = FALSE)

IwamotoSchizEffectSizes<-cbind.data.frame(Iwamoto_Model2_NoAnnotation$Coef.DiagnosisFactorSchizophrenia, Iwamoto_Annotation)
IwamotoBipolarEffectSizes<-cbind.data.frame(Iwamoto_Model2_NoAnnotation$Coef.DiagnosisFactorBipolar, Iwamoto_Annotation)

IwamotoStandardErrors<-read.csv("Iwamoto_Model2_StandardErrors.csv", stringsAsFactors=FALSE)
IwamotoSamplingVariances<-IwamotoStandardErrors^2

IwamotoSchizSamplingVariances<-cbind.data.frame(IwamotoSamplingVariances$DiagnosisFactorSchizophrenia, Iwamoto_Annotation)
IwamotoBipolarSamplingVariances<-cbind.data.frame(IwamotoSamplingVariances$DiagnosisFactorBipolar, Iwamoto_Annotation)

colnames(IwamotoSchizEffectSizes)
colnames(IwamotoSchizEffectSizes)[colnames(IwamotoSchizEffectSizes) == "Iwamoto_Model2_NoAnnotation$Coef.DiagnosisFactorSchizophrenia"]<-"Iwamoto_SchizEffectSize"

colnames(IwamotoSchizSamplingVariances)
colnames(IwamotoSchizSamplingVariances)[colnames(IwamotoSchizSamplingVariances) == "IwamotoSamplingVariances$DiagnosisFactorSchizophrenia"]<-"Iwamoto_SchizSamplingVariance"

colnames(IwamotoBipolarEffectSizes)
colnames(IwamotoBipolarEffectSizes)[colnames(IwamotoBipolarEffectSizes) == "Iwamoto_Model2_NoAnnotation$Coef.DiagnosisFactorBipolar"]<-"Iwamoto_BipolarEffectSize"

colnames(IwamotoBipolarSamplingVariances)
colnames(IwamotoBipolarSamplingVariances)[colnames(IwamotoBipolarSamplingVariances) == "IwamotoSamplingVariances$DiagnosisFactorBipolar"]<-"Iwamoto_BipolarSamplingVariance"

IwamotoSchizEffectSizesAndSamplingVariances<-join(IwamotoSchizEffectSizes, IwamotoSchizSamplingVariances, by="EntrezGeneID", type="left", match="all")
IwamotoBipolarEffectSizesAndSamplingVariances<-join(IwamotoBipolarEffectSizes, IwamotoBipolarSamplingVariances, by="EntrezGeneID", type="left", match="all")

#writing out the files for later use
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis")
write.csv(IwamotoSchizEffectSizesAndSamplingVariances, "IwamotoModel2SchizEffectSizesAndSamplingVariances.csv")
write.csv(IwamotoBipolarEffectSizesAndSamplingVariances, "IwamotoModel2BipolarEffectSizesAndSamplingVariances.csv")

##############################
#    Maycox Section Begins   #
##############################
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")
Maycox_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)
Maycox_Annotation_2<-cbind.data.frame(Maycox_Annotation$EntrezGeneID, Maycox_Annotation$GeneSymbol)
rm(Maycox_Annotation)
Maycox_Annotation<-Maycox_Annotation_2
rm(Maycox_Annotation_2)
colnames(Maycox_Annotation)[colnames(Maycox_Annotation) == "Maycox_Annotation$EntrezGeneID"]<-"EntrezGeneID"
colnames(Maycox_Annotation)[colnames(Maycox_Annotation) == "Maycox_Annotation$GeneSymbol"]<-"GeneSymbol"

Maycox_Model2_NoAnnotation<-read.csv("Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI.csv", header=TRUE, stringsAsFactors = FALSE)

MaycoxSchizEffectSizes<-cbind.data.frame(Maycox_Model2_NoAnnotation$Coef.DiagnosisFactorScz, Maycox_Annotation)

MaycoxStandardErrors<-read.csv("Maycox_Model2_StandardErrors.csv", stringsAsFactors=FALSE)
MaycoxSamplingVariances<-MaycoxStandardErrors^2

MaycoxSchizSamplingVariances<-cbind.data.frame(MaycoxSamplingVariances$DiagnosisFactorScz, Maycox_Annotation)

colnames(MaycoxSchizEffectSizes)
colnames(MaycoxSchizEffectSizes)[colnames(MaycoxSchizEffectSizes) == "Maycox_Model2_NoAnnotation$Coef.DiagnosisFactorScz"]<-"Maycox_SchizEffectSize"

colnames(MaycoxSchizSamplingVariances)
colnames(MaycoxSchizSamplingVariances)[colnames(MaycoxSchizSamplingVariances) == "MaycoxSamplingVariances$DiagnosisFactorScz"]<-"Maycox_SchizSamplingVariance"

MaycoxSchizEffectSizesAndSamplingVariances<-join(MaycoxSchizEffectSizes, MaycoxSchizSamplingVariances, by="EntrezGeneID", type="left", match="all")

#writing out the files for later use
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis")
write.csv(MaycoxSchizEffectSizesAndSamplingVariances, "MaycoxModel2SchizEffectSizesAndSamplingVariances.csv")


##############################
#   General Section Begins   #
##############################
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis")

IwamotoVsMaycoxSchizEffectSizes<-join(MaycoxSchizEffectSizes, IwamotoSchizEffectSizes, by="EntrezGeneID", type="left", match="all")
write.csv(IwamotoVsMaycoxSchizEffectSizes, "IwamotoVsMaycoxSchizEffectSizes.csv")
IwamotoVsMaycoxSchizSamplingVariances<-join(MaycoxSchizSamplingVariances, IwamotoSchizSamplingVariances, by="EntrezGeneID", type="left", match="all")
write.csv(IwamotoVsMaycoxSchizSamplingVariances, "IwamotoVsMaycoxSchizSamplingVariances.csv")