#Limma-voom analysis:
#More in-depth annotation is available on the github sample code.

library(limma)

##############################
#   Iwamoto Section Begins   #
##############################

#Transforming the appropriate variables into factors and releveling them
DiagnosisFactor<-relevel(as.factor(Diagnosis), ref="Control")
RateofDeathFactor<-relevel(as.factor(RateofDeath), ref="Sudden")
GenderFactor<-relevel(as.factor(Gender), ref="M")

## Simplest model: ~Diagnosis
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/Simplest_Model")

design <- model.matrix(~DiagnosisFactor, data=Iwamoto_SampleCharacteristics)

design


vfit<- lmFit(SignalSortedNoNA3, design)
str(vfit)

efit <- eBayes(vfit)

#getting the standard error and writing it to a CSV
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
SE<- sqrt(efit$s2.post) * efit$stdev.unscaled
write.csv(SE, "Iwamoto_Model1_StandardErrors.csv")
rm(SE)

write.fit(efit, adjust="BH", file="Limma_results_Model_onlyDiagnosis.txt")
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#number for "coef" determines which variable in design to look at
topTable(efit, coef=2)

## Model controlling for large sources of noise in data (determined by relationships with PC1 & PC2): ~Diagnosis+pH+RNADeg+RateOfDeath
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/Noise_Control_Model")

design <- model.matrix(~DiagnosisFactor+BrainpH+RNADegradPerSample+RateofDeathFactor, data=Iwamoto_SampleCharacteristics)

design


vfit<- lmFit(SignalSortedNoNA3, design)
str(vfit)

efit <- eBayes(vfit)

#getting the standard error and writing it to a CSV
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
SE<- sqrt(efit$s2.post) * efit$stdev.unscaled
write.csv(SE, "Iwamoto_Model2_StandardErrors.csv")
rm(SE)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/Noise_Control_Model")
write.fit(efit, adjust="BH", file="Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath.txt")
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#number for "coef" determines which variable in design to look at
topTable(efit, coef=2)

## Model controlling for large sources of noise and common reviewer demands: ~Diagnosis+pH+RNADeg+RateOfDeath+Age+PMI+Gender
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/Noise_Control_Reviewer_Demand_Model")

design <- model.matrix(~DiagnosisFactor+BrainpH+RNADegradPerSample+RateofDeathFactor+Age+PMI+GenderFactor, data=Iwamoto_SampleCharacteristics)

design


vfit<- lmFit(SignalSortedNoNA3, design)
str(vfit)

efit <- eBayes(vfit)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
SE<- sqrt(efit$s2.post) * efit$stdev.unscaled
write.csv(SE, "Iwamoto_Model3_StandardErrors.csv")
rm(SE)

write.fit(efit, adjust="BH", file="Limma_results_Model_Diagnosis_pH_RNADeg_RateofDeath_Age_PMI_Gender.txt")
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#number for "coef" determines which variable in design to look at
topTable(efit, coef=2)

##############################
#    Maycox Section Begins   #
##############################

#Transforming the appropriate variables into factors and releveling them
DiagnosisFactor<-relevel(as.factor(Diagnosis), ref="Control")

#Gender has one NA value, so I replaced it with "Not Available" in order to allow Limma to function. Ask Dr. Hagenauer if doing this is acceptable.
TempGender<-Gender
TempGender[10]<-"Not Available"
GenderFactor<-relevel(as.factor(TempGender), ref="Male")
rm(TempGender)

## Simplest model: ~Diagnosis+pH
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/Simplest_Model")

design <- model.matrix(~DiagnosisFactor+BrainpH, data=SampleCharacteristics_NoOutliers)

design


vfit<- lmFit(SignalSortedNoNA3, design)
str(vfit)

efit <- eBayes(vfit)

#getting the standard error and writing it to a CSV
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")
SE<- sqrt(efit$s2.post) * efit$stdev.unscaled
write.csv(SE, "Maycox_Model1_StandardErrors.csv")

write.fit(efit, adjust="BH", file="Limma_results_Model_Diagnosis_pH.txt")
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#number for "coef" determines which variable in design to look at
topTable(efit, coef=2)



## Model controlling for large sources of noise in data (determined by relationships with PC1 & PC2): ~Diagnosis+pH+Age+RNADeg+PMI
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/Noise_Control_Model")

design <- model.matrix(~DiagnosisFactor+BrainpH+Age+RNADegradPerSample+PMI, data=SampleCharacteristics_NoOutliers)

design


vfit<- lmFit(SignalSortedNoNA3, design)
str(vfit)

efit <- eBayes(vfit)

#getting the standard error and writing it to a CSV
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")
SE<- sqrt(efit$s2.post) * efit$stdev.unscaled
write.csv(SE, "Maycox_Model2_StandardErrors.csv")


setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/Noise_Control_Model")
write.fit(efit, adjust="BH", file="Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI.txt")
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#number for "coef" determines which variable in design to look at
topTable(efit, coef=2)



## Model controlling for large sources of noise and common reviewer demands: ~Diagnosis+pH+Age+RNADeg+PMI+Gender+ScanDate (???)
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/Noise_Control_Reviewer_Demand_Model")

design <- model.matrix(~DiagnosisFactor+BrainpH+Age+RNADegradPerSample+PMI+GenderFactor+ScanDateDayOnly, data=SampleCharacteristics_NoOutliers)

design


vfit<- lmFit(SignalSortedNoNA3, design)
str(vfit)

efit <- eBayes(vfit)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/Limma_DE_Analysis/All_Models_As_CSV")
SE<- sqrt(efit$s2.post) * efit$stdev.unscaled
write.csv(SE, "Maycox_Model3_StandardErrors.csv")

write.fit(efit, adjust="BH", file="Limma_results_Model_Diagnosis_pH_Age_RNADeg_PMI_Gender_ScanDateDayOnly.txt")
write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#number for "coef" determines which variable in design to look at
topTable(efit, coef=2)