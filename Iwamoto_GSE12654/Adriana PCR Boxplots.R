#code for Adriana PCR boxplots


#"GPHN" is left out of this list because it isn't present in the annotation for some reason. See below.
#topGenes<-c("HTR2B","DRD4", "MAOB", "SST", "MAPK1", "ABAT")

potentialtopGenes<-c("GPHN","DRD2", "GAPDH", "HTR6", "MAOA", "REEP5", "TFRC", "AQP4", "CALB1", "GABBR2", "GABRD", "GABRG1", "GAD1", "GFAP", "GJA1", "GNAQ","NSF", "PVALB", "S100B", "SLC1A2", "SLC1A3", "SLC6A1", "SC6A11", "SLC6A13", "SLC38A1", "TBP")

genesNotPresent<-""
topGenes<-""

#checking if genes are present in the annotation
for (i in potentialtopGenes) {
  if (length(which(RMAExpression_customCDFAnnotation2plus2[,3]==i))>0) {
    if (topGenes=="") {
      topGenes<-i
    } else {
      topGenes<-c(topGenes, i)
    }
  } else {
    if (genesNotPresent=="") {
      genesNotPresent<-i
    } else {
      genesNotPresent<-c(genesNotPresent, i)
    }
  }
}

genesNotPresent
#[1] "GPHN"    "GABRG1"  " NSF"    "SC6A11"  "SLC38A1"

topGenes
#[1] "DRD2"    "GAPDH"   "HTR6"    "MAOA"    "REEP5"   "TFRC"    "AQP4"    "CALB1"   "GABBR2"  "GABRD"   "GAD1"    "GFAP"    "GJA1"    "GNAQ"    "PVALB"   "S100B"   "SLC1A2"  "SLC1A3"  "SLC6A1"  "SLC6A13" "TBP"




#Iwamoto Simple Diagnosis boxplots
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Graphs/Adriana_PCR_Boxplots")

for (i in topGenes) {
  
  pdf(paste("Boxplot", i,"ByDiagnosis.pdf", sep="_"), width=8, height=6)
  
  boxplot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]==i),]~Iwamoto_SampleCharacteristics$Diagnosis, xlab="Diagnosis", ylab="", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,5,4), main=i, outline=FALSE)
  stripchart(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]==i),]~Iwamoto_SampleCharacteristics$Diagnosis, vertical = TRUE, 
             method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
  mtext(paste(i, "Log(2) Signal", sep=" "), cex=1.3, side=2, line=2.5)
  dev.off()
  
}
rm(i)

#Iwamoto Diagnosis boxplots with Explanatory variables
LifetimeAlcohol_Heavy<-LifetimeAlcohol=="Heavy in past"|LifetimeAlcohol=="Heavy in present"
LifetimeDrugs_Heavy<-LifetimeDrugs=="Heavy in past"|LifetimeDrugs=="Heavy in present"
SuicideYes<-SuicideStatus=="Yes"
SmokingYes<-SmokingatTimeofDeath=="Yes"

explanatoryVariables<-cbind(LifetimeAlcohol_Heavy, LifetimeDrugs_Heavy, PsychoticFeature, SmokingYes, SuicideYes)
explanatoryVariableNames<-c("LifetimeAlcohol_Heavy", "LifetimeDrugs_Heavy", "PsychoticFeature", "SmokingatTimeofDeath", "SuicideStatus")


theWidth<-40

for (a in 1:length(explanatoryVariables[1,])) {
  
  if (a<=2) {
    theWidth<-10
  } else {
    theWidth<-10
  }
  
  for (b in topGenes) {
    pdf(paste("Boxplot_", b,"_ByDiagnosis_by", explanatoryVariableNames[a], ".pdf", sep=""), width=theWidth, height=6)
    
    boxplot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]==b),]~explanatoryVariables[,a]+Iwamoto_SampleCharacteristics$Diagnosis, xlab="Diagnosis", ylab="", las=1, cex.axis=0.5, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4), main=b, outline=FALSE)
    #axis(1, at=1:40, las = 2, cex.axis = 0.8)
    stripchart(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]==b),]~explanatoryVariables[,a]+Iwamoto_SampleCharacteristics$Diagnosis, vertical = TRUE, 
               method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
    mtext(paste(b, "Log(2) Signal", sep=" "), cex=1.3, side=2, line=2.5)
    dev.off()
  }
}
rm(a)
rm(b)


#Maycox Simple Diagnosis boxplots
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Maycox_GSE17612/RemovedOutliers/Graphs/Adriana_PCR_Boxplots")

for (i in topGenes) {
  
  pdf(paste("Boxplot", i,"ByDiagnosis.pdf", sep="_"), width=4, height=6)
  
  boxplot(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]==i),]~SampleCharacteristics_NoOutliers$Diagnosis, xlab="Diagnosis", ylab="", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,4), main=i, outline=FALSE)
  stripchart(SignalSortedNoNA3[which(RMAExpression_customCDFAnnotation2plus2[,3]==i),]~SampleCharacteristics_NoOutliers$Diagnosis, vertical = TRUE, 
             method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
  mtext(paste(i, "Log(2) Signal", sep=" "), cex=1.3, side=2, line=2.5)
  dev.off()
  
}











