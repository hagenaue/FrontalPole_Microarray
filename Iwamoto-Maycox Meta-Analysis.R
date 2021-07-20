#Iwamoto-Maycox Meta-Analysis
library(metafor)
library(plyr)
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis")

IwamotoVsMaycoxSchizEffectSizes<-read.csv("IwamotoVsMaycoxSchizEffectSizes.csv", header=T, row.names=1, stringsAsFactors = FALSE, na.strings="NA")

#I am removing the EntrezGeneID from this data.frames because otherwise i believe they will be interpretted as numbers
IwamotoVsMaycoxSchizEffectSizes_NoAnnotation<-subset(IwamotoVsMaycoxSchizEffectSizes, select=-EntrezGeneID)

IwamotoVsMaycoxSchizSamplingVariances<-read.csv("IwamotoVsMaycoxSchizSamplingVariances.csv", header=T, row.names=1, stringsAsFactors = FALSE, na.strings="NA")
IwamotoVsMaycoxSchizSamplingVariances_NoAnnotation<-subset(IwamotoVsMaycoxSchizSamplingVariances, select=-EntrezGeneID)

#Making a matrix to store results:
RMAOutputForIwamotoVsMaycoxSchizResults<-matrix(NA, nrow(IwamotoVsMaycoxSchizEffectSizes_NoAnnotation), 3)
colnames(RMAOutputForIwamotoVsMaycoxSchizResults)<-c("b", "se", "pval")

for(i in c(1:nrow(IwamotoVsMaycoxSchizEffectSizes_NoAnnotation))) {
  print(i)
  TempYi<-as.numeric(unname(IwamotoVsMaycoxSchizEffectSizes_NoAnnotation[i,c(1,3)])[1,])
  TempVi<-as.numeric(unname(IwamotoVsMaycoxSchizSamplingVariances_NoAnnotation[i,c(1,3)])[1,])
  if(sum(is.na(TempYi))<1) {
    RMAOutput<-rma(yi=TempYi, vi<-TempVi)
    RMAOutputForIwamotoVsMaycoxSchizResults[i,1]<-RMAOutput$b
    RMAOutputForIwamotoVsMaycoxSchizResults[i,2]<-RMAOutput$se
    RMAOutputForIwamotoVsMaycoxSchizResults[i,3]<-RMAOutput$pval
    print("successful")
  } else{
    print("unsuccessful")
  }
}

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Iwamoto_GSE12654/Limma_DE_Analysis/All_Models_As_CSV")
Iwamoto_Annotation<-read.csv("RMAExpression_customCDFAnnotation2plus2.csv", stringsAsFactors=FALSE)
Iwamoto_Annotation_2<-cbind.data.frame(Iwamoto_Annotation$EntrezGeneID, Iwamoto_Annotation$GeneSymbol)
rm(Iwamoto_Annotation)
Iwamoto_Annotation<-Iwamoto_Annotation_2
rm(Iwamoto_Annotation_2)
colnames(Iwamoto_Annotation)[colnames(Iwamoto_Annotation) == "Iwamoto_Annotation$EntrezGeneID"]<-"EntrezGeneID"
colnames(Iwamoto_Annotation)[colnames(Iwamoto_Annotation) == "Iwamoto_Annotation$GeneSymbol"]<-"GeneSymbol"

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis")

Annotation<-Iwamoto_Annotation
rm(Iwamoto_Annotation)


RMAOutputForIwamotoVsMaycoxSchizResults_noNA<-as.data.frame(cbind(RMAOutputForIwamotoVsMaycoxSchizResults, IwamotoVsMaycoxSchizEffectSizes$EntrezGeneID))
#renaming collumns
colnames(RMAOutputForIwamotoVsMaycoxSchizResults_noNA)[4]<-"EntrezGeneID"

#removing NAs
RMAOutputForIwamotoVsMaycoxSchizResults_noNA<-na.omit(RMAOutputForIwamotoVsMaycoxSchizResults_noNA)

RMAOutputForIwamotoVsMaycoxSchizResults<-na.omit(RMAOutputForIwamotoVsMaycoxSchizResults)

Annotated_RMAOutputForIwamotoVsMaycoxSchizResults_noNA<-join(RMAOutputForIwamotoVsMaycoxSchizResults_noNA, Annotation, by="EntrezGeneID", type="inner", match="all")

#Running FDR corrections on the p-values:
library(multtest)
TempPvalAdj<-mt.rawp2adjp(RMAOutputForIwamotoVsMaycoxSchizResults[,3], proc=c("BH"))
RMAOutputForIwamotoVsMaycoxSchizResults<-cbind(RMAOutputForIwamotoVsMaycoxSchizResults, TempPvalAdj$adjp[order(TempPvalAdj$index),2])
colnames(RMAOutputForIwamotoVsMaycoxSchizResults)[4]<-"BH_adjPval"



RMAOutputForIwamotoVsMaycoxSchizResults<-cbind(RMAOutputForIwamotoVsMaycoxSchizResults, Annotated_RMAOutputForIwamotoVsMaycoxSchizResults_noNA$EntrezGeneID, Annotated_RMAOutputForIwamotoVsMaycoxSchizResults_noNA$GeneSymbol)
colnames(RMAOutputForIwamotoVsMaycoxSchizResults)[5]<-"EntrezGeneID"
colnames(RMAOutputForIwamotoVsMaycoxSchizResults)[6]<-"GeneSymbol"

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis/Output")
write.csv(RMAOutputForIwamotoVsMaycoxSchizResults, "RMAOutputForIwamotoVsMaycoxSchizResults.csv")

