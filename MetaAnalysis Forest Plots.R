#Making forest plots for top meta-analysis genes
library(metafor)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis/Output")
IwaCox_SchizMetaResults<-read.csv("RMAOutputForIwamotoVsMaycoxSchizResults.csv", stringsAsFactors = FALSE)
IwaCox_SchizMetaResults<-subset(IwaCox_SchizMetaResults, select=-X)
colnames(IwaCox_SchizMetaResults)<-paste("IwamotoVsMaycoxMetaResults_", colnames(IwaCox_SchizMetaResults), sep="")
colnames(IwaCox_SchizMetaResults)[5]<-"EntrezGeneID"
colnames(IwaCox_SchizMetaResults)[6]<-"GeneSymbol"

meta<-IwaCox_SchizMetaResults
rm(IwaCox_SchizMetaResults)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/Forest_Plots")



#checking for Adriana's top genes in the meta-analysis results (HTR2B, DRD4, SST, MAPK1, ABAT)

temp<-subset(meta, meta$GeneSymbol=="HTR2B" | meta$GeneSymbol=="DRD4" | meta$GeneSymbol=="SST" | meta$GeneSymbol=="MAPK1" | meta$GeneSymbol=="ABAT")
sum(meta$GeneSymbol=="HTR2B")
