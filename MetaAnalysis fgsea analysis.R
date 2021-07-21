#meta-analysis FGSEA analysis
library(fgsea)
library(plyr)

#reading in metaanalysis results
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis/Output")
IwaCox_SchizMetaResults<-read.csv("RMAOutputForIwamotoVsMaycoxSchizResults.csv", stringsAsFactors = FALSE)
IwaCox_SchizMetaResults<-subset(IwaCox_SchizMetaResults, select=-X)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/fGSEA")

MetaResults_noNA<-IwaCox_SchizMetaResults[is.na(IwaCox_SchizMetaResults$GeneSymbol)==FALSE,]

sum(duplicated(MetaResults_noNA$GeneSymbol))
# [1] 0

MetaResults_Betas_forGSEA<-tapply(X=MetaResults_noNA$b, INDEX=MetaResults_noNA$GeneSymbol, FUN=mean)
names(MetaResults_Betas_forGSEA)<-names(table(MetaResults_noNA$GeneSymbol))

length(MetaResults_Betas_forGSEA)
# [1] 8357

MetaResults_Betas_forGSEARanked<-MetaResults_Betas_forGSEA[order(MetaResults_Betas_forGSEA)]
head(MetaResults_Betas_forGSEARanked)

#Reading in updated GMT
GMT<-gmtPathways("c5.go.v7.3.symbols_BrainInABlender_Psych.gmt.txt")
str(GMT)

setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/fGSEA/MetaAnalysis")

MetaAnalysis_Schiz_FGSEAResults<-fgsea(GMT, MetaResults_Betas_forGSEARanked, nperm=50000, minSize = 10, maxSize = 1000)

#I keep getting the following error when trying to run the above line of code:
# Error in serialize(data, node$con) : error writing to connection
# In addition: Warning message:
#   In serialize(data, node$con) :
#   'package:stats' may not be available when loading


str(MetaAnalysis_Schiz_FGSEAResults)

MetaAnalysis_Schiz_FGSEAResults$leadingEdge<-vapply(MetaAnalysis_Schiz_FGSEAResults$leadingEdge, paste, collapse= ",", character(1L))

write.csv(MetaAnalysis_Schiz_FGSEAResults, "IwamotoVsMaycoxMetaAnalysis_SchizFGSEAResults.csv")

