#Iwamoto Vs Maycox Meta Analysis Vs Adriana
library(plyr)

#Reading in non-meta analysis results
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox")

AdrianaResults_vs_IwamotoAndMaycoxResults_noMetaAnalysisResults<-read.csv("AdrianaResults_vs_IwamotoAndMaycox.csv", stringsAsFactors = FALSE)

#reading in the meta analysis results
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwaCox_MetaAnalysis/Output")
IwaCox_SchizMetaResults<-read.csv("RMAOutputForIwamotoVsMaycoxSchizResults.csv", stringsAsFactors = FALSE)
IwaCox_SchizMetaResults<-subset(IwaCox_SchizMetaResults, select=-X)
colnames(IwaCox_SchizMetaResults)<-paste("IwamotoVsMaycoxMetaResults_", colnames(IwaCox_SchizMetaResults), sep="")
colnames(IwaCox_SchizMetaResults)[5]<-"EntrezGeneID"
colnames(IwaCox_SchizMetaResults)[6]<-"GeneSymbol"

#Joining the datasets
setwd("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/IwamotoVsMaycoxMetaAnalysisVsAdriana")

AllDataSets<-join(AdrianaResults_vs_IwamotoAndMaycoxResults_noMetaAnalysisResults, IwaCox_SchizMetaResults, by="GeneSymbol", type="inner", match="all")
write.csv(AllDataSets, "AdrianaResults_Plus_IwamotoAndMaycoxResults_Plus_MetaAnalysisResults.csv")


#Schiz Scatterplots using Log2FC from Adriana's data and effect size estimate from the IwamotoVsMaycox meta-analysis results
pdf("Scatterplot_IwamotoMaycoxMetaAnalysisvsAdriana_Schiz_LogFC.pdf", width=5, height=5)
plot(AllDataSets$IwamotoVsMaycoxMetaResults_b~AllDataSets$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta, xlab="Adriana Schiz Log2FC", ylab="IwamotoVsMaycox Meta-Analysis Schiz Effect Size")
TempLine<-lm(AllDataSets$IwamotoVsMaycoxMetaResults_b~AllDataSets$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta)
abline(TempLine, col=2, lwd=2)
dev.off()

summary.lm(TempLine)

# Call:
# lm(formula = AllDataSets$IwamotoVsMaycoxMetaResults_b ~ AllDataSets$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.19194 -0.05946 -0.01010  0.04619  0.47910 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                                  0.01314    0.01060   1.239   0.2184  
# AllDataSets$Adriana_Diagnosis_Schiz_PostHocSummary_MLM_Beta  0.16026    0.06625   2.419   0.0175 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.103 on 93 degrees of freedom
# Multiple R-squared:  0.0592,	Adjusted R-squared:  0.04908 
# F-statistic: 5.852 on 1 and 93 DF,  p-value: 0.01751