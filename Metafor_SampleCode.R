#Determining which of the SubjVarVsCellType effects are replicating across the four datasets (Pritzker, Narayan, Lanz, Barnes) using a make-shift little meta-analysis:
#Megan Hagenauer 10-9-2017

library(metafor)
setwd("~/Documents/Microarray Gen/Cell Type Paper/Replication")

#Betas = log2FC or coefficient from limma model
CellTypeVsSubjVarBetas<-read.csv("CellTypeVsSubjVarBetas_ForMetafor.csv", header=T, row.names=1, stringsAsFactors = F, na.strings="NA")
str(CellTypeVsSubjVarBetas)
# 'data.frame':	90 obs. of  4 variables:
#   $ Pritzker.DLPFC.     : num  -0.43655 0.071 -0.00399 0.00201 0.05953 ...
# $ Barnes.et.al..2011. : num  -0.3676 NA 0.000647 0.001487 0.2372 ...
# $ Narayan.et.al..2008.: num  -0.856764 NA -0.005818 -0.000648 0.278124 ...
# $ Lanz.2014           : num  -0.46067 NA 0.00212 0.00949 0.11853 ...

#Variances = standard error from limma model (SE) squared
CellTypeVsSubjVarVariances<-read.csv("CellTypeVsSubjVarVariances_ForMetafor.csv", header=T, row.names=1, stringsAsFactors = F, na.strings="NA")
str(CellTypeVsSubjVarVariances)
# 'data.frame':	90 obs. of  4 variables:
#   $ Pritzker.DLPFC.     : num  6.46e-03 9.73e-04 5.14e-06 2.13e-06 2.30e-03 ...
# $ Barnes.et.al..2011. : num  1.12e-01 NA 1.99e-04 1.55e-05 1.86e-02 ...
# $ Narayan.et.al..2008.: num  9.01e-02 NA 2.15e-05 1.28e-05 2.31e-02 ...
# $ Lanz.2014           : num  2.82e-02 NA 6.13e-05 2.08e-05 8.09e-03 ...


#Example run:

#So I think this is probably unhappy because my data has names. :(
TempYi<-as.numeric(unname(CellTypeVsSubjVarBetas[1,])[1,])
TempV<-as.numeric(unname(CellTypeVsSubjVarVariances[1,])[1,])
rma.mv(yi=TempYi, V=TempV)


# Multivariate Meta-Analysis Model (k = 4; method: REML)
# 
# Variance Components: none
# 
# Test for Heterogeneity: 
#   Q(df = 3) = 1.9085, p-val = 0.5916
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     
# -0.4599  0.0690  -6.6682  <.0001  -0.5951  -0.3247  ***
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


RMAOutput<-rma.mv(yi=TempYi, V=TempV)
forest(RMAOutput)

RMAOutput$b
#intrcpt -0.4598921
RMAOutput$se
#[1] 0.06896819

RMAOutput$pval
#[1] 2.590011e-11

RMAOutputForAllCellTypeVsVar<-matrix(NA, nrow(CellTypeVsSubjVarBetas), 3)
colnames(RMAOutputForAllCellTypeVsVar)<-c("b", "se", "pval")
row.names(RMAOutputForAllCellTypeVsVar)<-row.names(CellTypeVsSubjVarBetas)

#We should change the if/then statement to require at least 2 datasets (or >1) to have data to run the meta-analysis.
for(i in c(1:nrow(CellTypeVsSubjVarBetas))){
  print(i)
  TempYi<-as.numeric(unname(CellTypeVsSubjVarBetas[i,])[1,])
  TempV<-as.numeric(unname(CellTypeVsSubjVarVariances[i,])[1,])
  if(sum(is.na(TempYi))<3){
  RMAOutput<-rma.mv(yi=TempYi, V=TempV)
  RMAOutputForAllCellTypeVsVar[i,1]<-RMAOutput$b
  RMAOutputForAllCellTypeVsVar[i,2]<-RMAOutput$se
  RMAOutputForAllCellTypeVsVar[i,3]<-RMAOutput$pval
  }else{}
}



library(multtest)
  TempPvalAdj<-mt.rawp2adjp(RMAOutputForAllCellTypeVsVar[,3], proc=c("BH"))
  
  RMAOutputForAllCellTypeVsVar<-cbind(RMAOutputForAllCellTypeVsVar, TempPvalAdj$adjp[order(TempPvalAdj$index),2])
  colnames(RMAOutputForAllCellTypeVsVar)[4]<-"BH_adjPval"
  
  write.csv(RMAOutputForAllCellTypeVsVar, "RMAOutputForAllCellTypeVsVar.csv")
  
