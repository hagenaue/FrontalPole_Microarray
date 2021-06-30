##############################

#Limma-voom analysis:

library(limma)

#Change this code to set working directory to where you want to output your differential expression results:
setwd("~/Documents/Microarray Gen/Angela_HRLR_EE_Stress/HC_Remove6samples/JustEESD_CovDissectionCountsPerSample")

#Define the statistical model that you will be using.
#Include each variable that you are interested in (e.g., Diagnosis) as well as important confounds and co-variates (e.g., Age, Hypoxia, etc) as defined by our previous analysis examinig variable-variable correlations and PCA-variable correlations
#data=your metadata data.frame that contains the same subjects as your gene expression dataframe (i.e., after previous outlier removal) in the same order
#In each of the categorical variables should be defined as factors and given appropriate reference levels (e.g., control subjects as the reference for diagnosis, the largest batch or gender group as the reference level for those variables, etc)
design <- model.matrix(~Enrichment+SocialDefeat+date_of_dissection+CountsPerSample, data=Sample_MetaData_HC_NoOutliers_Ordered)

design

#Double check this code - look up in manual (copied from RNA-Seq code)
#Differential expression analysis: Basically this function runs a linear model (defined by design) for all rows of data in your gene expression matrix (all "genes" or transcripts) 
vfit <- lmFit(YourGeneExpressionDataNoOutliers, design)
str(vfit)

#This code applies an empirical Bayes "shrinkage" correction to the logFC/t-stat(? - double-check) results from your differential expression analysis - basically helps protect your results from the influence of extreme outliers/noise
efit <- eBayes(vfit)

#This writes out your results into a tab-delimited text file that can be opened easily in Excel
#Also provides a false detection rate correction for the p-values (FDR or q-value or Benjamini-Hochberg corrected p-value or BH)
#Make file= the file name that fits your statistical model
write.fit(efit, adjust="BH", file="Limma_results_Model_onlyEESD_Dissect_CountsPerSample.txt")

#Writing out the annotation for the rows/genes to go with the results (if you have it...)
write.csv(HC_RNASeq_Annotation_noLowHits, "HC_RNASeq_Annotation_noLowHits.csv")

#Quick glance at results:
#This function tells you how many results you have that meet a particular statistical threshold
#This is nothing special, you can also observe this by sorting your results in Excel
dt<-decideTests(efit)
summary(decideTests(efit))
str(efit)

#Pulls out the results for the top genes associated with your variable of interest
#You need to tell it which variable to attend to - the number is based on the column in your design in your design matrix
#This is nothing special, you can also observe this by sorting your results in Excel
toptable(efit, coef=2)
