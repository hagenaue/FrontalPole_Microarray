#Lanz data download and analysis (GSE53987)
#Megan Hagenauer, 1/2017-7/2017
#Affymetrix Human Genome U133 Plus 2.0 Array

setwd("~/Documents/Microarray Gen/Lanz_GSE53987")

library(GEOquery)
gse <- getGEO("GSE53987", GSEMatrix = FALSE)
#lots of warnings. I wonder why.

head(Meta(gse))
str(Meta(gse))
Meta(GSMList(gse)$GSM1304852)$characteristics_ch1
# [1] "age: 52"                         "gender: M"                       "race: W"                        
# [4] "pmi: 23.5"                       "ph: 6.7"                         "rin: 6.3"                       
# [7] "tissue: hippocampus"             "disease state: bipolar disorder"


sub("gender: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[2])
as.numeric(sub("age: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[1]))
as.numeric(sub("pmi: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[4], fixed = T))
as.numeric(sub("ph: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[5]))
sub("disease state: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[8], fixed=T)
sub("tissue: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[7])
as.numeric(sub("rin: ","", Meta(GSMList(gse)$GSM1304852)$characteristics_ch1[6]))

SampleID<-as.matrix(names(GSMList(gse)))
Gender<-matrix("a", nrow=205, ncol=1)
Age<-matrix(0, nrow=205, ncol=1)
BrainpH<-matrix(0, nrow=205, ncol=1)
PMI<-matrix(0, nrow=205, ncol=1)
Diagnosis2<-matrix("a", nrow=205, ncol=1)
Tissue<-matrix("a", nrow=205, ncol=1)
RIN<-matrix(0, nrow=205, ncol=1)

length(GSMList(gse))

i<-1
for(i in c(1:205)){
  Gender[i]<-sub("gender: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[2])
  Age[i]<-as.numeric(sub("age: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[1]))
  PMI[i]<-as.numeric(sub("pmi: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[4], fixed = T))
  BrainpH[i]<-as.numeric(sub("ph: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[5]))
  Diagnosis2[i]<-sub("disease state: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[8], fixed=T)
  Tissue[i]<-sub("tissue: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[7])
  RIN[i]<-as.numeric(sub("rin: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[6]))
}
Gender<-relevel(as.factor(Gender), ref="M")

Diagnosis<-Diagnosis2
Diagnosis<-relevel(as.factor(Diagnosis), ref="control")

Lanz_SampleCharacteristics<-data.frame(SampleID, Gender, Age, BrainpH, PMI, Diagnosis, Tissue, RIN, stringsAsFactors=F)

head(Lanz_SampleCharacteristics)

#Looks good.

write.csv(Lanz_SampleCharacteristics, "Lanz_SampleCharacteristics.csv")

setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC")


#Just for the DLPFC-Only Analyses:
temp<-Lanz_SampleCharacteristics[Tissue=="Pre-frontal cortex (BA46)",]
write.csv(temp, "Lanz_SampleCharacteristics_DLPFC.csv")

Lanz_SampleCharacteristics<-temp
SampleID<-temp$SampleID
Gender<-temp$Gender
Age<-temp$Age
PMI<-temp$PMI
BrainpH<-temp$BrainpH
Diagnosis2<-temp$Diagnosis
Tissue<-temp$Tissue
RIN<-temp$RIN


rm(gse)


#Now let's re-run RMA on their data:

library(org.Hs.eg.db)
library(plyr)
library(affy)

#This is where I obtained the updated custom .cdf for defining the probesets:
http://nmg-r.bioinformatics.nl/NuGO_R.html

#cdf and the chip.

# hgu133plus2hsentrezg.db_19.0.2.tar.gz
# hgu133plus2hsentrezgcdf_19.0.0.tar.gz
# hgu133plus2hsentrezgprobe_19.0.0.tar.gz

#Already installed:
#install.packages(pkgs = c("hgu133plus2hsentrezg.db", "hgu133plus2hsentrezgcdf", "hgu133plus2hsentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

#Changed working directory to where the .cel files are located
#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/RawData/GSE53987_RAW")

#For DLPFC-only analyses:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/RawData/JustDLPFC")

data2<-ReadAffy(cdfname ="hgu133plus2hsentrezg")
str(data2)
data2
# size of arrays=1164x1164 features (105 kb)
# cdf=hgu133plus2hsentrezg (19764 affyids)
# number of samples=205
# number of genes=19764
# annotation=hgu133plus2hsentrezg
# notes=

#Just DLPFC
# AffyBatch object
# size of arrays=1164x1164 features (44 kb)
# cdf=hgu133plus2hsentrezg (19764 affyids)
# number of samples=68
# number of genes=19764
# annotation=hgu133plus2hsentrezg
# notes=


eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDFplus2.txt")
RMAExpression_customCDFplus2<-read.delim("data_customCDFplus2.txt", sep="\t")
str(RMAExpression_customCDFplus2)
#'data.frame':	19764 obs. of  60 variables:
write.csv(RMAExpression_customCDFplus2, "RMAExpression_customCDFplus2.csv")


ScanDate<-protocolData(data2)$ScanDate
#Yep, there are definitely different scan dates here. And different date formats :(
#Oh wait - I think they are mostly just divided up by region
cbind(Tissue, ScanDate)
#Yep. Alright, never mind.

#Old code:
# library(reshape2)
# ScanDate_Split<-colsplit(ScanDate, pattern=" ", c("ScanDate", "ScanTime"))
# table(ScanDate_Split$ScanDate)


#####################################################
## This is an alternative version of the analysis that I ran to see if controlling for RNA degradation improved the results (especially since the PMI is so long for these analyses). Controlling for RIN improves our analyses of the Illumina data, but those results are more sensitive to the reliability of individual probes, whereas Affy uses probesets.
#This dataset actually has RIN values that accompany the data - I'll be curious to see how strongly they predict this.

source("https://bioconductor.org/biocLite.R")
biocLite("AffyRNADegradation")
library(AffyRNADegradation)

tongs <- GetTongs(data2, chip.idx = 4)
PlotTongs(tongs)
tongs <- GetTongs(data2, chip.idx = 5)
PlotTongs(tongs)

rna.deg<- RNADegradation(data2, location.type = "index")
RNADegradPerSample<-d(rna.deg)
str(RNADegradPerSample)

png("RNADegradPerSampleVsRIN.png")
plot(RNADegradPerSample~RIN)
dev.off()
#Definitely not perfect I wonder how much of that is due to the custom .cdf. I should come back later and try this without the custom .cdf.
# summary.lm(lm(RNADegradPerSample~RIN))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.136806   0.021866   6.256 2.29e-09 ***
#   RIN         0.037989   0.002841  13.372  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03419 on 203 degrees of freedom
# Multiple R-squared:  0.4683,	Adjusted R-squared:  0.4657 
# F-statistic: 178.8 on 1 and 203 DF,  p-value: < 2.2e-16

png("RNADegradPerSampleVsPMI.png")
plot(RNADegradPerSample~PMI)
dev.off()
#Interesting - the correlation with PMI is almost non-existent.
png("RINVsPMI.png")
plot(RIN~PMI)
dev.off()

#I'm going to test whether the custom cdf is screwing up the RNAdeg calculation:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/RawData/JustDLPFC")

data2<-ReadAffy()
# AffyBatch object
# size of arrays=1164x1164 features (44 kb)
# cdf=HG-U133_Plus_2 (54675 affyids)
# number of samples=68
# number of genes=54675
# annotation=hgu133plus2
# notes=
#eset2 <- rma(data2)
rna.deg<- RNADegradation(data2, location.type = "index")
RNADegradPerSample2<-d(rna.deg)
str(RNADegradPerSample2)
write.csv(RNADegradPerSample2, "RNADegradSample2.csv")

png("RNADegradPerSampleNoCustomCDFVsRIN.png")
plot(RNADegradPerSample2~RIN)
dev.off()

png("RNADegradPerSampleNoCustomCDFVsPMI.png")
plot(RNADegradPerSample2~PMI)
dev.off()

#Conclusion: Calculating RNADegradation using the original .cdf and the custom .cdf makes little difference.
rm(data2)

###################
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC")
Lanz_SampleCharacteristics<-data.frame(Lanz_SampleCharacteristics, RNADegradPerSample)
write.csv(Lanz_SampleCharacteristics, "Lanz_SampleCharacteristics_DLPFC.csv")

#Old code:
#' data3<-afbatch(rna.deg)
#' eset2 <- rma(data3)
#' write.exprs(eset2,file="data_customCDFplus2_RNADegCNTRL.txt")
#' RMAExpression_customCDFplus2_RNADegCNTRL<-read.delim("data_customCDFplus2_RNADegCNTRL.txt", sep="\t")
#' str(RMAExpression_customCDFplus2_RNADegCNTRL)
#' #'data.frame':	19764 obs. of  60 variables:
#' write.csv(RMAExpression_customCDFplus2_RNADegCNTRL, "RMAExpression_customCDFplus2_RNADegCNTRL.csv")
#' RMAExpression_customCDFplus2<-RMAExpression_customCDFplus2_RNADegCNTRL

#Old comments:
#Alright - after snooping at the PCA vs. RNA degradation plots, it looks like they were actually aggravated by this correction. I suspect it is because we are using a custom .cdf and the correction that they are applying is based on the index of the probe within the probeset - i.e. probably references the original probeset and not the custom .cdf mapping of probe to transcript.
#Here's some quotes from the users manual to back up my suspicion:
# "Instead of using the probe index within the probeset as argument of the degra- dation degree, one can use the actual probe locations within the transcript. We have pre-computed the distance of each probe to the 3’ end of its target transcript for all Affymetrix 3’ expression arrays. These probe location files are available under the URL http://www.izbi.uni-leipzig.de/downloads_ links/programs/rna_integrity.php."
#"It is possible to use custom probe locations, for example if one wishes to analyze custom built microarrays or if one relies on alternative probe annotations."
#My thought: this sort of problem is likely to bias the results for any particular probeset, but the average amount of degradation across probesets per individual sample is likely to be unaltered - I think I should probably just include this as a covariate in my analyses instead of using their corrected Affybatch values.

#####################################################


head(RMAExpression_customCDFplus2)
RMAExpression_EntrezIDplus2<-sub("_at", "", RMAExpression_customCDFplus2[,1])
head(RMAExpression_EntrezIDplus2)
RMAExpression_customCDFAnnotationplus2<-data.frame(RMAExpression_customCDFplus2[,1], RMAExpression_EntrezIDplus2, stringsAsFactors = F )
colnames(RMAExpression_customCDFAnnotationplus2)<-c("ProbesetID", "EntrezGeneID")

library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1 
# 59887 

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

RMAExpression_customCDFAnnotation2plus2<-join(RMAExpression_customCDFAnnotationplus2, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")


sum(is.na(RMAExpression_customCDFAnnotation2plus2[,3])==F)
#[1] 19528
dim(RMAExpression_customCDFAnnotation2plus2)
#[1] 19764     3
#So almost all of the results have gene symbols.

write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

SignalSortedNoNA3<-as.matrix(RMAExpression_customCDFplus2[,-1])

#Double-checking that the Sample info and Expression data are in the same order:
cbind(SampleID, colnames(SignalSortedNoNA3))
# Yep, they are both in ascending numeric order.

#Old code:
#Lanz_SampleCharacteristics<-data.frame(Lanz_SampleCharacteristics, ScanDate_Split)


####Quality Control################

setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/QC")

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",]

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]

str(RMAExpression_customCDFplus2)

png("Boxplot_RMAExpression_customCDFplus2.png", width=2000, height=400)
boxplot(SignalSortedNoNA3)
dev.off()
#Looks quantile normalized and the signal values now have an appropriate range. Hurray!

png("XIST_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]~Gender, col=2)
dev.off()
#Full datasets: Ooh - there are a couple of gender switches.
#DLPFC: Only one obvious mix-up
#let's double check that with another gene

png("RPS4Y1_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="RPS4Y1",][169,]~Gender, col=2)
dev.off()
#Full dataset: Yep, there are a whole slew of subjects with mistaken identity. :(
#DLPFC: Only one obvious mix-up

#Let's double check one more time:

png("DDX3Y_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="DDX3Y",][175,]~Gender, col=2)
dev.off()
#Full dataset: Yep, there are at least 7 subjects with mistaken identities.
#DLPFC: Still one mix-up.

png("RPS4Y1vsXIST_GenderCheck.png")
plot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="RPS4Y1",][169,]~SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,], col=as.numeric(Gender)+1)
dev.off()

png("DDX3YvsXIST_GenderCheck.png")
plot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="DDX3Y",][175,]~SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,], col=as.numeric(Gender)+1)
dev.off()


Lanz_SampleCharacteristics[Lanz_SampleCharacteristics$Gender=="M" & SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]>7,]

# SampleID Gender Age BrainpH  PMI                 Diagnosis      Tissue RIN
# GSM1304898_2_MDD_HPC.CEL.gz GSM1304898      M  44    6.50 19.3 major depressive disorder hippocampus 6.3
# GSM1304914_4_SCZ_HPC.CEL.gz GSM1304914      M  47    6.58 28.9             schizophrenia hippocampus 6.7
# RNADegradPerSample
# GSM1304898_2_MDD_HPC.CEL.gz          0.3903135
# GSM1304914_4_SCZ_HPC.CEL.gz          0.3179166

Lanz_SampleCharacteristics[Lanz_SampleCharacteristics$Gender=="F" & SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]<7,]


# SampleID Gender Age BrainpH  PMI
# GSM1304859_18_BP_HPC.CEL.gz                           GSM1304859      F  42     6.5 31.2
# GSM1304886_7_C_HPC.CEL.gz                             GSM1304886      F  58     6.6 18.8
# GSM1304889_11_MDD_HPC.CEL.gz                          GSM1304889      F  53     6.7 11.9
# GSM1304917_7_SCZ_HPC.CEL.gz                           GSM1304917      F  44     6.2 18.7
# GSM1304953_7_C_B46.CEL.gz                             GSM1304953      F  58     6.6 18.8
# GSM1305011_EA11002_238283_HG-U133_PLUS_2_STRC7.CEL.gz GSM1305011      F  58     6.6 18.8
# Diagnosis                    Tissue RIN
# GSM1304859_18_BP_HPC.CEL.gz                                    bipolar disorder               hippocampus 5.6
# GSM1304886_7_C_HPC.CEL.gz                                               control               hippocampus 7.2
# GSM1304889_11_MDD_HPC.CEL.gz                          major depressive disorder               hippocampus 8.1
# GSM1304917_7_SCZ_HPC.CEL.gz                                       schizophrenia               hippocampus 6.4
# GSM1304953_7_C_B46.CEL.gz                                               control Pre-frontal cortex (BA46) 8.7
# GSM1305011_EA11002_238283_HG-U133_PLUS_2_STRC7.CEL.gz                   control      Associative striatum 8.6
# RNADegradPerSample
# GSM1304859_18_BP_HPC.CEL.gz                                    0.3055808
# GSM1304886_7_C_HPC.CEL.gz                                      0.3853569
# GSM1304889_11_MDD_HPC.CEL.gz                                   0.3887742
# GSM1304917_7_SCZ_HPC.CEL.gz                                    0.3977488
# GSM1304953_7_C_B46.CEL.gz                                      0.4350562
# GSM1305011_EA11002_238283_HG-U133_PLUS_2_STRC7.CEL.gz          0.4239887

#Were any of these outliers actually removed in the original study? The hippocampal data is a mess - 6 samples out of 68 with *known* sample mix-ups (almost 10%!) - and there are probably more hiding in there that we can't easily see.
#Interesting. 

#Let's check PCA and see if there are any other obvious problem samples:

################################
# #Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/PCA")

setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/PCA")

pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")


PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, RMAExpression_customCDFAnnotation2plus2)
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

# #Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()


# #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Lanz_SampleCharacteristics$Tissue))
dev.off()
#Yes, there are three outliers identified through PCA too.

#After outlier removal:
png("10 PC1 vs PC2_byDiagnosis.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Diagnosis))
dev.off()


# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Lanz_SampleCharacteristics$Tissue))
dev.off()

#After outlier removal:
png("10 PC3 vs PC4_byGender.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Gender))
dev.off()

png("10 PC3 vs PC4_byDiagnosis.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers", col=as.factor(Diagnosis))
dev.off()

SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

PCAoutput<-cbind(SubjectFactorVariables, SubjectContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#Looks like there may be one outlier sample in the DLPFC data.

#######################################

#Visualize the sample-sample correlations using a heatmap:
png("09 Sample Sample Correlations Heatmap.png")
image(cor(SignalSortedNoNA3), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=1000, height=600)
boxplot(data.frame(cor(SignalSortedNoNA3)), cex=0.25, las=3, par(cex.axis=0.75, mar=c(15,4,4,4)), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()
#Using the full dataset: Beautiful data - no outliers.  Too bad they mixed up their sample identities. :(
#Using just the DLPFC: there is one obvious outlier. The data as a whole has a median sample-sample correlation of 0.993. GSM1304979 has a median sample-sample correlation of ~0.9775 and a range of sample-sample correlations that doesn't even overlap with the other samples.


#Outlier removal: 

#Full Dataset:
#KnownBadSamples<-c(SampleID[Lanz_SampleCharacteristics$Gender=="F" & SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]<7,], SampleID[Lanz_SampleCharacteristics$Gender=="M" & SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]>7,])

#DLPFC
KnownBadSamples<-c("GSM1304953", "GSM1304979")

SampleCharacteristics_NoOutliers<-Lanz_SampleCharacteristics[(SampleID%in%KnownBadSamples)==F,]
dim(SampleCharacteristics_NoOutliers)
#Full dataset: [1] 197   9
#DLPFC:[1] 66  9

SignalSortedNoNA3NoOutliers<-SignalSortedNoNA3[,(SampleID%in%KnownBadSamples)==F]
dim(SignalSortedNoNA3NoOutliers)
#[1] 19764    66

RNADegradPerSampleNoOutliers<-RNADegradPerSample[(SampleID%in%KnownBadSamples)==F]


#Redefining the variables to only include good data:
SignalSortedNoNA3<-SignalSortedNoNA3NoOutliers

Gender<-SampleCharacteristics_NoOutliers$Gender
Age<-SampleCharacteristics_NoOutliers$Age
PMI<-SampleCharacteristics_NoOutliers$PMI
BrainpH<-SampleCharacteristics_NoOutliers$BrainpH
Diagnosis<-SampleCharacteristics_NoOutliers$Diagnosis
Tissue<-SampleCharacteristics_NoOutliers$Tissue
RIN<-SampleCharacteristics_NoOutliers$RIN
RNADegradPerSample<-RNADegradPerSampleNoOutliers
#Old code:
#ScanDate<-SampleCharacteristics_NoOutliers$ScanDate
#ScanDate<-as.factor(ScanDate)

#I then re-ran QC and PCA with no outliers:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/PCA_NoOutliers")


png("RNADegradPerSampleVsRIN_noOutliers.png")
plot(RNADegradPerSample~RIN)
dev.off()

#This is from the full dataset
summary.lm(lm(RNADegradPerSample~RIN))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.140174   0.022536    6.22 2.97e-09 ***
#   RIN         0.037707   0.002922   12.91  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0338 on 195 degrees of freedom
# Multiple R-squared:  0.4607,	Adjusted R-squared:  0.4579 
# F-statistic: 166.6 on 1 and 195 DF,  p-value: < 2.2e-16

#That doesn't really improve things - the sample mix up must have occurred before they measured RIN.

#############################################

#Characterizing subjects and looking for relationships between subject variables:
#Changed Working directory
#Full dataset:
#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/SubjectVarVsSubjectVar")
#DLPFC:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/SubjVarVsSubjVar")

#Full Dataset:
#SubjectFactorVariables<-cbind(Diagnosis, Gender, Tissue)
#colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender", "Tissue")

#DLPFC
SubjectFactorVariables<-cbind(Diagnosis, Gender)
colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender")

SubjectContinuousVariables<-cbind(BrainpH, PMI, Age, RNADegradPerSample, RIN)
colnames(SubjectContinuousVariables)<-c("BrainpH", "PMI", "Age", "RNADegradPerSample", "RIN")

for (i in 1:length(SubjectContinuousVariables[1,])){
  png(paste(paste("Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep="."))	
  hist(SubjectContinuousVariables[, i], col=i+1)
  dev.off()		
}
#Wide age distribution, but mostly middle-aged
#wide pH distribution (5.8-7.4), but typically around 6.4-6.5
#PMI ranges from 5-35 hrs, typically ~20


#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}


#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}


#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("14 Cross Tabs Between Subject Factors.txt")
out<-c(
  capture.output(
    
    summary(Diagnosis),
    summary(Gender),
    
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(ContingencyTable)
        print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
      }		
    }
  )
)
cat(out, file="14 Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)


library(car)

StatisticalRelationshipsIV<-file("14 Statistical Relationships between Subject Variables.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. 
    vif(lm(SignalSortedNoNA3[1,]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample))
    
  ),
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(chisq.test(ContingencyTable))		
      }		
    }
  )
  
)
cat(out, file="14 Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)


#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  ),
  
#Using chi-square to examine the statistical relationships between the categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
        if(chisq.test(ContingencyTable)$p.value<0.05){
          print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))
        }else{}
      }
    }
  )
)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)

############################################

SubjectPCA<-SubjectPCA[(SampleID%in%KnownBadSamples)==F,]

#then looked at relationships with subject variables.

#Full dataset: setwd("~/Documents/Microarray Gen/Lanz_GSE53987/SubjectVarVsPCA")
#DLPFC:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/PCAvsSubjVar")

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectPCA[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
    RegressionLine<-lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}



#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectPCA[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectPCA)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}




#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("15 Statistical Relationships between Subject Variables and PCA.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
    summary.lm(lm(SubjectPCA[,1]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample))
  ),
  
  capture.output(
    summary.lm(lm(SubjectPCA[,2]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample))
  ),
  
  capture.output(
    summary.lm(lm(SubjectPCA[,3]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample))
  ),
  
  capture.output(
    summary.lm(lm(SubjectPCA[,4]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSample))
  ),
  
  
  #Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between PCA and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  )
  
)
cat(out, file="15 Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)



#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  )
)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)

#Full Dataset: Clearly tissue matters the most (duh)
#RNAdegradation seems better than RIN.

#####################

#Running Cell Type analyses:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/BrainInABlender")

install.packages("devtools")

library("devtools")

install_github("hagenaue/BrainInABlender")


library("BrainInABlender")

#DLPFC:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/BrainInABlender")

temp<-data.frame(as.character(RMAExpression_customCDFAnnotation2plus2[,3]), SignalSortedNoNA3, stringsAsFactors=F)
write.csv(temp, "UserInput.csv")
str(temp)
table(table(RMAExpression_customCDFAnnotation2plus2[,3]))
SirUnMixALotOutput<-Sir_UnMixALot(userInput=temp, dataColumns=c(2:67), geneColumn=1, species="human")

PublicationSpecific_CellTypeIndex<-SirUnMixALotOutput$PublicationSpecific_CellTypeIndex
AveragePrimary_CellTypeIndex<-SirUnMixALotOutput$AveragePrimary_CellTypeIndex

str(PublicationSpecific_CellTypeIndex)

#Wow- the output is really beautiful. Beautiful clusters.


temp<-cbind(as.matrix(SubjectContinuousVariables), t(PublicationSpecific_CellTypeIndex), t(AveragePrimary_CellTypeIndex)) 
str(temp)
SubjectContinuousVariables<-temp

########################################################################
#Then I re-ran the subjvar and PCA code.
#For the full dataset:
#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/SubjVarVsCellType")
#For the dlpfc:

row.names(AveragePrimary_CellTypeIndex)
for(i in c(1:10)){
  print(" 
        ")
  print(row.names(AveragePrimary_CellTypeIndex)[i])
  print(summary.lm(lm(AveragePrimary_CellTypeIndex[i,]~BrainpH+PMI+Diagnosis+Gender+Age+RNADegradPerSampleNoOutliers)))
}

#A version that doesn't include RNA degradation (to make it parallel our Affy analyses)
for(i in c(1:10)){
  print(" 
        ")
  print(row.names(AveragePrimary_CellTypeIndex)[i])
  print(summary.lm(lm(AveragePrimary_CellTypeIndex[i,]~BrainpH+PMI+Diagnosis+Gender+Age)))
}

#I need to re-run the PCA comparisons too.


write.csv(cor(cbind(as.numeric(SubjectFactorVariables[,2]), SubjectContinuousVariables), SubjectPCA), "Cor_SubjVar_vs_PCA.csv")




for(i in c(9:12)){

  print(paste((colnames(PCAoutput)[i]), "vsConfounds"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+BrainpH+PMI+Gender+Age+RNADegradPerSampleNoOutliers)))
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligos"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All)))
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoRNAdeg"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+RNADegradPerSampleNoOutliers))) 
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoConfounds"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+BrainpH+PMI+Gender+Age+RNADegradPerSampleNoOutliers))) 
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoAgePH"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+RNADegradPerSampleNoOutliers+Age+BrainpH))) 
 
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoMicroglia"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All))) 
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoMicroglia"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All+RNADegradPerSampleNoOutliers)))
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoMicrogliaConfounds"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All+BrainpH+PMI+Gender+Age+RNADegradPerSampleNoOutliers))) 
  
  print(paste((colnames(PCAoutput)[i]), "vsAstroOligoMicrogliaEndothelial"))
  print(summary.lm(lm(PCAoutput[,i]~Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All+Endothelial_All+RNADegradPerSampleNoOutliers))) 
 
  
}

  # colnames(PCAoutput)
  # [1] "Diagnosis"          "Gender"             "Tissue"             "BrainpH"            "PMI"               
  # [6] "Age"                "RNADegradPerSample" "RIN"                "PC1noOutliers"      "PC2noOutliers"     
  # [11] "PC3noOutliers"      "PC4noOutliers" 
  


###############################

#Full dataset:
#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/LM4_Basic")

#Just DLPFC:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/LM4_Basic")

#I need to output results from a basic linear model for comparison:
#Note - the p-values for this are going to be inflated because it should really be a hierarchical model with Tissue as a within subjects variable.

BrainpHcentered<-BrainpH-median(BrainpH)
Agecentered<-Age-median(Age)
PMIcentered<-PMI-median(PMI)
Tissue<-as.factor(Tissue)
Tissue<-relevel(Tissue, ref="Pre-frontal cortex (BA46)")
Lanz_SampleCharacteristics$Gender #Correct levels
levels(Gender)
Lanz_SampleCharacteristics$Diagnosis #Correct levels
levels(Diagnosis)

GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 9)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 9)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "Schiz", "RNADeg")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "Schiz", "RNADeg")

GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 9)
colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "Schiz", "RNADeg")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Tstat)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}

GeneByCellTypeSubjVar2_Pvalues2<-cbind(GeneByCellTypeSubjVar2_Pvalues, RMAExpression_customCDFAnnotation2plus2)
GeneByCellTypeSubjVar2_Betas2<-cbind(GeneByCellTypeSubjVar2_Betas, RMAExpression_customCDFAnnotation2plus2)

write.csv(GeneByCellTypeSubjVar2_Pvalues2, "GeneByCellTypeSubjVar2_Pvalues.csv")
write.csv(GeneByCellTypeSubjVar2_Betas2, "GeneByCellTypeSubjVar2_Betas.csv")

GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)

write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")



for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  png(paste(paste("17 Histogram of Raw Pvalues for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), "png", sep="."))	
  hist(GeneByCellTypeSubjVar2_Pvalues[,i], breaks=100, col=i, main=paste("Raw P-values for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
  abline(a=(length(GeneByCellTypeSubjVar2_Pvalues[,1])/100), b=0)
  dev.off()		
}	


GeneByCellTypeSubjVar2_PvaluesAdj<-matrix(0, length(GeneByCellTypeSubjVar2_Pvalues[,1]), length(GeneByCellTypeSubjVar2_Pvalues[1,]))
colnames(GeneByCellTypeSubjVar2_PvaluesAdj)<-colnames(GeneByCellTypeSubjVar2_Pvalues)
row.names(GeneByCellTypeSubjVar2_PvaluesAdj)<-row.names(GeneByCellTypeSubjVar2_Pvalues)


library(multtest)

for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  
  #Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
  TempPvalAdj<-mt.rawp2adjp(GeneByCellTypeSubjVar2_Pvalues[,i], proc=c("BH"))
  GeneByCellTypeSubjVar2_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]
  
}

GeneByCellTypeSubjVar2_PvaluesAdj2<-cbind(GeneByCellTypeSubjVar2_PvaluesAdj, RMAExpression_customCDFAnnotation2plus2)
write.csv(GeneByCellTypeSubjVar2_PvaluesAdj2, "GeneByCellTypeSubjVar2_PvaluesAdj.csv")

GeneByCellTypeSubjVar2DF<-as.data.frame(cbind(GeneByCellTypeSubjVar2_Betas, GeneByCellTypeSubjVar2_Tstat, GeneByCellTypeSubjVar2_Pvalues, GeneByCellTypeSubjVar2_PvaluesAdj))

temp<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2DF)
write.csv(temp, "GeneByCellTypeSubjVar2DF.csv" )


setwd("~/Documents/Microarray Gen/Lanz_GSE53987/LM4_Prevelent")

row.names(AveragePrimary_CellTypeIndex)
Astrocyte_All<-AveragePrimary_CellTypeIndex[1,]
Oligodendrocyte_All<-AveragePrimary_CellTypeIndex[8,]
Microglia_All<-AveragePrimary_CellTypeIndex[3,]
Endothelial_All<-AveragePrimary_CellTypeIndex[2,]
RBC_All<-AveragePrimary_CellTypeIndex[10,]
Neuron_All<-AveragePrimary_CellTypeIndex[5,]
Neuron_Interneuron<-AveragePrimary_CellTypeIndex[6,]
Neuron_Projection<-AveragePrimary_CellTypeIndex[7,]
Oligodendrocyte_Immature<-AveragePrimary_CellTypeIndex[9,]
Mural_All<-AveragePrimary_CellTypeIndex[4,]

summary.lm(lm(PC1noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+RNADegradPerSampleNoOutliers))

summary.lm(lm(PC2noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+RNADegradPerSampleNoOutliers))

summary.lm(lm(PC3noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+RNADegradPerSampleNoOutliers))

summary.lm(lm(PC4noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+RNADegradPerSampleNoOutliers))

#Full dataset - never outputted:
#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/LM4_Prevelent")

#DLPFC:
setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/LM4_PrevelentCellTypes")

GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 14)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 14)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "RNADeg")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "RNADeg")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)

GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 14)
colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "RNADeg")


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}



GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)

write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")




#setwd("~/Documents/Microarray Gen/Lanz_GSE53987/LM4_Everything")
#I don't know if the sample size for this dataset can actually handle this.

#DLPFC: 

GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 18)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 18)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "RNADeg")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "RNADeg")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)

GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 18)
colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "BP", "MDD", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "RNADeg")

for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Endothelial_All+Microglia_All+Mural_All+Neuron_All+Neuron_Projection+Neuron_Interneuron+Oligodendrocyte_All+RBC_All+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}




GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)

write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")



#I'm pretty sure that the signature captured by several fo those cell types is not the cell type of interest.


########################

#Running a version of BrainInABlender that doesn't remove Doyle:

Sir_UnMixALot_wDoyle(userInput=temp, dataColumns=c(2:67), geneColumn=1, species="human")

#######################

#Making pretty PC figures for the paper:

load("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC/Workspace_20170323.RData")

setwd("~/Documents/Microarray Gen/Lanz_GSE53987/DLPFC")

plot(SubjectPCA[,1]~Astrocyte_All)

pdf("Lanz_PC1vsAstrocyte_forPaper.pdf", width=4, height=4, pointsize=10)
plot(SubjectPCA[,1]~Astrocyte_All, xlab="Astrocyte Index", ylab="PC1", font.lab=2, lwd=1, cex.lab=1.3, cex.axis=1)
abline(lm(SubjectPCA[,1]~Astrocyte_All), lwd=4, col="dodgerblue3")
mtext(paste("r-squared: ", signif(summary.lm(lm(SubjectPCA[,1]~Astrocyte_All))$r.squared, digits=3), sep=""), side=3, line=2,cex=1.5, font=2)
mtext(paste("p-value: ", signif(summary.lm(lm(SubjectPCA[,1]~Astrocyte_All))$coefficients[2,4], digits=3), sep=""), side=3, line=0.5,cex=1.5, font=2)
dev.off()

