#Outputting boxplots for Adriana's top genes:


#Sample code uses "ABAT" as the example gene of interest - replace with whichever of the genes you are exploring.

#Replace "SignalSortedNoNA3" with whatever the name is of your log2 gene expression matrix with the outliers removed
#This code assumes that the gene symbols are the row.names for the log2 gene expression matrix - if not, the code will need to be changed to point correctly to where the gene annotation for each row is located.

#Replace "SubjectInfo" with whatever the name is of your subject metadata data.frame with the outliers removed
#"Exploratory Variable" in the code should be replaced with whichever variable you are looking at in Iwamoto (e.g., SuicideStatus)

#Depending on the number of diagnoses in the dataset, you will have a different number of boxes in the plot, so the col=c(...) will change to reflect the number of boxes, with the number-code indicating the diagnosis of the box:
#Colors
#Controls are currently colored red (2)
#BP is currently colored green (3)
#Schiz is currently colored blue (4)
#We don't have a color yet for MDD - maybe make teal? (5)

#If there are more diagnoses, you may need to change the width of the pdf to be larger
#If there are less diagnoses, you may need to change the width of the pdf to be smaller


#Simple diagnosis boxplot:

pdf("Boxplot_ABAT_ByDiagnosis.pdf", width=4, height=6)

boxplot(SignalSortedNoNA3[row.names(SignalSortedNoNA3)=="ABAT",]~SubjectInfo$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="ABAT", outline=FALSE)
stripchart(SignalSortedNoNA3[row.names(SignalSortedNoNA3)=="ABAT",]~SubjectInfo$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext("ABAT Log(2) Signal", cex=1.3, side=2, line=2.5)
dev.off()


#Diagnosis with Explanatory Variable boxplot:

pdf("Boxplot_ABAT_ByDiagnosis_byExploratoryVariable.pdf", width=10, height=6)

boxplot(SignalSortedNoNA3[row.names(SignalSortedNoNA3)=="ABAT",]~SubjectInfo$ExploratoryVariable+SubjectInfo$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4), main="ABAT", outline=FALSE)
stripchart(SignalSortedNoNA3[row.names(SignalSortedNoNA3)=="ABAT",]~SubjectInfo[,i]+SubjectInfo$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext("ABAT Log(2) Signal", cex=1.3, side=2, line=2.5)
dev.off()
