#This is some sample code for making volcano plots:

#In this code:
# metaAnalysisOutputFDR is the data.frame that includes the effect size, p-values, and FDR
# estimate is the effect size estimate from the meta-analysis. For an individual study (e.g., Iwamoto), this would be the log2 fold change (coefficient) in the limma output
# pval means nominal p-value
# padj or FDR means p-value adjusted for false detection rate (BH adjusted)

#This plot emphasizes results with an padj cut off of <0.05 - we may want to up that to 0.1

tiff(paste0(outDir, "VolcanoPlotAdult.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

par(mai=c(1.02, 1,0.9,0.40))
with(metaAnalysisOutputFDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                 xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))
# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(metaAnalysisOutputFDR, abs(estimate)>1), points(estimate, -log10(pval), 
                                                            pch=19, col="red", cex=0.6))
with(subset(metaAnalysisOutputFDR, BH< .05 ), points(estimate, -log10(pval), 
                                                     pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()
