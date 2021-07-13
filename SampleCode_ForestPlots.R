#This is some example code for making forest plots to illustrate the results for the top genes in the meta-analysis:

library(metafor)

#In this example code, 
# "meta" is the data.frame that contains the effect sizes and variances from each study.
# the effect sizes are all include "d" in their column names, and the sampling variances all include "var d" in their column names.
# "Etv4" and "Sox9" are the genes of interest.
# the code following "slab" lists the study names

#Other things that will need changing:
#Unfortunately, we have a different number of studies included in our meta-analyses, so...
# the file size (width, height) in the tiff function is probably too big
#... and the x and y coordinates for the placement of text() will also be incorrect

#That means that overall, there will probably need to be *quite a bit of tweaking* to get this to look correct.

#Replace 
temp<-subset(meta, meta$Symbol=="Etv4" | meta$Symbol=="Sox9")

for(i in 1:length(temp$Symbol)){
  
  tiff(paste0(outGOIPlots, "Forest Plot ", 
              paste((temp$Symbol)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp$`d Affy F4`[i], temp$`d RNAseq F29`[i], temp$`d NC F34`[i], 
            temp$`d F37`[i], temp$`d F43`[i])
  var<-c(temp$`var d Affy F4`[i], temp$`var d RNAseq F29`[i], temp$`var d NC F34`[i], 
         temp$`var d F37`[i], temp$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29", 
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37",
                                         "MBNI_RNASeq_F43"), 
             xlim=c(-18, 15), cex=1.6)
  text(-11, 3.6, "Investigator & Study", cex=1.6)
  text(-14, 4.6, "bLR", cex=1.8)
  mtext(paste("Adult\n", temp$Symbol[i]), line=-1.5, cex=2.5)
  text(8, 3.6, "Cohen's D [95% CI]", cex=1.6)
  text(11, 4.6, "bHR", cex=1.8)
  
  dev.off()
}
