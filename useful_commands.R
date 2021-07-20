if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
nBiocManager::install("oligoClasses")
BiocManager::install("oligo")

BiocManager::install("AffyRNADegradation")

BiocManager::install("fgsea")
BiocManager::install("multtest")

ainstall.packages("C:/Users/Frosty/Desktop/Research/Summer 2021 Frontal Pole Research/FrontalPole_Microarray/affy_1.68.0.tar.gz", repos = NULL, type = "source")


if(!require(installr)) {
  install.packages("installr"); 
  require(installr)
}




options(error=recover)



#removes everything in r environment:
rm(list = ls())


sum(duplicated(IwamotoSchizEffectSizes$EntrezGeneID))