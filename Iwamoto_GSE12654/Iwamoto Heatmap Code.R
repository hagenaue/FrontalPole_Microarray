library(car)

car::Anova(lm(PC1noOutliers~RNADegradPerSample+ScanDateDayOnly))

car::Anova(lm(PC1noOutliers~ScanDateDayOnly+RNADegradPerSample))

car::Anova(lm(PC2noOutliers~RateofDeath+ScanDateDayOnly))

plot(PC1noOutliers~PC2noOutliers, col=as.factor(ScanDateDayOnly), pch=18, cex=1.5)

Temp<-SignalSortedNoNA3
colnames(Temp)<-ScanDateDayOnly
heatmap(cor(Temp))


plot(PC2noOutliers~PC3noOutliers, col=as.factor(ScanDateDayOnly), pch=18, cex=1.5)


TempResiduals<-SignalSortedNoNA3

i<-1

for(i in c(1:nrow(SignalSortedNoNA3))) {
  TempModel<-lm(SignalSortedNoNA3[i,]~RNADegradPerSample+RateofDeath+BrainpH+PMI)
  TempResiduals[i,]<-TempModel$residuals
  rm(TempModel)
}

Temp<-TempResiduals

colnames(Temp)<-ScanDateDayOnly
heatmap(cor(Temp))


colnames(Temp)<-Diagnosis
heatmap(cor(Temp))

Temp<-SignalSortedNoNA3
colnames(Temp)<-Diagnosis
heatmap(cor(Temp))