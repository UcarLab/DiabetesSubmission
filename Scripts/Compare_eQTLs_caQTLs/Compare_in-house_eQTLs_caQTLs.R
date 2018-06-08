library(lattice); library(ggplot2)
caQTL_pvalue_cutoff_Value <- 0.02276441
setwd("/Users/skn/Dropbox/Diabetes_FinalSubmission_RemoveBatchEffectArtifacts/Z.Scripts/Compare_eQTLs_caQTLs/")
# Reading in-house eQTLs
Lead_caQTLs_eQTLs <- read.table("Islets19_WithCovariates_Rasqual_eQTL_Results_Output_WithWindows_50kb_AllSNPs_WithBatchEffects_UnwantedRowsRemoved_Filtered.txt", header=TRUE)
# Reading RASQUAL caQTLs
Lead_caQTLs_caQTLs <- read.table("WithoutWindows_All_caQTLs_temp.txt", header=TRUE)
# Merging the two files
caQTLs_eQTLs <- unique(merge(Lead_caQTLs_caQTLs, Lead_caQTLs_eQTLs, by="rsID"))
caQTLs_eQTLs$caQTL_or_Not <- as.factor((caQTLs_eQTLs$AdjustedPVal_1<caQTL_pvalue_cutoff_Value)*1)
rm(Lead_caQTLs_eQTLs)
# QQPlot
set.seed(4732874)
x <- -log10(sort(caQTLs_eQTLs$PValue.y[which(caQTLs_eQTLs$AdjustedPVal_1<caQTL_pvalue_cutoff_Value)], decreasing = FALSE))
y_prime <- sample(caQTLs_eQTLs$PValue.y[which(caQTLs_eQTLs$AdjustedPVal_1>caQTL_pvalue_cutoff_Value)], size=length(x)); y <- -log10(sort(y_prime, decreasing = FALSE))
e = -log10(1:length(x)/length(x))
plot(e,x,pch=16,cex=1, col="black",
     xlab=expression(Expected~~-log[10](italic(p))),
     ylab=expression(Observed~~-log[10](italic(p))),
     xlim=c(0,max(e)), ylim=c(0,max(c(x,y))))
lines(e,e,col="red")
points(e, y, col="blue", pch=17)