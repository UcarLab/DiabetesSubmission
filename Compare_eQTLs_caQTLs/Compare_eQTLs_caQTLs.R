library(lattice); library(ggplot2)
caQTL_pvalue_cutoff_Value <- 0.02276441
setwd("/Users/skn/Dropbox/Diabetes Submission/Z.Revision/Z.Scripts/Compare_eQTLs_caQTLs/")
# Reading Varshney eQTLs (All)
Lead_caQTLs_eQTLs <- read.table("islet_meta.results.annotated.qvalue.postfreq.bestEqtl.tbl", header=TRUE)
hist(as.numeric(Lead_caQTLs_eQTLs$dist_tss), breaks=1000)
Lead_caQTLs_eQTLs <- Lead_caQTLs_eQTLs[order(abs(Lead_caQTLs_eQTLs$dist_tss), decreasing = FALSE),]
nrow(Lead_caQTLs_eQTLs[which(abs(Lead_caQTLs_eQTLs$dist_tss)<10000),])

# Reading Varshney eQTLs
Lead_caQTLs_eQTLs <- read.table("islet_meta.results.annotated.qvalue.postfreq.mod.leadCaQTLsnpsOnly.tbl")
colnames(Lead_caQTLs_eQTLs) <- c("SNP",	"gene",	"Allele1",	"Allele2",	"Weight",	"Zscore",	"Direction",	"HetISq",	"HetChiSq",	"HetDf",	"logHetP",	"log_pval",	"pvalue",	"q_storey",	"chrom",	"pos",	"dist",	"dist_tss",	"gene_type",	"gene_name",	"gene_chrom")
# Reading RASQUAL caQTLs
Lead_caQTLs_caQTLs <- read.table("WithoutWindows_All_caQTLs_temp.txt", header=TRUE)
# Merging the two files
caQTLs_eQTLs <- unique(merge(Lead_caQTLs_caQTLs, Lead_caQTLs_eQTLs, by.x="rsID", by.y="SNP"))
caQTLs_eQTLs$caQTL_or_Not <- as.factor((caQTLs_eQTLs$AdjustedPVal_1<caQTL_pvalue_cutoff_Value)*1)

# Batch effect associated caQTLs
BatchEffectAssociated_FeatureIDs <- read.table("BatchEffectAssociated_FeatureIDs.txt", header=FALSE)
caQTLs_eQTLs$caQTL_or_Not[which(caQTLs_eQTLs$FeatureID %in% BatchEffectAssociated_FeatureIDs$V1)] <- 0

# Proportion of Concordant caQTL-eQTL pairs
proportionConcordant <- vector()
proportionConcordant_Count <- vector()
proportionConcordant_xaxis <- vector()
for(i in 1:10)
{
    caQTLs_eQTLs_SigOnly <- caQTLs_eQTLs[which(caQTLs_eQTLs$caQTL_or_Not==1 & caQTLs_eQTLs$pvalue<(1/10^i)),]
    proportionConcordant[i] <- (nrow(caQTLs_eQTLs_SigOnly[which(caQTLs_eQTLs_SigOnly$EffectSize<0.5 & caQTLs_eQTLs_SigOnly$Zscore<0),]) + nrow(caQTLs_eQTLs_SigOnly[which(caQTLs_eQTLs_SigOnly$EffectSize>0.5 & caQTLs_eQTLs_SigOnly$Zscore>0),]))/nrow(caQTLs_eQTLs_SigOnly)
    proportionConcordant_Count[i] <- nrow(caQTLs_eQTLs_SigOnly)
    proportionConcordant_xaxis[i] <- 1/10^i
}
plot(-log10(proportionConcordant_xaxis), proportionConcordant, type="o", xlab="-log10(eQTL pvalue)", ylab="Proportion of Concordant caQTL-eQTL pairs")

proportionConcordant <- vector()
proportionConcordant_Count <- vector()
proportionConcordant_xaxis <- vector()
for(i in 1:10)
{
  caQTLs_eQTLs_SigOnly <- caQTLs_eQTLs[which(caQTLs_eQTLs$PValue<(1/10^i) & caQTLs_eQTLs$q_storey<0.1),]
  proportionConcordant[i] <- (nrow(caQTLs_eQTLs_SigOnly[which(caQTLs_eQTLs_SigOnly$EffectSize<0.5 & caQTLs_eQTLs_SigOnly$Zscore<0),]) + nrow(caQTLs_eQTLs_SigOnly[which(caQTLs_eQTLs_SigOnly$EffectSize>0.5 & caQTLs_eQTLs_SigOnly$Zscore>0),]))/nrow(caQTLs_eQTLs_SigOnly)
  proportionConcordant_Count[i] <- nrow(caQTLs_eQTLs_SigOnly)
  proportionConcordant_xaxis[i] <- 1/10^i
}
plot(-log10(proportionConcordant_xaxis), proportionConcordant, type="o", xlab="-log10(caQTL pvalue)", ylab="Proportion of Concordant caQTL-eQTL pairs")

# Examining direction-of-effect
caQTLs_eQTLs_SigOnly <- caQTLs_eQTLs[which(caQTLs_eQTLs$caQTL_or_Not==1 & caQTLs_eQTLs$q_storey<0.1),]
plot(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore, xlab="Effect Size of caQTLs", ylab="Z-Score of eQTLs", xlim=c(0.2,0.8), ylim=c(-10,12),
     main=paste("Looking at significant caQTLs and significant eQTLs. Cor", cor(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore), sep=":"))
abline(h=0, col="red"); abline(v=0.5, col="red")

# Examining direction-of-effect
caQTLs_eQTLs_SigOnly <- caQTLs_eQTLs[which(caQTLs_eQTLs$caQTL_or_Not!=1 & caQTLs_eQTLs$q_storey<0.1),]
plot(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore, xlab="Effect Size of caQTLs", ylab="Z-Score of eQTLs", xlim=c(0.2,0.8), ylim=c(-10,12),
     main=paste("Looking at non-significant caQTLs and significant eQTLs. Cor", cor(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore), sep=":"))
abline(h=0, col="red"); abline(v=0.5, col="red")

# Examining direction-of-effect (Only Significant caQTLs)
caQTLs_eQTLs_SigOnly <- caQTLs_eQTLs[which(caQTLs_eQTLs$caQTL_or_Not==1),]
plot(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore, xlab="Effect Size of caQTLs", ylab="Z-Score of eQTLs", 
     main=paste("Looking at significant caQTLs only. Cor", cor(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore), sep=":"))
abline(h=0, col="red"); abline(v=0.5, col="red")

# Examining direction-of-effect (Non-Significant caQTLs)
caQTLs_eQTLs_SigOnly <- caQTLs_eQTLs[which(caQTLs_eQTLs$caQTL_or_Not!=1),]
plot(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore, xlab="Effect Size of caQTLs", ylab="Z-Score of eQTLs", 
     main=paste("Looking at non-significant caQTLs only. Cor", cor(caQTLs_eQTLs_SigOnly$EffectSize, caQTLs_eQTLs_SigOnly$Zscore), sep=":"))
abline(h=0, col="red"); abline(v=0.5, col="red")

# QQPlot
set.seed(4732874)
x <- -log10(sort(caQTLs_eQTLs$pvalue[which(caQTLs_eQTLs$AdjustedPVal_1<caQTL_pvalue_cutoff_Value)], decreasing = FALSE))
y_prime <- sample(caQTLs_eQTLs$pvalue[which(caQTLs_eQTLs$AdjustedPVal_1>caQTL_pvalue_cutoff_Value)], size=length(x)); y <- -log10(sort(y_prime, decreasing = FALSE))
e = -log10(1:length(x)/length(x))
plot(e,x,pch=16,cex=1, col="black",
     xlab=expression(Expected~~-log[10](italic(p))),
     ylab=expression(Observed~~-log[10](italic(p))),
     xlim=c(0,max(e)), ylim=c(0,max(c(x,y))))
lines(e,e,col="red")
points(e, y, col="blue", pch=17)

# Examining Luc or T2D-associated caQTLs 
Luc_T2D_Associated <- read.table("/Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/All_Islets_caQTL_Analysis/Rasqual_Results_Output/RASQUAL_Results/SNPs_Only_Within_Peaks/Without_Windows_WithBatchEffects/T2D-Associated-Significant-caQTLs.txt", header=FALSE)
plot(caQTLs_eQTLs[which(caQTLs_eQTLs$FeatureID %in% Luc_T2D_Associated$V1),"EffectSize"], caQTLs_eQTLs[which(caQTLs_eQTLs$FeatureID %in% Luc_T2D_Associated$V1),"Zscore"])
abline(h=0, col="red"); abline(v=0.5, col="red")

caQTLs_eQTLs_SigOnly[which(caQTLs_eQTLs_SigOnly$EffectSize<0.5 & caQTLs_eQTLs_SigOnly$Zscore>0),]
