StretchEnhancers_Background_PvalueCalculation <- function(TestPeakSet_GR)
{
  setwd("/Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/StretchEnhancers/Z.Unique_StretchEnhancers/")
  Tissues <- read.table("/Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/TissuesStretchEnhancers.txt")
  Background_SEs <- read.table("0.All_StretchEnhancers_Merged.bed", header=FALSE)
  Background_SEs_GR <- GRanges(seqnames=Background_SEs$V1, ranges=IRanges(start=Background_SEs$V2, end=Background_SEs$V3))
  
  Tissue_SEs_Overlap_NonOverlap <- data.frame()
  Background_SEs_Overlap_NonOverlap <- data.frame()
  Background_SEs_Overlap_NonOverlap_ROWNAMES <- vector()
  
  for(i in 1:length(Tissues$V1))
  {
    Background_SEs_temp <- Background_SEs
    readFile <- paste("/Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/StretchEnhancers/Z.Unique_StretchEnhancers/", Tissues$V1[i], sep="")
    Tissue_StretchEnhancers <- read.table(readFile, header=FALSE)
    Tissue_StretchEnhancers_GR <- GRanges(seqnames=Tissue_StretchEnhancers$V1, ranges=IRanges(start=Tissue_StretchEnhancers$V2, end=Tissue_StretchEnhancers$V3))
    
    Overlaps <- findOverlaps(Tissue_StretchEnhancers_GR, TestPeakSet_GR, select="first")
    Tissue_StretchEnhancers[,ncol(Tissue_StretchEnhancers)+1] <- (!is.na(Overlaps))*1
    Tissue_SEs_Overlap_NonOverlap[i,1] <- sum(Tissue_StretchEnhancers[,ncol(Tissue_StretchEnhancers)])
    Tissue_SEs_Overlap_NonOverlap[i,2] <- nrow(Tissue_StretchEnhancers) - sum(Tissue_StretchEnhancers[,ncol(Tissue_StretchEnhancers)])
    
    Overlaps <- findOverlaps(Background_SEs_GR, Tissue_StretchEnhancers_GR, select="first")
    Background_SEs_temp[,ncol(Background_SEs_temp)+1] <- (!is.na(Overlaps))*1
    Background_SEs_temp <- Background_SEs_temp[which(Background_SEs_temp[,ncol(Background_SEs_temp)]==0),]
    Background_SEs_temp_GR <- GRanges(seqnames=Background_SEs_temp$V1, ranges=IRanges(start=Background_SEs_temp$V2, end=Background_SEs_temp$V3))
    
    Overlaps <- findOverlaps(Background_SEs_temp_GR, TestPeakSet_GR, select="first")
    Background_SEs_temp[,ncol(Background_SEs_temp)+1] <- (!is.na(Overlaps))*1
    Background_SEs_Overlap_NonOverlap[i,1] <- sum(Background_SEs_temp[,ncol(Background_SEs_temp)])
    Background_SEs_Overlap_NonOverlap[i,2] <- nrow(Background_SEs_temp) - sum(Background_SEs_temp[,ncol(Background_SEs_temp)])
    Background_SEs_Overlap_NonOverlap_ROWNAMES[i] <-  gsub("\\..*","",Tissues$V1[i])
  }
  rownames(Background_SEs_Overlap_NonOverlap) <- Background_SEs_Overlap_NonOverlap_ROWNAMES
  colnames(Background_SEs_Overlap_NonOverlap) <- c("Overlap", "NonOverlap")
  rownames(Tissue_SEs_Overlap_NonOverlap) <- Background_SEs_Overlap_NonOverlap_ROWNAMES
  colnames(Tissue_SEs_Overlap_NonOverlap) <- c("Overlap", "NonOverlap")
  
  # Calculate P-Values using Fisher's Exact Test
  pvalues <- nrow(Tissue_SEs_Overlap_NonOverlap)
  for(i in 1:nrow(Tissue_SEs_Overlap_NonOverlap))
  {
    ContingencyTable <- rbind(Tissue_SEs_Overlap_NonOverlap[i,], Background_SEs_Overlap_NonOverlap[i,])
    FisherTest <- fisher.test(ContingencyTable, alternative="greater")
    pvalues[i] <- FisherTest$p.value
  }
  names(pvalues) <- rownames(Tissue_SEs_Overlap_NonOverlap)
  return(sort(pvalues, decreasing=FALSE))
}
