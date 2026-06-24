# This script compares LOH tract sizes across shared (>18) and isolated (<6) 
# tracts.
# output:
#   LOHSizeSigDiff.tsv: tsv with two columns: statistical test and 
#                       proportion of corrected p-values significant (< 0.05)

#### LIBRARY ####
library(twosamples)

#### SET SEED ####
set.seed(1)

#### READ IN PARAMS ####
args <- commandArgs(trailingOnly = TRUE)
outSummary <- args[1]

#### LOAD IN DATA ####
chrs <- read.delim(paste0("data/chrIndex.tsv"),
                   sep="\t",row.names = NULL,header = T)
# colnames(genomeCN)[5:23] <- gsub(".","-",colnames(genomeCN)[5:23],
#                                  fixed=T)

#cutoff GQ
parentalSpecies <- "LM"
cutoffGQ <- 10
cutoffDP <- 6

LOHRegions <- as.data.frame(matrix(NA,nrow=0,ncol=4))

#build datafra
for(i in unique(chrs$CHR)){
  curLOHRegions <- read.delim(paste0("output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
  LOHRegions <- rbind(LOHRegions,curLOHRegions)
}

LOHRegions$intSize <- LOHRegions$end - LOHRegions$start

#### SEPARATE OUT OVERLAPPING AND NONOVERLAPPING ####
LOHRegionsIsolated <- LOHRegions[which(LOHRegions$group <= 6 & 
                                         LOHRegions$intSize > 1000),]
LOHRegionsShared <- LOHRegions[which(LOHRegions$group >= 18 &
                                       LOHRegions$intSize > 1000),]

#### SAMPLE AND COMPARE ####
sigDiff <- list()
sigDiffMWU <- list()
sigDiffCVM <- list()
for(i in 1:1000){
  #sample intSizes
  curIsolatedSample <- sample(LOHRegionsIsolated$intSize,size=65,replace=T)
  
  #two sided ks test
  pValKS <- ks.test(LOHRegionsShared$intSize,curIsolatedSample)[2]
  pvalMWU <- wilcox.test(curIsolatedSample,LOHRegionsShared$intSize)[3]
  pvalCVM <- cvm_test(curIsolatedSample,LOHRegionsShared$intSize,nboots=1000)[2]
  sigDiff[[length(sigDiff) + 1]] <- pValKS
  sigDiffMWU[[length(sigDiffMWU) + 1]] <- pvalMWU
  sigDiffCVM[[length(sigDiffCVM) + 1]] <- pvalCVM
}
sigDiff <- unlist(sigDiff)
sigDiffMWU <- unlist(sigDiffMWU)
sigDiffCVM <- unlist(sigDiffCVM)
propSigKS <- sum(p.adjust(sigDiff,method="BH")<=0.05)/1000
propSigMWU <- sum(p.adjust(sigDiffMWU,method="BH")<=0.05)/1000
propSigCVM <- sum(p.adjust(sigDiffCVM,method="BH")<=0.05)/1000

#### SAVE ####
propSigDf <- as.data.frame(matrix(c("Kolmogorov-Smirnow", propSigKS,
                       "Mann-Whitney U", propSigMWU,
                       "Cramer-Von Mises", propSigCVM),nrow=3,byrow=T)) 

colnames(propSigDf) <- c("test", "propPValsSig")
#save
write.table(propSigDf, 
            file=outSummary, 
            quote=FALSE, sep='\t',row.names=FALSE)










