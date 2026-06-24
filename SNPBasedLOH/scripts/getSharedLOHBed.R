# This script generates a BED file from shared LOH tracts.
# Outputs: 
#     LOHRegionsSharedCombined.bed: three column bed file
#     LOHAncSharedCombined.bed: five column bed file, column 4: ancestry 
#                               retained across tract, column 5:            
#                               length of tract

#### READ IN PARAMS ####
args <- commandArgs(trailingOnly = TRUE)
outRegions <- args[1]
outAnc <- args[2]

#### LOAD IN DATA ####
chrs <- read.delim(paste0("data/chrIndex.tsv"),
                   sep="\t",row.names = NULL,header = T)
# colnames(genomeCN)[5:23] <- gsub(".","-",colnames(genomeCN)[5:23],
#                                  fixed=T)

#cutoff GQ
parentalSpecies <- "LM"
cutoffGQ <- 10
cutoffDP <- 6

#set up dfs to store
SNPall <- as.data.frame(matrix(NA,nrow=0,ncol=3))
LOHRegions <- as.data.frame(matrix(NA,nrow=0,ncol=4))
LOHAnc <- as.data.frame(matrix(NA,nrow=0,ncol=4))

for(i in unique(chrs$CHR)){
  if(inherits(try(read.delim(paste0("output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")), "try-error") ||
     inherits(try(read.delim(paste0("output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")), "try-error") ||
     inherits(try(read.delim(paste0("output/SNPall/",i,"SNPAll_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")), "try-error"))
  {
    next()
  } else {
    curLOHRegions <- read.delim(paste0("output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
    curLOHAnc <- read.delim(paste0("output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
    curSNPall <- read.delim(paste0("output/SNPall/",i,"SNPAll_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
    LOHRegions <- rbind(LOHRegions,curLOHRegions)
    LOHAnc <- rbind(LOHAnc,curLOHAnc)
    SNPall <- rbind(SNPall,curSNPall)
  }
}

#subset to just 19 sample individuals
LOHRegionsShared <- LOHRegions[which(LOHRegions[,4] >= 18),]
LOHRegionsShared$intSize <- LOHRegionsShared$end - LOHRegionsShared$start

#get intervals
sharedAncIndex <- c()
for(i in unique(LOHAnc$chr)){
  if(nrow(LOHRegionsShared[which(LOHRegionsShared[,1] == i),]) == 0){
    next
  } else {
    curChrLOHRegionsShared <- LOHRegionsShared[which(LOHRegionsShared[,1] == 
                                                       i),]
    curChrLOHSharedIntervals <- c()
    for(j in 1:nrow(curChrLOHRegionsShared)){
      curChrLOHSharedIntervals <- c(curChrLOHSharedIntervals,
                                    curChrLOHRegionsShared[j,2]:curChrLOHRegionsShared[j,3])
    }
    
    sharedAncIndex <- c(sharedAncIndex,
                        which(LOHAnc[,1] == i & (LOHAnc[,2] %in% curChrLOHSharedIntervals &
                                                   LOHAnc[,3] %in% curChrLOHSharedIntervals)))
  }
}

#subset anc
LOHAncShared <- LOHAnc[sharedAncIndex,]

#get index list
LOHAncShared$intSize <- LOHAncShared$end - LOHAncShared$start

#save shared regions
LOHRegionsSharedLong <- LOHRegionsShared[which(LOHRegionsShared$intSize >= 1000),1:3]

write.table(LOHRegionsSharedLong,
            outRegions,
            quote=FALSE, sep='\t',row.names=F,col.names = F)

#save shared ancestry
LOHAncSharedLong <- LOHAncShared[which(LOHAncShared$intSize >= 1000),]

write.table(LOHRegionsSharedLong,
            outAnc,
            quote=FALSE, sep='\t',row.names=F,col.names = F)



