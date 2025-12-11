#### LOAD PACKAGES ####
library(karyoploteR)
library(viridisLite)
library(grDevices)

#### LOAD IN DATA ####
chrs <- read.delim(paste0("../../SNPBasedLOH/data/chrIndex.tsv"),
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
  if(inherits(try(read.delim(paste0("../../SNPBasedLOH/output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")), "try-error") ||
     inherits(try(read.delim(paste0("../../SNPBasedLOH/output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")), "try-error") ||
     inherits(try(read.delim(paste0("../../SNPBasedLOH/output/SNPall/",i,"SNPAll_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")), "try-error"))
  {
    next()
  } else {
    curLOHRegions <- read.delim(paste0("../../SNPBasedLOH/output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
    curLOHAnc <- read.delim(paste0("../../SNPBasedLOH/output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
    curSNPall <- read.delim(paste0("../../SNPBasedLOH/output/SNPall/",i,"SNPAll_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
    LOHRegions <- rbind(LOHRegions,curLOHRegions)
    LOHAnc <- rbind(LOHAnc,curLOHAnc)
    SNPall <- rbind(SNPall,curSNPall)
  }
}

LOHRegions$intSize <- LOHRegions$end - LOHRegions$start

#### SETUP FOR PLOTTING ####
#create custom karyotype
pforKaryotype <- toGRanges(data.frame(chr=c(chrs$CHR), start=c(chrs$START), end=c(chrs$END)))
# pforKaryotype <- toGRanges(data.frame(chr=c("chr7_24"), start=c(1), end=c(chrs$END[7])),plot.type=2)
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=2)
kpAddCytobandsAsLine(pforKaryoteBase)

#generate color scales
ancCols <- c("#0000FF","#7A3990","#F37221")

#loop through and add color to anc df
for(i in 1:nrow(LOHAnc)){
  LOHAnc$col[i] <- ancCols[which(c("plat","mix","pmex") == LOHAnc$group[i])]
}

#get shared LOH
LOHRegionsShared <- LOHRegions[which(LOHRegions[,4] >= 18),]

#get shared intervals
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
LOHAncShared$intSize <- LOHAncShared$end - LOHAncShared$start

#convert to grange object
LOHGAncShared <- toGRanges(LOHAncShared)

#plot base karyotype
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=0.75,font=2)
kpAddCytobandsAsLine(pforKaryoteBase)
kpAddBaseNumbers(pforKaryoteBase,tick.dist=5000000,minor.tick.dist = 1000000,cex=0.5,font=2)

#plot regions
kpPlotRegions(pforKaryoteBase,data=LOHGAncShared,col=LOHGAncShared$col,
              data.panel=1,border=NULL,avoid.overlapping = F)


