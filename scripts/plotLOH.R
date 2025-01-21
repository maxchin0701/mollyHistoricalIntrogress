#### LOAD PACKAGES ####
library(karyoploteR)
library(viridisLite)
library(grDevices)

#### LOAD IN DATA ####
chrs <- read.delim(paste0("../data/chrIndex.tsv"),
           sep="\t",row.names = NULL,header = T)
#cutoff GQ
parentalSpecies <- "GLMV"
cutoffGQ <- 20
cutoffDP <- 6

#set up dfs to store
SNPall <- as.data.frame(matrix(NA,nrow=0,ncol=3))
LOHRegions <- as.data.frame(matrix(NA,nrow=0,ncol=4))
LOHAnc <- as.data.frame(matrix(NA,nrow=0,ncol=4))

for(i in unique(chrs$CHR)){
  # if(inherits(try(read.delim(paste0("../output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,".tsv"),sep="\t")), "try-error") ||
  #    inherits(try(read.delim(paste0("../output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,".tsv"),sep="\t")), "try-error") ||
  #    )
  curLOHRegions <- read.delim(paste0("../output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
  curLOHAnc <- read.delim(paste0("../output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
  curSNPall <- read.delim(paste0("../output/SNPall/",i,"SNPAll_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
  LOHRegions <- rbind(LOHRegions,curLOHRegions)
  LOHAnc <- rbind(LOHAnc,curLOHAnc)
  SNPall <- rbind(SNPall,curSNPall)
}

#### SETUP FOR PLOTTING ####

#create custom karyotype
pforKaryotype <- toGRanges(data.frame(chr=c(chrs$CHR), start=c(chrs$START), end=c(chrs$END)))
# pforKaryotype <- toGRanges(data.frame(chr=c("chr7_24"), start=c(1), end=c(chrs$END[7])),plot.type=2)
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=2)
kpAddCytobandsAsLine(pforKaryoteBase)

#generate color scales
nSampCols <- viridis(19,alpha=1,begin=0,end=1)
ancCols <- rocket(3,alpha=0.75,begin=0.25,end=1)

#loop through and add color color to regions df
for(i in 1:nrow(LOHRegions)){
  #LOHRegions$col[i] <- nSampCols[which(sort(unique(LOHRegions$group)) == LOHRegions$group[i])]
  LOHRegions$col[i] <- nSampCols[which(1:19 == LOHRegions$group[i])]
}

#loop through and add color to anc df
for(i in 1:nrow(LOHAnc)){
  LOHAnc$col[i] <- ancCols[which(c("plat","mix","pmex") == LOHAnc$group[i])]
}

#convert to grange object
SNPallGRegion <- toGRanges(SNPall)
LOHGRegions <- toGRanges(LOHRegions)
LOHGAnc <- toGRanges(LOHAnc)

#### PLOT ALL LOH REGIONS ####

#open cairo device
pdf(file=paste0("../figures/LOHPlotGQ",cutoffGQ,"DP",cutoffDP,"_",parentalSpecies,".pdf"),
          width=14.17,height=18.10)
#plot base karyotype
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=2)
kpAddCytobandsAsLine(pforKaryoteBase)

#plot regions
kpPlotRegions(pforKaryoteBase,data=LOHGRegions,col=LOHGRegions$col,border=NULL,avoid.overlapping = F,r0=0,r1=1)
#kpPlotRegions(pforKaryoteBase,data=SNPall,col="black",border=NULL,avoid.overlapping=F,r0=0.5,r1=1)
kpPlotRegions(pforKaryoteBase,data=LOHGAnc,col=LOHGAnc$col,data.panel=2,border=NULL,avoid.overlapping = F)
dev.off()


#### PLOT STRICTLY SHARED LOH REGIONS #####
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

#convert to grange object
LOHGRegionsShared <- toGRanges(LOHRegionsShared)
LOHGAncShared <- toGRanges(LOHAncShared)

#### PLOT SHARED LOH REGIONS ####

#open cairo device
pdf(file=paste0("../figures/LOHSharedPlotGQ",cutoffGQ,"DP",cutoffDP,"_",parentalSpecies,".pdf"),
    width=14.17,height=18.10)
#plot base karyotype
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=2)
kpAddCytobandsAsLine(pforKaryoteBase)

kpPlotRegions(pforKaryoteBase,data=LOHGRegionsShared,col=LOHGRegionsShared$col,border=NULL,avoid.overlapping = F,r0=0,r1=1)
#kpPlotRegions(pforKaryoteBase,data=SNPall,col="black",border=NULL,avoid.overlapping=T,r0=0.5,r1=1)
kpPlotRegions(pforKaryoteBase,data=LOHGAncShared,col=LOHGAncShared$col,data.panel=2,border=NULL,avoid.overlapping = F)
dev.off()
