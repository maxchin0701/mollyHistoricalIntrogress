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

#### SET UP COPY NUMBER DF ####
# keepCN <- c()
# for(i in 1:nrow(LOHRegions)){
#   keepCN <- c(keepCN,which(genomeCN$chr==LOHRegions$chr[i] & 
#           ((genomeCN$start<=LOHRegions$start[i] &
#              genomeCN$end>=LOHRegions$start[i]) | 
#           (genomeCN$start>=LOHRegions$start[i] &
#              genomeCN$end<=LOHRegions$end[i]) |
#           (genomeCN$start<=LOHRegions$end[i] &
#              genomeCN$end>=LOHRegions$end[i]))))
# }
# lohCN <- genomeCN[keepCN,]  
# 
# sampWideAvgL<- list()
# for(i in 1:nrow(lohCN)){
#   sampWideAvgL[[length(sampWideAvgL) + 1]] <- mean(as.numeric(lohCN[i,5:23]))
# }
# 
# sampWideAvgL<- list()
# for(i in 1:nrow(lohCN)){
#   sampWideAvgL[[length(sampWideAvgL) + 1]] <- mean(as.numeric(lohCN[i,5:23]))
# }
# 
# 
# lohCN$sampWideAvg <- do.call(rbind,sampWideAvgL)

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

#### PLOT ALL LOH REGIONS WITH CN ####
#open cairo device
pdf(file=paste0("../figures/LOHCNPlotGQ",cutoffGQ,"DP",cutoffDP,"_",parentalSpecies,".pdf"),
    width=14.17,height=18.10)
#plot base karyotype
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=2)
kpAddCytobandsAsLine(pforKaryoteBase)

#plot regions
kpPlotRegions(pforKaryoteBase,data=LOHGRegions,col=LOHGRegions$col,border=NULL,avoid.overlapping = F,r0=0,r1=1)
#kpPlotRegions(pforKaryoteBase,data=SNPall,col="black",border=NULL,avoid.overlapping=F,r0=0.5,r1=1)
kpLines(pforKaryoteBase,chr=lohCN$chr,x=lohCN$mean,y=lohCN$sampWideAvg,ymin=0,ymax=4,
        data.panel=2,r0=1,r1=0)
kpAxis(pforKaryoteBase,data.panel = 2,ymin=0,ymax=4,r0=1,r1=0,cex=0.75)
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

#save shared regions
LOHRegionsSharedLong <- LOHRegionsShared[which(LOHRegionsShared$intSize >= 1000),1:3]

write.table(LOHRegionsSharedLong,
      "../output/LOHRegionsSharedCombined/LOHRegionsSharedCombined.bed",
      quote=FALSE, sep='\t',row.names=F,col.names = F)

#### PLOT SHARED LOH REGIONS CN ####

#open cairo device
pdf(file=paste0("../figures/LOHCNSharedPlotGQ",cutoffGQ,"DP",cutoffDP,"_",parentalSpecies,".pdf"),
    width=14.17,height=18.10)
#plot base karyotype
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=2)
kpAddCytobandsAsLine(pforKaryoteBase)

kpPlotRegions(pforKaryoteBase,data=LOHGRegionsShared,col=LOHGRegionsShared$col,border=NULL,avoid.overlapping = F,r0=0,r1=1)
#kpPlotRegions(pforKaryoteBase,data=SNPall,col="black",border=NULL,avoid.overlapping=T,r0=0.5,r1=1)
# for(i in 1:nrow(LOHRegionsShared)){
#   curLOHCN <- lohCN[which(lohCN$chr==LOHRegionsShared$chr[i] &
#                 ((lohCN$start<=LOHRegionsShared$start[i] &
#                     lohCN$end>=LOHRegionsShared$start[i]) |
#                    (lohCN$start>=LOHRegionsShared$start[i] &
#                       lohCN$end<=LOHRegionsShared$end[i]) |
#                    (lohCN$start<=LOHRegionsShared$end[i] &
#                       lohCN$end>=LOHRegionsShared$end[i]))),]
#   kpLines(pforKaryoteBase,chr=curLOHCN$chr,x=curLOHCN$mean,y=curLOHCN$sampWideAvg,ymin=0,ymax=4,
#           data.panel=2,r0=1,r1=0)
# }
kpLines(pforKaryoteBase,chr=lohCN$chr,x=lohCN$mean,y=lohCN$sampWideAvg,ymin=0,ymax=4,
        data.panel=2,r0=1,r1=0)
kpAxis(pforKaryoteBase,data.panel = 2,ymin=0,ymax=4,r0=1,r1=0,cex=0.75)
dev.off()

