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
  curLOHRegions <- read.delim(paste0("../../SNPBasedLOH/output/LOHRegions/",i,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
  LOHRegions <- rbind(LOHRegions,curLOHRegions)
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
colorGrad <- colorRampPalette(c("#0000FF", "#F37221"))
nSampCols <- colorGrad(19)

#loop through and add color color to regions df
for(i in 1:nrow(LOHRegions)){
  #LOHRegions$col[i] <- nSampCols[which(sort(unique(LOHRegions$group)) == LOHRegions$group[i])]
  LOHRegions$col[i] <- nSampCols[which(1:19 == LOHRegions$group[i])]
}

#subset to just large regions
LOHRegions <- LOHRegions[which(LOHRegions$intSize >= 1000),]

#convert to grange object
LOHGRegions <- toGRanges(LOHRegions)

#### PLOT ALL LOH REGIONS ####

#open cairo device
pdf(file=paste0("suppFigAllTracts.pdf"),
    width=6.49,height=5.15)
#plot base karyotype
pforKaryoteBase <- plotKaryotype(genome = pforKaryotype,chromosomes="all",
                                 ideogram.plotter = NULL,
                                 plot.type=2,cex=0.75,font=2)
kpAddCytobandsAsLine(pforKaryoteBase)
kpAddBaseNumbers(pforKaryoteBase,tick.dist=5000000,minor.tick.dist = 1000000,cex=0.5,font=2)
kpPlotRegions(pforKaryoteBase,data=LOHGRegions,col=LOHGRegions$col,border=NULL,avoid.overlapping = F,r0=0,r1=1)
dev.off()
