#### PREPARE OBJECTS ####

#set genome size, distribution to simulate
chrSizes <- read.delim(paste0("data/chrIndex.tsv"),
                       sep="\t",row.names = NULL,header = T)
#load in data
sharedLOHEmpirical <- read.delim("output/LOHRegionsSharedCombined/LOHRegionsSharedCombined.bed",
                                 sep='\t',header=F)
colnames(sharedLOHEmpirical) <- c("chr","start","end")
sharedLOHEmpirical$intSize <- sharedLOHEmpirical$end - sharedLOHEmpirical$start

LOHSharedIntMean <- mean(sharedLOHEmpirical$intSize)
LOHSharedIntTotal <- sum(sharedLOHEmpirical$intSize)
nReps<-1000
simLOH <- list()

#### SIMULATE INTERVALS ####
#loop through replicates
for(i in 1:nReps){
  simLOHList <- list()
  curTotLength <- 0
  while(curTotLength <= LOHSharedIntTotal){
    chr <- sample(size=1,chrSizes[,1])
    intSize <- round(rexp(1,1/LOHSharedIntMean))
    intStart <- sample(1:(chrSizes[which(chrSizes[,1]==chr),3] - intSize),1)
    simLOHList[[length(simLOHList)+1]] <- cbind(chr,intStart,intStart+intSize)
    curTotLength <- curTotLength+intSize
  }
  simLOHDf <- do.call(rbind,simLOHList)
  simLOH[[length(simLOH)+1]] <- simLOHDf
}

#save
for(i in 1:length(simLOH)){
  write.table(simLOH[[i]],
              file=paste0("output/LOHRegionsSharedSim/simLOHShared_",i,".bed"),
              quote=FALSE, sep='\t',row.names=F,col.names = F)
}













