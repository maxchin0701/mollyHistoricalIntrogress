# This script summarizes copy number variation across simulated and empirical 
# LOH tracts. 
# output:
#   LOHSummaryCN.tsv: tsv with two columns: category of tract (simulated vs 
#                     empirical) and hpdInterval/mean CNV across tract


#### LIBRARIES ####
library(coda)

#### SET SEED ####
set.seed(1)

#### READ IN PARAMS ####
args <- commandArgs(trailingOnly = TRUE)
inLOHRegions <- args[1]
inGCNV <- args[2]
outSummary <- args[3]

#### LOAD IN DATA ####
#load in empirical data
sharedLOHEmpirical <- read.delim(inLOHRegions,
                                 sep='\t',header=F)
colnames(sharedLOHEmpirical) <- c("chr","start","end")
sharedLOHEmpirical$intSize <- sharedLOHEmpirical$end - sharedLOHEmpirical$start

#Load in copy number variation
genomeCN <- read.delim(inGCNV,
                       sep='\t',header=T)
colnames(genomeCN)[5:23] <- gsub(".","-",colnames(genomeCN)[5:23],
                                 fixed=T)

#### PROCESS GCNV DF ####
#calculate average across samples
sampWideAvgL<- list()
for(i in 1:nrow(genomeCN)){
  #print(paste0("interval ",i))
  sampWideAvgL[[length(sampWideAvgL) + 1]] <- mean(as.numeric(genomeCN[i,5:23]))
}

#append to df
genomeCN$sampWideAvg <- do.call(rbind,sampWideAvgL)
colnames(genomeCN)[24] <- c("sampAvg")

#### PROCESS SIM DATA ####

#set up list to store
simLOHCNV <- list()

#for loop through sim data
for(i in 1:1000){
  #print message
  print(paste0("processing iteration ",i))
  simIter <- i
  
  #read in sim regions
  LOHSim <- read.delim(paste0("output/LOHRegionsSharedSim/simLOHShared_",simIter,".bed"),
                                   sep='\t',header=F)
  
  #sample tract
  tract <- sample(1:nrow(LOHSim),1)
  
  #subset CNV
  curSimCN <- genomeCN[which(genomeCN$chr==LOHSim[tract,1] &
            ((genomeCN$start<=LOHSim[tract,2] &
               genomeCN$end>=LOHSim[tract,2]) |
            (genomeCN$start>=LOHSim[tract,2] &
               genomeCN$end<=LOHSim[tract,3]) |
            (genomeCN$start<=LOHSim[tract,3] &
               genomeCN$end>=LOHSim[tract,3]))),]
  
  #calculate average across region
  simLOHCNV[[length(simLOHCNV)+1]] <- mean(curSimCN$sampAvg)
  
  remove(LOHSim)
}

hpdSimCNV <- HPDinterval(as.mcmc(unlist(simLOHCNV)))

#### PROCESS EMPIRICAL DATA ####
#object to store
empLOHCNV <- list()

#loop
for(i in 1:nrow(sharedLOHEmpirical)){
  #print message
  print(paste0("processing empirical tract ",i))
  
  #subset CNV
  curEmpCN <- genomeCN[which(genomeCN$chr==sharedLOHEmpirical[i,1] &
                               ((genomeCN$start<=sharedLOHEmpirical[i,2] &
                                   genomeCN$end>=sharedLOHEmpirical[i,2]) |
                                  (genomeCN$start>=sharedLOHEmpirical[i,2] &
                                     genomeCN$end<=sharedLOHEmpirical[i,3]) |
                                  (genomeCN$start<=sharedLOHEmpirical[i,3] &
                                     genomeCN$end>=sharedLOHEmpirical[i,3]))),]
  
  #calculate average across region
  empLOHCNV[[length(empLOHCNV)+1]] <- mean(curEmpCN$sampAvg)
}

#### SAVE ####
LOHCNSummary <- as.data.frame(matrix(c("sim",rep("emp",65),
         paste0(hpdSimCNV[1]," - ",hpdSimCNV[2]), unlist(empLOHCNV)),ncol=2,nrow=66))

colnames(propSigDf) <- c("category", "tractAveragedCN")
#save
write.table(LOHCNSummary, 
            file=outSummary, 
            quote=FALSE, sep='\t',row.names=FALSE)













