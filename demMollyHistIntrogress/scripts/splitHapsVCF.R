#### LOAD PACKAGES ####
library(vcfR)

#### FUNCTION ####
#split sample
splitSamp <- function(sampleColumn){
  
  #list for haplotypes
  haps <- list()
  
  #split sample
  haps[[1]] <- sapply(strsplit(sampleColumn,":"),function(x) paste(c(strsplit(x[1],
                                                                              split=c("[/\\|]"))[[1]][1],
                                                              x[2:length(x)]),collapse=":"))
  haps[[2]] <- sapply(strsplit(sampleColumn,":"),function(x) paste(c(strsplit(x[1],
                                                                              split=c("[/\\|]"))[[1]][2],
                                                                     x[2:length(x)]),collapse=":"))
  
  #return both haplotypes
  return(haps)
}


#### LOAD IN VCF ####
#take in command line arguments
args <- commandArgs(trailingOnly = TRUE)
vcfPath <- args[1]
splitVCFPath <- args[2]

#read vcf
#dat <-  read.vcfR(paste0("data/haps/chrSymParHybCombinedFilt.vcf.gz"))[1:100000,]
dat <-  read.vcfR(vcfPath)
hybHapMap <- read.delim(paste0("data/hybHapMap.tsv"),
                  sep="\t",row.names = NULL,header = T)

#### SPLIT HYBRID SAMPLES #### 

#list to store
splitHybGts <- list()

#loop through hybrid samples
# MAKE SURE TO CHANGE #
for(i in 2:20){
  
  #get samp
  samp <- colnames(dat@gt)[i]
  
  #print message
  print(paste0("processing ",samp))
  
  #split hybrid genotypes
  curSplitHybGts <- splitSamp(dat@gt[,i])
  
  #add to overall list
  splitHybGts[[length(splitHybGts) + 1]] <- curSplitHybGts[[1]]
  names(splitHybGts)[[length(splitHybGts)]] <- paste0(samp,"_hap",hybHapMap[which(hybHapMap$SAMP == samp),2])
  splitHybGts[[length(splitHybGts) + 1]] <- curSplitHybGts[[2]]
  names(splitHybGts)[[length(splitHybGts)]] <- paste0(samp,"_hap",hybHapMap[which(hybHapMap$SAMP == samp),3])

}

#bind split haplotypes to par samples
newGTs <- as.matrix(cbind(dat@gt[,1],
                          do.call(cbind,splitHybGts),
                          dat@gt[,21:26]))
colnames(newGTs)[1] <- "FORMAT"

#replace gts
datNew <- dat
datNew@gt <- as.matrix(newGTs)

#### SAVE VCF ####
#save
write.vcf(datNew,splitVCFPath)

#unzip
#system(paste0("gunzip data/splitHaps/",chr,"PhasedSplit.vcf.gz")) 






