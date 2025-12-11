#### LOAD PACKAGES ####
library(stringr)

#### READ IN PARAMS ####
outFile <- snakemake@output[[1]]
#outFile <- "data/gCNV/gCNV/processedGCNV.tsv"

#### OBJECT TO STORE ####
genomeCN <- list()

#### LOOP THROUGH CHRS ####
#chrDirs <- list.dirs(path=paste0("../data/gCNV/gCNV"),recursive=F)
chrDirs <- list.dirs(path="data/gCNV/gCNV",recursive=F)
chrDirs <- grep("calls",chrDirs,value=T)
chrs <- matrix(unlist(strsplit(chrDirs,"V_|[s-]+")),ncol=3,byrow=T)[,2]

#loop through chrs
for(i in 1:length(chrs)){
  #current chr and dir
  curChr <- chrs[i]
  curChrDir <- chrDirs[i]
  
  ##### READ IN CHR BED FILE ####
  # chrIntervals <- read.delim(paste0("../data/refGenome/intervalLists/PFor_",curChr,".bed"),
  #                         sep="\t",row.names = NULL,header = F)
  chrIntervals <- read.delim(paste0("data/refGenome/intervalLists/PFor_",curChr,".bed"),
                             sep="\t",row.names = NULL,header = F)
  colnames(chrIntervals) <- c("chr","start","end")
  chrIntervals$intMean <- chrIntervals$start + (chrIntervals$end - chrIntervals$start)/2
  
  #### COLLECT SAMPLE NAMES ####
  samps <- c()
  sampDirs <- list.dirs(path=curChrDir,recursive=F)

  for(j in sampDirs){
    samps <- c(samps,
               names(read.delim(paste0(j,"/sample_name.txt"))))
  }
  print(samps)
  #replace dots with dash
  samps <- gsub(".","-",samps,fixed=T)
  print(samps)
  print(i)
  #extract pfor samples
  pForIndices <- which(matrix(unlist(strsplit(samps,split="gi|ve|[_]+")),
                              ncol=2,byrow=T)[,1] == "Pfo")
  
  #list for chromosome CN
  chrCN <- list()
  chrCN[[1]] <- chrIntervals
  
  #### LOOP THROUGH SAMPLES ####
  for(j in pForIndices){
    curSampDir <- sampDirs[j] 
    #read in raw copy counts
    muDCR <- read.delim(paste0(curSampDir,"/mu_denoised_copy_ratio_t.tsv"),sep="\t",row.names = NULL)
    muDCR <- as.numeric(muDCR[4:nrow(muDCR),1])
    
    chrCN[[length(chrCN)+1]] <- muDCR
    names(chrCN)[[length(chrCN)]] <- samps[j]
  }
  
  #### ADD TO GENOME WIDE LIST ####
  #bind to genomeCN
  genomeCN[[length(genomeCN)+1]] <- do.call(cbind,chrCN)
  names(genomeCN)[[length(genomeCN)]] <- curChr
}

#### RBIND AND RENAME COLUMNS ####
genomeCNDf <- do.call(rbind,
                      genomeCN)
colnames(genomeCNDf)[1:4] <- c("chr","start","end","mean")

#### SAVE ####
write.table(genomeCNDf, file=outFile, 
            quote=FALSE, sep='\t',row.names=FALSE)









