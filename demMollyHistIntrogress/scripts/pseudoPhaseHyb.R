#### LOAD PACKAGES ####
library(vcfR)

#### FUNCTIONS ####
#replace PS if preexisting
replacePhaseSet <- function(gt,hapStart){
  curGT <- paste(strsplit(gt,
                          split=":")[[1]][-length(strsplit(gt,
                                                           split=":")[[1]])],collapse=":")
  
  #enter new GT
  return(paste(curGT,hapStart,sep=':'))
  
}

#new PS if not preexisting
createPhaseSet <- function(gt,hapStart){
  #enter new GT
  return(paste(gt,hapStart,sep=':'))
}

#phase if refAll=mexAll
phaseMatchedRef <- function(gt){
  curGTString <- gt
  curGT <- strsplit(gt,split=":")[[1]][1]
  if(curGT=="1|0"){
    newGT <- "0|1"
  } else {
    newGT <- gsub("/","|",curGT)
  }
  
  return(paste(c(newGT,
                 strsplit(gt,split=":")[[1]][-1]),
               collapse = ':'))
  
}

#phase if refAll!=mexAll
phaseMismatchedRef <- function(gt){
  curGTString <- gt
  curGT <- strsplit(gt,split=":")[[1]][1]
  if(curGT=="0|0" ||
     curGT=="0/0"){
    newGT <- "1|1"
  } else if(curGT=="1|1" ||
            curGT=="1/1"){
    newGT <- "0|0"
  } else if(curGT=="1|0"){
    newGT <- "0|1"
  } else {
    newGT <- gsub("/","|",curGT)
  }
  
  return(paste(c(newGT,
                 strsplit(gt,split=":")[[1]][-1]),
               collapse = ':'))
  
}


#### LOAD IN DATA ####

#take in command line arguments
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]

#read vcf
dat <-  read.vcfR(paste0("data/hapCallParHybs/",chr,".vcf.gz"))
#dat <-  read.vcfR("../data/hapCalls/chr1.vcf.gz")[1:500000,]

#define parental samples
plat <- 20:25
pmex <- 26:30

#### PULL OUT GENOTYPES ####
gtsRaw <- extract.gt(dat,element="GT",return.alleles = T)

#### IDENTIFY FIXED DIFFERENCES ####
#list for fixed sites
fixedDiff <- list()

#loop through loci
for(i in 1:nrow(gtsRaw)){
  #get latipinna characters at locus
  locCharSp1 <- unique(unlist(strsplit(as.vector(gtsRaw[i,plat]),split=c("[/\\|]"))))
  #get latipinna allele(s)
  allSp1 <- locCharSp1[which(locCharSp1 
                             %in% c("A","C","T","G"))]
  
  #get Mexicana characters at locus
  locCharSp2 <- unique(unlist(strsplit(as.vector(gtsRaw[i,pmex]),split=c("[/\\|]"))))
  #get mexicana allele(s)
  allSp2 <- locCharSp2[which(locCharSp2 
                             %in% c("A","C","T","G"))]
  
  #get formosa characters at locus
  locCharSp3 <- unique(unlist(strsplit(as.vector(gtsRaw[i,-c(pmex,plat)]),split=c("[/\\|]"))))
  #get mexicana allele(s)
  allSp3 <- locCharSp3[which(locCharSp3 
                             %in% c("A","C","T","G"))]
  
  #check for fixed diffs
  if(length(allSp1) == 1 &&
     length(allSp2) == 1 &&
     all(allSp3 %in% c(allSp1,allSp2)) &&
     length(allSp3) == 2){
    if(allSp1 != allSp2){
      fixedDiff[[length(fixedDiff)+1]] <- i
    }
  }
}

#subset loci
fixedDiffDat <- dat[unlist(fixedDiff),]
remove(dat)
fixedDiffGts <- gtsRaw[unlist(fixedDiff),]
remove(gtsRaw)

#### PHASE VCF ####

#first locus
hapStart <- strsplit(rownames(fixedDiffGts),split="_")[[1]][2]

#lists to store
refAllList <- list()
altAllList <- list()
gtList <- list()

#subset filtered vcf
for(j in 1:nrow(fixedDiffDat)){
  print(paste0("processing locus ",j))
  #get ref,parental alleles
  refAll <- fixedDiffDat@fix[j,4]
  pmexAll <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,pmex]),
                                    split=c("[/\\|]"))))
  platAll <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,plat]),
                                    split=c("[/\\|]"))))

  #assign reference allele to pmex and alt allele to plat
  refAllList[[j]] <- pmexAll[which(pmexAll
                                         %in% c("A","C","T","G"))]
  altAllList[[j]] <- platAll[which(platAll
                                         %in% c("A","C","T","G"))]

  #if the previous reference allele was mexicana
  if(refAll == fixedDiffDat@fix[j,4]){
    #phase
    newGts <- sapply(fixedDiffDat@gt[j,2:ncol(fixedDiffDat@gt)],
                                                         phaseMatchedRef)
  } else {
    #phase
    newGts <- sapply(fixedDiffDat@gt[j,2:ncol(fixedDiffDat@gt)],
                                                         phaseMismatchedRef)
  }

  #add in phase set data
  if(!("PS" %in% as.vector(strsplit(fixedDiffDat@gt[j,1],
                                  split=":"))[[1]])){
    #function to phase set of
    gtList[[j]] <- c(paste0(fixedDiffDat@gt[j,1],":PS"),
                sapply(newGts,
                       createPhaseSet,hapStart=hapStart))
  } else {
    #print("updating phase")
    gtList[[j]] <- c(fixedDiffDat@gt[j,1],
                     sapply(newGts,
                            replacePhaseSet,hapStart=hapStart))
  }
}

#modify gt object
fixedDiffDat@fix[,4] <- do.call(rbind,refAllList)
fixedDiffDat@fix[,5] <- do.call(rbind,altAllList)
fixedDiffDat@gt <- do.call(rbind,gtList)
colnames(fixedDiffDat@gt)[1] <- "FORMAT"

#### WRITE VCF FILES ####
write.vcf(fixedDiffDat,file=paste0("data/refPanel/",chr,"WhatshapRefPanel.vcf.gz"))


