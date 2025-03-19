#### LOAD PACKAGES ####
library(vcfR)
library(coda)

#### READ IN DATA ####
chr <- snakemake@wildcards[["chr"]]
#chr <- "chr7_24"
parSpecies <- snakemake@wildcards[["pars"]]
#parSpecies <- "LM"

#print message
print(paste0("Processing ",chr))

#### READ IN AND PRUNE CHROMOSOME VCF ###
dat <-  read.vcfR(snakemake@input[[1]])
#dat <- read.vcfR(paste0("../data/hapCalls/",chr,"CombinedFilt.vcf.gz"))

#### DEFINE PARENTAL SPECIES SAMPS####
if(parSpecies == "LM"){
  dat <- dat[,-c(21,31)]
  plat <- 20:24
  pmex <- 25:28
} else {
  plat <- 20:25
  pmex <- 26:30
}

#extract allele depths and genotype qualities
sampDepths <- extract.gt(dat,element="DP")
genQual <- extract.gt(dat,element="GQ")
unfiltGts <- extract.gt(dat,element="GT",return.alleles = T)

#define cutoffs
cutoffGQ <- as.numeric(snakemake@wildcards[["gqcut"]])
#cutoffGQ <- 20
#cutoffDP <- 8
cutoffDP <- as.numeric(snakemake@wildcards[["dpcut"]])
#cutoffDP <- 10
#maxDP

#print message
print(paste0("Filtering SNPs based on depth and genotype quality. Quality cutoff: ", cutoffGQ," Depth cutoff: ", cutoffDP))

#setup variable for loci to remove
dropSNPs <- list()

#loop through samp cols
for(j in 1:ncol(sampDepths)){
  dropSNPs[[length(dropSNPs)+1]] <- which(as.numeric(sampDepths[,j])<cutoffDP)
  dropSNPs[[length(dropSNPs)+1]] <- which(as.numeric(genQual[,j])<cutoffGQ)
  dropSNPs[[length(dropSNPs)+1]] <- which(unfiltGts[,j] == ".")
  dropSNPs[[length(dropSNPs)+1]] <- which(as.numeric(sampDepths[,j])>4*mean(as.numeric(sampDepths[,j]),na.rm=T))
}

#get final list of snps to drop
dropSNPs <- unique(unlist(dropSNPs))

if(nrow(sampDepths) - length(dropSNPs) == 0){
  #print message
  print("No loci retained after GQ and DP filtering. Saving empty files")
  
  file.create(file=paste0("output/LOHRegions/",chr,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
  file.create(file=paste0("output/LOHAnc/",chr,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
  file.create(file=paste0("output/SNPall/",chr,"SNPall_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
  file.create(file=paste0("output/fixDiffVCF/",chr,"FixDiff_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".vcf.gz"))
  
} else {
  #print message
  print(paste0(nrow(sampDepths) - length(dropSNPs), " loci kept after depth, missing data, and GQ filtering"))
  
  #remove large objects
  rm(sampDepths,genQual,unfiltGts)
  
  #drop snps from vcf file
  datFilt <- dat[-dropSNPs,]
  
  #remove unnecessary objects
  rm(dropSNPs,dat)
  
  #### DETERMINE FIXED HETEROZYGOUS SITES ####
  
  #print message
  print(paste0("Determining fixed differences between parental species"))
  
  #get genotypes from datFilt
  gts <- extract.gt(datFilt,element="GT",return.alleles=T)
  #sampDepthsFilt <- extract.gt(datFilt,element="DP")
  #genQualFilt <- extract.gt(datFilt,element="GQ")
  
  #list for fixed sites
  fixedDiff <- list()
  
  #loop through loci
  for(j in 1:nrow(gts)){
    #get latipinna characters at locus
    locCharSp1 <- unique(unlist(strsplit(as.vector(gts[j,plat]),split=c("|"))))
    #get latipinna allele(s)
    allSp1 <- locCharSp1[which(locCharSp1 
                 %in% c("A","C","T","G"))]
    
    #get Mexicana characters at locus
    locCharSp2 <- unique(unlist(strsplit(as.vector(gts[j,pmex]),split=c("|"))))
    #get mexicana allele(s)
    allSp2 <- locCharSp2[which(locCharSp2 
                            %in% c("A","C","T","G"))]
    
    #check for fixed diffs
    if(length(allSp1) == 1 &&
       length(allSp2) == 1){
      if(allSp1 != allSp2){
        fixedDiff[[length(fixedDiff)+1]] <- j
      }
    }
  }
  
  if(length(fixedDiff) == 0){
    #print message
    print("No fixed differences called after GQ and DP filtering. Saving empty files")
    
    file.create(file=paste0("output/LOHRegions/",chr,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
    file.create(file=paste0("output/LOHAnc/",chr,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
    file.create(file=paste0("output/SNPall/",chr,"SNPall_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
    file.create(file=paste0("output/fixDiffVCF/",chr,"FixDiff_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".vcf.gz"))
    
  } else {
    #print message
    print(paste0(length(fixedDiff), " loci fixed different between parental species"))
    
    #### LOOP THROUGH AND IDENTIFY LOH EVENTS ####
    print(paste0("Identifying LOH events in P. formosa samples"))
    
    #subset loci to just fixed diffs
    if("matrix" %in% class(gts[unlist(fixedDiff),])){
      #fixedDiffDp <- apply(sampDepthsFilt[unlist(fixedDiff),],2,as.numeric)
      #fixedDiffGQ <- apply(genQualFilt[unlist(fixedDiff),],2,as.numeric)
      fixedDiffGts <- gts[unlist(fixedDiff),]
    } else {
      #fixedDiffDp <- as.numeric(sampDepthsFilt[unlist(fixedDiff),])
      #fixedDiffGQ <- as.numeric(genQualFilt[unlist(fixedDiff),])
      fixedDiffGts <- matrix(gts[unlist(fixedDiff),],ncol=ncol(gts))
      rownames(fixedDiffGts) <- rownames(gts)[unlist(fixedDiff)]
      colnames(fixedDiffGts) <- colnames(gts)
    }
    #initialize objects to store
    start <- list()
    end <- list()
    nSamps <- list()
    nSNPs <- list()
    samps <- list()
    startAnc <- list()
    endAnc <- list()
    anc <- list()
    startAll <- list()
    endAll <- list()
    
    
    #boolean for if prev site also had LOH
    prevLocLOH <- FALSE
    prevSamps <- c()
    consecLOH <- 0
    
    for(j in 1:nrow(fixedDiffGts)){
      
      #print message
      print(paste0("processing locus ", j))
      
      #get parental alleles
      allSp1 <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,plat]),
                                                        split=c("|"))))[1]
      allSp2 <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,pmex]),
                                       split=c("|"))))[1]
      
      #extract position
      curPos <- as.numeric(strsplit(rownames(fixedDiffGts),"_")[[j]][
        length(strsplit(rownames(fixedDiffGts),"_")[[j]])])
      
      #add position to start and end
      startAll[[length(startAll)+1]] <- curPos - 1
      endAll[[length(endAll)+1]] <- curPos + 1
      
      #extract pfor genotypes
      pforGTs <- strsplit(as.vector(fixedDiffGts[j,-c(pmex,plat)]),split=c("|"))
      
      #variables for storing LOH samps and ancestry
      sampLOH <- list()
      ancLOH <- list()
      
      #loop through pfor samples
      for(k in 1:length(pforGTs)){
        #check if LOH has occurred for sample
        if(length(unique(pforGTs[[k]])) == 2){

          #add to sample LOH list
          sampLOH[[length(sampLOH) + 1]] <- colnames(fixedDiffGts)[k]
          
          #get ancestry
          if(unique(pforGTs[[k]])[1] == allSp1){
            ancLOH[[length(ancLOH) + 1]] <- "plat"
          } else {
            ancLOH[[length(ancLOH) + 1]] <- "pmex"
          }
        }
      }
      
      #skip rest if no samples have LOH
      if(length(sampLOH) == 0){
        prevLocLOH <- FALSE
        next
      }
      
      #summarize ancestry across samples
      if(length(unique(ancLOH)) == 1){
        consensusAnc <- unique(ancLOH)
      } else {
        consensusAnc <- list("mix")
      }
      
      #add to lists
      if(prevLocLOH == T &&
         (length(intersect(sampLOH,prevSamps)) == length(sampLOH) &&
          length(intersect(sampLOH,prevSamps)) == length(prevSamps)) &&
         unlist(consensusAnc) == unlist(anc[[length(anc)]]) &&
         consecLOH >= 1){
        
        #check if this is initiating a block
        if(consecLOH == 1){
          start <- start[-length(start)]
          end <- end[-length(end)]
          nSamps <- nSamps[-length(nSamps)]
          nSNPs <- nSNPs[-length(nSNPs)]
          nSNPs[[length(nSNPs)]] <- 2
          samps <- samps[-length(samps)]
          samps[[length(samps)]] <- paste(sampLOH,collapse=",")
          startAnc <- startAnc[-length(startAnc)]
          endAnc <- endAnc[-length(endAnc)]
          anc <- anc[-length(anc)]
        }
        
        #add to end of most recent LOH region
        end[[length(end)]] <- curPos + 1
        nSNPs[[length(nSNPs)]] <- nSNPs[[length(nSNPs)]] + 1
        
        #add to end of most recent anc region
        endAnc[[length(endAnc)]] <- curPos + 1
        
        #iterate consecLOH
        consecLOH <- consecLOH + 1
        
      # } else if (prevLocLOH == T &&
      #            (length(intersect(sampLOH,prevSamps)) == length(sampLOH) &&
      #              length(intersect(sampLOH,prevSamps)) == length(prevSamps)) &&
      #            unlist(consensusAnc) != unlist(anc[[length(anc)]])){
      # 
      #   #add to end of most recent LOH region
      #   end[[length(end)]] <- curPos + 1
      # 
      #   #start new anc region
      #   startAnc[[length(startAnc)+1]] <- endAnc[[length(endAnc)]] + 1
      #   endAnc[[length(endAnc)+1]] <- curPos + 1
      #   anc[[length(anc)+1]] <- unlist(consensusAnc)
      # 
      # } else if (prevLocLOH == T &&
      #            (length(intersect(sampLOH,prevSamps)) != 0) &&
      #            unlist(consensusAnc) == unlist(anc[[length(anc)]])){
      #   
      #   #start new LOH region immediately following previous one
      #   start[[length(start)+1]] <- end[[length(end)]] + 1
      #   end[[length(end)+1]] <- curPos + 1
      #   nSamps[[length(nSamps)+1]] <- length(sampLOH)
      #   
      #   #add to end of most recent anc region
      #   endAnc[[length(endAnc)]] <- curPos + 1
      #   
      # } else if (prevLocLOH == T &&
      #            (length(intersect(sampLOH,prevSamps)) != 0) &&
      #            unlist(consensusAnc) != unlist(anc[[length(anc)]])){
      #   
      #   #start new LOH region immediately following previous one
      #   start[[length(start)+1]] <- end[[length(end)]] + 1
      #   end[[length(end)+1]] <- curPos + 1
      #   nSamps[[length(nSamps)+1]] <- length(sampLOH)
      #   
      #   #start new anc region
      #   startAnc[[length(startAnc)+1]] <- endAnc[[length(endAnc)]] + 1
      #   endAnc[[length(endAnc)+1]] <- curPos + 1
      #   anc[[length(anc)+1]] <- unlist(consensusAnc)
        
      } else {
        
        if(prevLocLOH == T &&
           (length(intersect(sampLOH,prevSamps)) == length(sampLOH) &&
            length(intersect(sampLOH,prevSamps)) == length(prevSamps)) &&
           unlist(consensusAnc) == unlist(anc[[length(anc)]])){
          consecLOH <- consecLOH + 1
        } else {
          consecLOH <- 0
        }
        
        #start new LOH region from curPos
        start[[length(start)+1]] <- curPos - 1
        end[[length(end)+1]] <- curPos + 1
        nSamps[[length(nSamps)+1]] <- length(sampLOH)
        nSNPs[[length(nSNPs)+1]] <- 1
        samps[[length(samps)+1]] <- paste(sampLOH,collapse=",")
        
        #start new anc region
        startAnc[[length(startAnc)+1]] <- curPos - 1
        endAnc[[length(endAnc)+1]] <- curPos + 1
        anc[[length(anc)+1]] <- unlist(consensusAnc)
      }
      
      #update
      prevLocLOH <- TRUE
      prevSamps <- sampLOH
      
    }
    
    #### PHASE VCF ####
    
    #subset filtered vcf
    fixedDiffVCF <- datFilt[unlist(fixedDiff),]
    for(j in 1:nrow(fixedDiffVCF)){
      refAll <- fixedDiffVCF@fix[j,4]
      pmexAll <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,pmex]),
                                        split=c("|"))))[1]
      
      if(refAll != pmexAll){
        for(k in 2:29){
          curGTString <- fixedDiffVCF@gt[j,k]
          curGT <- strsplit(fixedDiffVCF@gt[j,k],split=":")[[1]][1]
          if(curGT=="0/0"){
            newGT <- "1/1"
          } else if(curGT=="0|0"){
            newGT <- "1|1"
          } else if(curGT=="1/1"){
            newGT <- "0/0"
          } else if(curGT=="1|1"){
            newGT <- "0|0"
          } else {
            newGT <- curGT
          }
          
          fixedDiffVCF@gt[j,k] <- paste(c(newGT,
                                          strsplit(fixedDiffVCF@gt[j,k],split=":")[[1]][-1]), 
                                        collapse = ':')
        }
      }
      fixedDiffVCF@fix[j,4] <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,pmex]),
                                                      split=c("|"))))[1]
      fixedDiffVCF@fix[j,5] <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,plat]),
                                                      split=c("|"))))[1]
      
      
    }
    
    nMismatch <- 0
    for(j in 1:nrow(fixedDiffVCF)){
      refAll <- fixedDiffVCF@fix[j,4]
      pmexAll <- unique(unlist(strsplit(as.vector(fixedDiffGts[j,pmex]),
                             split=c("|"))))[1]
      print(paste0("locus ",j,": reference allele is ",
                   refAll,", pmex allele is ", pmexAll))
      if(refAll != pmexAll){
        nMismatch <- nMismatch + 1
      }
    }
    
    #print message
    print("Done processing VCF. Saving region and ancestry files now.")
    
    #### CHECK IF EMPTY ####
    if(length(start) == 0){
      #print message
      print("No LOH regions detected after GQ and DP filtering. Saving empty files")
      
      #### COMBINE LISTS ####
      SNPall <- as.data.frame(cbind(chr,do.call(rbind,startAll),
                                    do.call(rbind,endAll)))
      
      file.create(file=paste0("output/LOHRegions/",chr,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
      file.create(file=paste0("output/LOHAnc/",chr,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"))
      write.table(SNPall, file=paste0("output/SNPall/",chr,"SNPall_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"), 
                  quote=FALSE, sep='\t',row.names=FALSE)
      write.vcf(datFilt[unlist(fixedDiff),], file=paste0("output/fixDiffVCF/",chr,"FixDiff_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".vcf.gz"))
      
    } else {
      #### COMBINE LISTS ####
      SNPall <- as.data.frame(cbind(chr,do.call(rbind,startAll),
                                    do.call(rbind,endAll)))
      
      LOHRegions <- as.data.frame(cbind(chr,do.call(rbind,start),
                                  do.call(rbind,end),
                                  do.call(rbind,nSamps),
                                  do.call(rbind,nSNPs),
                                  do.call(rbind,samps)))
      
      LOHAnc <- as.data.frame(cbind(chr,
                                    do.call(rbind,startAnc),
                                    do.call(rbind,endAnc),
                                    do.call(rbind,anc)))
      
      #### PREPARE FOR SAVING ####
      colnames(SNPall) <- c("chr","start","end")
      colnames(LOHRegions) <- c("chr","start","end","group","nSNPS","samps")
      colnames(LOHAnc) <- c("chr","start","end","group")
      
      #converting columns
      SNPall$start <- as.numeric(SNPall$start)
      SNPall$end <- as.numeric(SNPall$end)
      LOHRegions$start <- as.numeric(LOHRegions$start)
      LOHRegions$end <- as.numeric(LOHRegions$end)
      LOHRegions$group <- as.numeric(LOHRegions$group)
      LOHRegions$nSNPS <- as.numeric(LOHRegions$nSNPS)
      LOHRegions$samps <- as.numeric(LOHRegions$samps)
      LOHAnc$start <- as.numeric(LOHAnc$start)
      LOHAnc$end <- as.numeric(LOHAnc$end)
      
      #### SAVE ####
      write.table(LOHRegions, file=paste0("output/LOHRegions/",chr,"LOHRegions_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"), 
                  quote=FALSE, sep='\t',row.names=FALSE)
      write.table(LOHAnc, file=paste0("output/LOHAnc/",chr,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"), 
                  quote=FALSE, sep='\t',row.names=FALSE)
      write.table(SNPall, file=paste0("output/SNPall/",chr,"SNPall_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".tsv"), 
                  quote=FALSE, sep='\t',row.names=FALSE)
      write.vcf(fixedDiffVCF, file=paste0("output/fixDiffVCF/",chr,"FixDiff_GQ",cutoffGQ,"_DP",cutoffDP,"_",parSpecies,".vcf.gz"))
    }
  }
}





