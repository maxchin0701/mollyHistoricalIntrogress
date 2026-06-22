#### LOAD PACKAGES ####
library(twosamples)

#empirical exon
datExonEmp <- read.delim(paste0("output/selectionStats/LOHRegionsEmpStatsExon.tsv"),
                         sep="\t",row.names = NULL,header = T)
datExonEmp$iter <- 0
datExonEmp$class <- "Empirical"

#sim
datSimExonList <- list()
for(i in 1:1000){
  #exon
  curDatExonSim <- read.delim(paste0("output/selectionStats/simStatsOut/simRegionsStatsExon_",i,".tsv"),
                              sep="\t",row.names = NULL,header = T)
  curDatExonSim$iter <- i
  datSimExonList[[i]] <- curDatExonSim
}

#bind to form df
datExonSim <- do.call(rbind,datSimExonList)
datExonSim$class <- "Simulated"

#combine empirical and simulated
datExonAll <- rbind(datExonEmp,datExonSim)

#exon
#non exon
sigDiffKSExonPlatD <- list()
sigDiffMWUExonPlatD <- list()
sigDiffCVMExonPlatD <- list()
sigDiffKSExonPmexD <- list()
sigDiffMWUExonPmexD <- list()
sigDiffCVMExonPmexD <- list()

sigDiffKSExonPlatH <- list()
sigDiffMWUExonPlatH <- list()
sigDiffCVMExonPlatH <- list()
sigDiffKSExonPmexH <- list()
sigDiffMWUExonPmexH <- list()
sigDiffCVMExonPmexH <- list()

sigDiffKSExonPlatE <- list()
sigDiffMWUExonPlatE <- list()
sigDiffCVMExonPlatE <- list()
sigDiffKSExonPmexE <- list()
sigDiffMWUExonPmexE <- list()
sigDiffCVMExonPmexE <- list()


#loop through simulations
for(i in 1:1000){
  #get simulated and empirical data (D)
  simDatLatD <- datExonAll[which(datExonAll$class == "Simulated" & 
                                      datExonAll$iter == i &
                                      datExonAll$pop == "Plat"),8]
  empDatLatD <- datExonAll[which(datExonAll$class == "Empirical" & 
                                      datExonAll$pop == "Plat"),8]
  simDatMexD <- datExonAll[which(datExonAll$class == "Simulated" & 
                                      datExonAll$iter == i &
                                      datExonAll$pop == "Pmex"),8]
  empDatMexD <- datExonAll[which(datExonAll$class == "Empirical" & 
                                      datExonAll$pop == "Pmex"),8]
  
  #two sided ks test (D)
  pValKSLatD <- ks.test(empDatLatD,simDatLatD)[2]
  pvalMWULatD <- wilcox.test(empDatLatD,simDatLatD)[3]
  pvalCVMLatD <- cvm_test(empDatLatD,simDatLatD,nboots=1000)[2]
  pValKSMexD <- ks.test(empDatMexD,simDatMexD)[2]
  pvalMWUMexD <- wilcox.test(empDatMexD,simDatMexD)[3]
  pvalCVMMexD <- cvm_test(empDatMexD,simDatMexD,nboots=1000)[2]
  sigDiffKSExonPlatD[[length(sigDiffKSExonPlatD) + 1]] <- pValKSLatD
  sigDiffMWUExonPlatD[[length(sigDiffMWUExonPlatD) + 1]] <- pvalMWULatD
  sigDiffCVMExonPlatD[[length(sigDiffCVMExonPlatD) + 1]] <- pvalCVMLatD
  sigDiffKSExonPmexD[[length(sigDiffKSExonPmexD) + 1]] <- pValKSMexD
  sigDiffMWUExonPmexD[[length(sigDiffMWUExonPmexD) + 1]] <- pvalMWUMexD
  sigDiffCVMExonPmexD[[length(sigDiffCVMExonPmexD) + 1]] <- pvalCVMMexD
  
  #get simulated and empirical data (H)
  simDatLatH <- datExonAll[which(datExonAll$class == "Simulated" & 
                                      datExonAll$iter == i &
                                      datExonAll$pop == "Plat"),12]
  empDatLatH <- datExonAll[which(datExonAll$class == "Empirical" & 
                                      datExonAll$pop == "Plat"),12]
  simDatMexH <- datExonAll[which(datExonAll$class == "Simulated" & 
                                      datExonAll$iter == i &
                                      datExonAll$pop == "Pmex"),12]
  empDatMexH <- datExonAll[which(datExonAll$class == "Empirical" & 
                                      datExonAll$pop == "Pmex"),12]
  
  #two sided ks test (H)
  pValKSLatH <- ks.test(empDatLatH,simDatLatH)[2]
  pvalMWULatH <- wilcox.test(empDatLatH,simDatLatH)[3]
  pvalCVMLatH <- cvm_test(empDatLatH,simDatLatH,nboots=1000)[2]
  pValKSMexH <- ks.test(empDatMexH,simDatMexH)[2]
  pvalMWUMexH <- wilcox.test(empDatMexH,simDatMexH)[3]
  pvalCVMMexH <- cvm_test(empDatMexH,simDatMexH,nboots=1000)[2]
  sigDiffKSExonPlatH[[length(sigDiffKSExonPlatH) + 1]] <- pValKSLatH
  sigDiffMWUExonPlatH[[length(sigDiffMWUExonPlatH) + 1]] <- pvalMWULatH
  sigDiffCVMExonPlatH[[length(sigDiffCVMExonPlatH) + 1]] <- pvalCVMLatH
  sigDiffKSExonPmexH[[length(sigDiffKSExonPmexH) + 1]] <- pValKSMexH
  sigDiffMWUExonPmexH[[length(sigDiffMWUExonPmexH) + 1]] <- pvalMWUMexH
  sigDiffCVMExonPmexH[[length(sigDiffCVMExonPmexH) + 1]] <- pvalCVMMexH
  
  #get simulated and empirical data (E)
  simDatLatE <- datExonAll[which(datExonAll$class == "Simulated" & 
                                      datExonAll$iter == i &
                                      datExonAll$pop == "Plat"),15]
  empDatLatE <- datExonAll[which(datExonAll$class == "Empirical" & 
                                      datExonAll$pop == "Plat"),15]
  simDatMexE <- datExonAll[which(datExonAll$class == "Simulated" & 
                                      datExonAll$iter == i &
                                      datExonAll$pop == "Pmex"),15]
  empDatMexE <- datExonAll[which(datExonAll$class == "Empirical" & 
                                      datExonAll$pop == "Pmex"),15]
  
  #two sided ks test (E)
  pValKSLatE <- ks.test(empDatLatE,simDatLatE)[2]
  pvalMWULatE <- wilcox.test(empDatLatE,simDatLatE)[3]
  pvalCVMLatE <- cvm_test(empDatLatE,simDatLatE,nboots=1000)[2]
  pValKSMexE <- ks.test(empDatMexE,simDatMexE)[2]
  pvalMWUMexE <- wilcox.test(empDatMexE,simDatMexE)[3]
  pvalCVMMexE <- cvm_test(empDatMexE,simDatMexE,nboots=1000)[2]
  sigDiffKSExonPlatE[[length(sigDiffKSExonPlatE) + 1]] <- pValKSLatE
  sigDiffMWUExonPlatE[[length(sigDiffMWUExonPlatE) + 1]] <- pvalMWULatE
  sigDiffCVMExonPlatE[[length(sigDiffCVMExonPlatE) + 1]] <- pvalCVMLatE
  sigDiffKSExonPmexE[[length(sigDiffKSExonPmexE) + 1]] <- pValKSMexE
  sigDiffMWUExonPmexE[[length(sigDiffMWUExonPmexE) + 1]] <- pvalMWUMexE
  sigDiffCVMExonPmexE[[length(sigDiffCVMExonPmexE) + 1]] <- pvalCVMMexE
}


#calculate prop sig nonexon (D)
sigDiffKSExonPlatD <- unlist(sigDiffKSExonPlatD)
sigDiffMWUExonPlatD <- unlist(sigDiffMWUExonPlatD)
sigDiffCVMExonPlatD <- unlist(sigDiffCVMExonPlatD)
sigDiffKSExonPmexD <- unlist(sigDiffKSExonPmexD)
sigDiffMWUExonPmexD <- unlist(sigDiffMWUExonPmexD)
sigDiffCVMExonPmexD <- unlist(sigDiffCVMExonPmexD)
propSigLatDKS <- sum(p.adjust(sigDiffKSExonPlatD,method="BH")<=0.05)/1000
propSigLatDMWU <- sum(p.adjust(sigDiffMWUExonPlatD,method="BH")<=0.05)/1000
propSigLatDCVM <- sum(p.adjust(sigDiffCVMExonPlatD,method="BH")<=0.05)/1000
propSigMexDKS <- sum(p.adjust(sigDiffKSExonPmexD,method="BH")<=0.05)/1000
propSigMexDMWU <- sum(p.adjust(sigDiffMWUExonPmexD,method="BH")<=0.05)/1000
propSigMexDCVM <- sum(p.adjust(sigDiffCVMExonPmexD,method="BH")<=0.05)/1000

#calculate prop sig nonexon (H)
sigDiffKSExonPlatH <- unlist(sigDiffKSExonPlatH)
sigDiffMWUExonPlatH <- unlist(sigDiffMWUExonPlatH)
sigDiffCVMExonPlatH <- unlist(sigDiffCVMExonPlatH)
sigDiffKSExonPmexH <- unlist(sigDiffKSExonPmexH)
sigDiffMWUExonPmexH <- unlist(sigDiffMWUExonPmexH)
sigDiffCVMExonPmexH <- unlist(sigDiffCVMExonPmexH)
propSigLatHKS <- sum(p.adjust(sigDiffKSExonPlatH,method="BH")<=0.05)/1000
propSigLatHMWU <- sum(p.adjust(sigDiffMWUExonPlatH,method="BH")<=0.05)/1000
propSigLatHCVM <- sum(p.adjust(sigDiffCVMExonPlatH,method="BH")<=0.05)/1000
propSigMexHKS <- sum(p.adjust(sigDiffKSExonPmexH,method="BH")<=0.05)/1000
propSigMexHMWU <- sum(p.adjust(sigDiffMWUExonPmexH,method="BH")<=0.05)/1000
propSigMexHCVM <- sum(p.adjust(sigDiffCVMExonPmexH,method="BH")<=0.05)/1000

#calculate prop sig nonexon (E)
sigDiffKSExonPlatE <- unlist(sigDiffKSExonPlatE)
sigDiffMWUExonPlatE <- unlist(sigDiffMWUExonPlatE)
sigDiffCVMExonPlatE <- unlist(sigDiffCVMExonPlatE)
sigDiffKSExonPmexE <- unlist(sigDiffKSExonPmexE)
sigDiffMWUExonPmexE <- unlist(sigDiffMWUExonPmexE)
sigDiffCVMExonPmexE <- unlist(sigDiffCVMExonPmexE)
propSigLatEKS <- sum(p.adjust(sigDiffKSExonPlatE,method="BH")<=0.05)/1000
propSigLatEMWU <- sum(p.adjust(sigDiffMWUExonPlatE,method="BH")<=0.05)/1000
propSigLatECVM <- sum(p.adjust(sigDiffCVMExonPlatE,method="BH")<=0.05)/1000
propSigMexEKS <- sum(p.adjust(sigDiffKSExonPmexE,method="BH")<=0.05)/1000
propSigMexEMWU <- sum(p.adjust(sigDiffMWUExonPmexE,method="BH")<=0.05)/1000
propSigMexECVM <- sum(p.adjust(sigDiffCVMExonPmexE,method="BH")<=0.05)/1000

#### SAVE ####
#D
propSigDfD <- as.data.frame(matrix(c(propSigLatDKS, propSigMexDKS,
                                     propSigLatDMWU, propSigMexDMWU,
                                     propSigLatDCVM, propSigMexDCVM),nrow=2,byrow=F)) 

colnames(propSigDfD) <- c("KS","MWU","CVM")
rownames(propSigDfD) <- c("latExon","mexExon")

#H
propSigDfH <- as.data.frame(matrix(c(propSigLatHKS, propSigMexHKS,
                                     propSigLatHMWU, propSigMexHMWU,
                                     propSigLatHCVM, propSigMexHCVM),nrow=2,byrow=F)) 

colnames(propSigDfH) <- c("KS","MWU","CVM")
rownames(propSigDfH) <- c("latExon","mexExon")

#H
propSigDfE <- as.data.frame(matrix(c(propSigLatEKS, propSigMexEKS,
                                     propSigLatEMWU, propSigMexEMWU,
                                     propSigLatECVM, propSigMexECVM),nrow=2,byrow=F)) 

colnames(propSigDfE) <- c("KS","MWU","CVM")
rownames(propSigDfE) <- c("latExon","mexExon")

#save
write.table(propSigDfD, 
            file="output/results/LOHDStatsExon.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)
write.table(propSigDfH, 
            file="output/results/LOHHStatsExon.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)
write.table(propSigDfE,
            file="output/results/LOHEStatsExon.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)


