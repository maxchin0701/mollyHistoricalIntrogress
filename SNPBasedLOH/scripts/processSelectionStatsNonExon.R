#### LOAD PACKAGES ####
library(twosamples)

##### LOAD IN DATA ####
#empirical nonexon
datNonExonEmp <- read.delim(paste0("output/selectionStats/LOHRegionsEmpStatsNonExon.tsv"),
                   sep="\t",row.names = NULL,header = T)
datNonExonEmp$iter <- 0
datNonExonEmp$class <- "Empirical"

#sim
datSimNonExonList <- list()
for(i in 1:1000){
  #non exon
  curDatNonExonSim <- read.delim(paste0("output/selectionStats/simStatsOut/simRegionsStatsNonExon_",i,".tsv"),
                          sep="\t",row.names = NULL,header = T)
  curDatNonExonSim$iter <- i
  datSimNonExonList[[i]] <- curDatNonExonSim
}

#bind to form df
datNonExonSim <- do.call(rbind,datSimNonExonList)
datNonExonSim$class <- "Simulated"

#combine empirical and simulated
datNonExonAll <- rbind(datNonExonEmp,datNonExonSim)

#### TEST SIGNIFICANCE ####
#non exon
sigDiffKSNonExonPlatD <- list()
sigDiffMWUNonExonPlatD <- list()
sigDiffCVMNonExonPlatD <- list()
sigDiffKSNonExonPmexD <- list()
sigDiffMWUNonExonPmexD <- list()
sigDiffCVMNonExonPmexD <- list()

sigDiffKSNonExonPlatH <- list()
sigDiffMWUNonExonPlatH <- list()
sigDiffCVMNonExonPlatH <- list()
sigDiffKSNonExonPmexH <- list()
sigDiffMWUNonExonPmexH <- list()
sigDiffCVMNonExonPmexH <- list()

sigDiffKSNonExonPlatE <- list()
sigDiffMWUNonExonPlatE <- list()
sigDiffCVMNonExonPlatE <- list()
sigDiffKSNonExonPmexE <- list()
sigDiffMWUNonExonPmexE <- list()
sigDiffCVMNonExonPmexE <- list()


#loop through simulations
for(i in 1:1000){
  #get simulated and empirical data (D)
  simDatLatD <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                  datNonExonAll$iter == i &
                                  datNonExonAll$pop == "Plat"),8]
  empDatLatD <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                  datNonExonAll$pop == "Plat"),8]
  simDatMexD <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                     datNonExonAll$iter == i &
                                     datNonExonAll$pop == "Pmex"),8]
  empDatMexD <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                     datNonExonAll$pop == "Pmex"),8]
  
  #two sided ks test (D)
  pValKSLatD <- ks.test(empDatLatD,simDatLatD)[2]
  pvalMWULatD <- wilcox.test(empDatLatD,simDatLatD)[3]
  pvalCVMLatD <- cvm_test(empDatLatD,simDatLatD,nboots=1000)[2]
  pValKSMexD <- ks.test(empDatMexD,simDatMexD)[2]
  pvalMWUMexD <- wilcox.test(empDatMexD,simDatMexD)[3]
  pvalCVMMexD <- cvm_test(empDatMexD,simDatMexD,nboots=1000)[2]
  sigDiffKSNonExonPlatD[[length(sigDiffKSNonExonPlatD) + 1]] <- pValKSLatD
  sigDiffMWUNonExonPlatD[[length(sigDiffMWUNonExonPlatD) + 1]] <- pvalMWULatD
  sigDiffCVMNonExonPlatD[[length(sigDiffCVMNonExonPlatD) + 1]] <- pvalCVMLatD
  sigDiffKSNonExonPmexD[[length(sigDiffKSNonExonPmexD) + 1]] <- pValKSMexD
  sigDiffMWUNonExonPmexD[[length(sigDiffMWUNonExonPmexD) + 1]] <- pvalMWUMexD
  sigDiffCVMNonExonPmexD[[length(sigDiffCVMNonExonPmexD) + 1]] <- pvalCVMMexD
  
  #get simulated and empirical data (H)
  simDatLatH <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                     datNonExonAll$iter == i &
                                     datNonExonAll$pop == "Plat"),12]
  empDatLatH <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                     datNonExonAll$pop == "Plat"),12]
  simDatMexH <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                     datNonExonAll$iter == i &
                                     datNonExonAll$pop == "Pmex"),12]
  empDatMexH <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                     datNonExonAll$pop == "Pmex"),12]
  
  #two sided ks test (H)
  pValKSLatH <- ks.test(empDatLatH,simDatLatH)[2]
  pvalMWULatH <- wilcox.test(empDatLatH,simDatLatH)[3]
  pvalCVMLatH <- cvm_test(empDatLatH,simDatLatH,nboots=1000)[2]
  pValKSMexH <- ks.test(empDatMexH,simDatMexH)[2]
  pvalMWUMexH <- wilcox.test(empDatMexH,simDatMexH)[3]
  pvalCVMMexH <- cvm_test(empDatMexH,simDatMexH,nboots=1000)[2]
  sigDiffKSNonExonPlatH[[length(sigDiffKSNonExonPlatH) + 1]] <- pValKSLatH
  sigDiffMWUNonExonPlatH[[length(sigDiffMWUNonExonPlatH) + 1]] <- pvalMWULatH
  sigDiffCVMNonExonPlatH[[length(sigDiffCVMNonExonPlatH) + 1]] <- pvalCVMLatH
  sigDiffKSNonExonPmexH[[length(sigDiffKSNonExonPmexH) + 1]] <- pValKSMexH
  sigDiffMWUNonExonPmexH[[length(sigDiffMWUNonExonPmexH) + 1]] <- pvalMWUMexH
  sigDiffCVMNonExonPmexH[[length(sigDiffCVMNonExonPmexH) + 1]] <- pvalCVMMexH
  
  #get simulated and empirical data (E)
  simDatLatE <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                     datNonExonAll$iter == i &
                                     datNonExonAll$pop == "Plat"),15]
  empDatLatE <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                     datNonExonAll$pop == "Plat"),15]
  simDatMexE <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                     datNonExonAll$iter == i &
                                     datNonExonAll$pop == "Pmex"),15]
  empDatMexE <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                     datNonExonAll$pop == "Pmex"),15]
  
  #two sided ks test (E)
  pValKSLatE <- ks.test(empDatLatE,simDatLatE)[2]
  pvalMWULatE <- wilcox.test(empDatLatE,simDatLatE)[3]
  pvalCVMLatE <- cvm_test(empDatLatE,simDatLatE,nboots=1000)[2]
  pValKSMexE <- ks.test(empDatMexE,simDatMexE)[2]
  pvalMWUMexE <- wilcox.test(empDatMexE,simDatMexE)[3]
  pvalCVMMexE <- cvm_test(empDatMexE,simDatMexE,nboots=1000)[2]
  sigDiffKSNonExonPlatE[[length(sigDiffKSNonExonPlatE) + 1]] <- pValKSLatE
  sigDiffMWUNonExonPlatE[[length(sigDiffMWUNonExonPlatE) + 1]] <- pvalMWULatE
  sigDiffCVMNonExonPlatE[[length(sigDiffCVMNonExonPlatE) + 1]] <- pvalCVMLatE
  sigDiffKSNonExonPmexE[[length(sigDiffKSNonExonPmexE) + 1]] <- pValKSMexE
  sigDiffMWUNonExonPmexE[[length(sigDiffMWUNonExonPmexE) + 1]] <- pvalMWUMexE
  sigDiffCVMNonExonPmexE[[length(sigDiffCVMNonExonPmexE) + 1]] <- pvalCVMMexE
}


#calculate prop sig nonexon (D)
sigDiffKSNonExonPlatD <- unlist(sigDiffKSNonExonPlatD)
sigDiffMWUNonExonPlatD <- unlist(sigDiffMWUNonExonPlatD)
sigDiffCVMNonExonPlatD <- unlist(sigDiffCVMNonExonPlatD)
sigDiffKSNonExonPmexD <- unlist(sigDiffKSNonExonPmexD)
sigDiffMWUNonExonPmexD <- unlist(sigDiffMWUNonExonPmexD)
sigDiffCVMNonExonPmexD <- unlist(sigDiffCVMNonExonPmexD)
propSigLatDKS <- sum(p.adjust(sigDiffKSNonExonPlatD,method="BH")<=0.05)/1000
propSigLatDMWU <- sum(p.adjust(sigDiffMWUNonExonPlatD,method="BH")<=0.05)/1000
propSigLatDCVM <- sum(p.adjust(sigDiffCVMNonExonPlatD,method="BH")<=0.05)/1000
propSigMexDKS <- sum(p.adjust(sigDiffKSNonExonPmexD,method="BH")<=0.05)/1000
propSigMexDMWU <- sum(p.adjust(sigDiffMWUNonExonPmexD,method="BH")<=0.05)/1000
propSigMexDCVM <- sum(p.adjust(sigDiffCVMNonExonPmexD,method="BH")<=0.05)/1000

#calculate prop sig nonexon (H)
sigDiffKSNonExonPlatH <- unlist(sigDiffKSNonExonPlatH)
sigDiffMWUNonExonPlatH <- unlist(sigDiffMWUNonExonPlatH)
sigDiffCVMNonExonPlatH <- unlist(sigDiffCVMNonExonPlatH)
sigDiffKSNonExonPmexH <- unlist(sigDiffKSNonExonPmexH)
sigDiffMWUNonExonPmexH <- unlist(sigDiffMWUNonExonPmexH)
sigDiffCVMNonExonPmexH <- unlist(sigDiffCVMNonExonPmexH)
propSigLatHKS <- sum(p.adjust(sigDiffKSNonExonPlatH,method="BH")<=0.05)/1000
propSigLatHMWU <- sum(p.adjust(sigDiffMWUNonExonPlatH,method="BH")<=0.05)/1000
propSigLatHCVM <- sum(p.adjust(sigDiffCVMNonExonPlatH,method="BH")<=0.05)/1000
propSigMexHKS <- sum(p.adjust(sigDiffKSNonExonPmexH,method="BH")<=0.05)/1000
propSigMexHMWU <- sum(p.adjust(sigDiffMWUNonExonPmexH,method="BH")<=0.05)/1000
propSigMexHCVM <- sum(p.adjust(sigDiffCVMNonExonPmexH,method="BH")<=0.05)/1000

#calculate prop sig nonexon (E)
sigDiffKSNonExonPlatE <- unlist(sigDiffKSNonExonPlatE)
sigDiffMWUNonExonPlatE <- unlist(sigDiffMWUNonExonPlatE)
sigDiffCVMNonExonPlatE <- unlist(sigDiffCVMNonExonPlatE)
sigDiffKSNonExonPmexE <- unlist(sigDiffKSNonExonPmexE)
sigDiffMWUNonExonPmexE <- unlist(sigDiffMWUNonExonPmexE)
sigDiffCVMNonExonPmexE <- unlist(sigDiffCVMNonExonPmexE)
propSigLatEKS <- sum(p.adjust(sigDiffKSNonExonPlatE,method="BH")<=0.05)/1000
propSigLatEMWU <- sum(p.adjust(sigDiffMWUNonExonPlatE,method="BH")<=0.05)/1000
propSigLatECVM <- sum(p.adjust(sigDiffCVMNonExonPlatE,method="BH")<=0.05)/1000
propSigMexEKS <- sum(p.adjust(sigDiffKSNonExonPmexE,method="BH")<=0.05)/1000
propSigMexEMWU <- sum(p.adjust(sigDiffMWUNonExonPmexE,method="BH")<=0.05)/1000
propSigMexECVM <- sum(p.adjust(sigDiffCVMNonExonPmexE,method="BH")<=0.05)/1000

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
            file="output/results/LOHDStatsNonExon.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)
write.table(propSigDfH, 
            file="output/results/LOHHStatsNonExon.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)
write.table(propSigDfE,
            file="output/results/LOHEStatsNonExon.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)

