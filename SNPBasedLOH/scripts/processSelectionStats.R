#### LOAD PACKAGES ####
library(twosamples)

##### LOAD IN DATA ####
#empirical nonexon
datNonExonEmp <- read.delim(paste0("../output/selectionStats/LOHRegionsEmpStatsNonExon.tsv"),
                   sep="\t",row.names = NULL,header = T)
datNonExonEmp$iter <- 0
datNonExonEmp$class <- "Empirical"

#empirical exon
datExonEmp <- read.delim(paste0("../output/selectionStats/LOHRegionsEmpStatsExon.tsv"),
                            sep="\t",row.names = NULL,header = T)
datExonEmp$iter <- 0
datExonEmp$class <- "Empirical"

#sim
datSimNonExonList <- list()
datSimExonList <- list()
for(i in 1:1000){
  #non exon
  curDatNonExonSim <- read.delim(paste0("../output/selectionStats/simStatsOut/simRegionsStatsNonExon_",i,".tsv"),
                          sep="\t",row.names = NULL,header = T)
  curDatNonExonSim$iter <- i
  datSimNonExonList[[i]] <- curDatNonExonSim
  
  #exon
  curDatExonSim <- read.delim(paste0("../output/selectionStats/simStatsOut/simRegionsStatsExon_",i,".tsv"),
                                 sep="\t",row.names = NULL,header = T)
  curDatExonSim$iter <- i
  datSimExonList[[i]] <- curDatExonSim
}

#bind to form df
datNonExonSim <- do.call(rbind,datSimNonExonList)
datNonExonSim$class <- "Simulated"

#bind to form df
datExonSim <- do.call(rbind,datSimExonList)
datExonSim$class <- "Simulated"

#combine empirical and simulated
datNonExonAll <- rbind(datNonExonEmp,datNonExonSim)
datExonAll <- rbind(datExonEmp,datExonSim)

#### TEST SIGNIFICANCE ####
#non exon
sigDiffKSNonExonPlat <- list()
sigDiffMWUNonExonPlat <- list()
sigDiffCVMNonExonPlat <- list()
sigDiffKSNonExonPmex <- list()
sigDiffMWUNonExonPmex <- list()
sigDiffCVMNonExonPmex <- list()
for(i in 1:1000){
  #get simulated and empirical data
  simDatLat <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                  datNonExonAll$iter == i &
                                  datNonExonAll$pop == "Plat"),15]
  empDatLat <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                  datNonExonAll$pop == "Plat"),15]
  simDatMex <- datNonExonAll[which(datNonExonAll$class == "Simulated" & 
                                     datNonExonAll$iter == i &
                                     datNonExonAll$pop == "Pmex"),15]
  empDatMex <- datNonExonAll[which(datNonExonAll$class == "Empirical" & 
                                     datNonExonAll$pop == "Pmex"),15]
  
  #two sided ks test
  pValKSLat <- ks.test(empDatLat,simDatLat)[2]
  pvalMWULat <- wilcox.test(empDatLat,simDatLat)[3]
  pvalCVMLat <- cvm_test(empDatLat,simDatLat,nboots=1000)[2]
  pValKSMex <- ks.test(empDatMex,simDatMex)[2]
  pvalMWUMex <- wilcox.test(empDatMex,simDatMex)[3]
  pvalCVMMex <- cvm_test(empDatMex,simDatMex,nboots=1000)[2]
  sigDiffKSNonExonPlat[[length(sigDiffKSNonExonPlat) + 1]] <- pValKSLat
  sigDiffMWUNonExonPlat[[length(sigDiffMWUNonExonPlat) + 1]] <- pvalMWULat
  sigDiffCVMNonExonPlat[[length(sigDiffCVMNonExonPlat) + 1]] <- pvalCVMLat
  sigDiffKSNonExonPmex[[length(sigDiffKSNonExonPmex) + 1]] <- pValKSMex
  sigDiffMWUNonExonPmex[[length(sigDiffMWUNonExonPmex) + 1]] <- pvalMWUMex
  sigDiffCVMNonExonPmex[[length(sigDiffCVMNonExonPmex) + 1]] <- pvalCVMMex
}

#exon
sigDiffKSExonPlat <- list()
sigDiffMWUExonPlat <- list()
sigDiffCVMExonPlat <- list()
sigDiffKSExonPmex <- list()
sigDiffMWUExonPmex <- list()
sigDiffCVMExonPmex <- list()
for(i in 1:1000){
  #get simulated and empirical data
  simDatLat <- datExonAll[which(datExonAll$class == "Simulated" & 
                                     datExonAll$iter == i &
                                     datExonAll$pop == "Plat"),15]
  empDatLat <- datExonAll[which(datExonAll$class == "Empirical" & 
                                     datExonAll$pop == "Plat"),15]
  simDatMex <- datExonAll[which(datExonAll$class == "Simulated" & 
                                     datExonAll$iter == i &
                                     datExonAll$pop == "Pmex"),15]
  empDatMex <- datExonAll[which(datExonAll$class == "Empirical" & 
                                     datExonAll$pop == "Pmex"),15]
  
  #two sided ks test
  pValKSLat <- ks.test(empDatLat,simDatLat)[2]
  pvalMWULat <- wilcox.test(empDatLat,simDatLat)[3]
  pvalCVMLat <- cvm_test(empDatLat,simDatLat,nboots=1000)[2]
  pValKSMex <- ks.test(empDatMex,simDatMex)[2]
  pvalMWUMex <- wilcox.test(empDatMex,simDatMex)[3]
  pvalCVMMex <- cvm_test(empDatMex,simDatMex,nboots=1000)[2]
  sigDiffKSExonPlat[[length(sigDiffKSExonPlat) + 1]] <- pValKSLat
  sigDiffMWUExonPlat[[length(sigDiffMWUExonPlat) + 1]] <- pvalMWULat
  sigDiffCVMExonPlat[[length(sigDiffCVMExonPlat) + 1]] <- pvalCVMLat
  sigDiffKSExonPmex[[length(sigDiffKSExonPmex) + 1]] <- pValKSMex
  sigDiffMWUExonPmex[[length(sigDiffMWUExonPmex) + 1]] <- pvalMWUMex
  sigDiffCVMExonPmex[[length(sigDiffCVMExonPmex) + 1]] <- pvalCVMMex
}

#calculate prop sig nonexon
sigDiffKSNonExonPlat <- unlist(sigDiffKSNonExonPlat)
sigDiffMWUNonExonPlat <- unlist(sigDiffMWUNonExonPlat)
sigDiffCVMNonExonPlat <- unlist(sigDiffCVMNonExonPlat)
sigDiffKSNonExonPmex <- unlist(sigDiffKSNonExonPmex)
sigDiffMWUNonExonPmex <- unlist(sigDiffMWUNonExonPmex)
sigDiffCVMNonExonPmex <- unlist(sigDiffCVMNonExonPmex)
sum(p.adjust(sigDiffKSNonExonPlat,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffMWUNonExonPlat,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffCVMNonExonPlat,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffKSNonExonPmex,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffMWUNonExonPmex,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffCVMNonExonPmex,method="BH")<=0.05)/1000


#calculate prop sig exon
sigDiffKSExonPlat <- unlist(sigDiffKSExonPlat)
sigDiffMWUExonPlat <- unlist(sigDiffMWUExonPlat)
sigDiffCVMExonPlat <- unlist(sigDiffCVMExonPlat)
sigDiffKSExonPmex <- unlist(sigDiffKSExonPmex)
sigDiffMWUExonPmex <- unlist(sigDiffMWUExonPmex)
sigDiffCVMExonPmex <- unlist(sigDiffCVMExonPmex)
sum(p.adjust(sigDiffKSExonPlat,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffMWUExonPlat,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffCVMExonPlat,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffKSExonPmex,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffMWUExonPmex,method="BH")<=0.05)/1000
sum(p.adjust(sigDiffCVMExonPmex,method="BH")<=0.05)/1000

# #### PLOT ####
ggplot(datNonExonAll[which(datNonExonAll$pop == "Plat"),], aes(x=H,group=iter,color=class,alpha=as.factor(class))) +
  geom_line(stat="density")+
  scale_alpha_manual(values = c(1,0.05))

