#### LIBRARY ####
library(coda)
library(ggplot2)

#### LOAD IN DATA ####
LOHRegions <- read.delim(paste0("../../SNPBasedLOH/output/LOHRegions/LOHRegionsAll.tsv"),
                         sep="\t")

#### SEPARATE OUT OVERLAPPING AND NONOVERLAPPING ####
LOHRegionsIsolated <- LOHRegions[which(LOHRegions$group <= 5 & 
                                         LOHRegions$intSize > 1000),]
LOHRegionsShared <- LOHRegions[which(LOHRegions$group >= 18 &
                                       LOHRegions$intSize > 1000),]

#### SAMPLE AND COMPARE ####
sigDiff <- list()
for(i in 1:1000){
  #sample intSizes
  curIsolatedSample <- sample(LOHRegionsIsolated$intSize,size=65,replace=T)
  
  #two sided ks test
  pVal<- ks.test(LOHRegionsShared$intSize,curIsolatedSample)[2]
  sigDiff[[length(sigDiff) + 1]] <- pVal 
}
sigDiff <- unlist(sigDiff)
sum(sigDiff<=0.05)/1000
