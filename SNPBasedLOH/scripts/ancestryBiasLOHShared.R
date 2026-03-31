#### LIBRARY ####
library(coda)
library(ggplot2)

#### LOAD IN DATA ####
LOHAncShared <- read.delim(paste0("../../SNPBasedLOH/output/LOHRegionsSharedCombined/LOHAncSharedCombined.bed"),
                         sep="\t",header=F)

#### SIMULATE WITH BINOM ####
sum(rbinom(1000,65,prob=0.5) < length(which(LOHAncShared[,4] == "pmex")))/1000

#### BINOMIAL TEST ####
binom.test(length(which(LOHAncShared[,4] == "pmex")),65)