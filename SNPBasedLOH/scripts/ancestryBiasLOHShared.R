#### LIBRARY ####
library(coda)
library(broom)

#### LOAD IN DATA ####
LOHAncShared <- read.delim(paste0("output/LOHRegionsSharedCombined/LOHAncSharedCombined.bed"),
                         sep="\t",header=F)

#### SIMULATE WITH BINOM ####
sum(rbinom(1000,65,prob=0.5) > length(which(LOHAncShared[,4] == "pmex")))/1000

#### BINOMIAL TEST ####
binomOut <- binom.test(length(which(LOHAncShared[,4] == "pmex")),65)

#### SAVE ####
#convert to neat format
binomOutTidy <- tidy(binomOut)

#save
write.table(binomOutTidy, 
            file="output/results/LOHAncBinomTest.tsv", 
            quote=FALSE, sep='\t',row.names=FALSE)
