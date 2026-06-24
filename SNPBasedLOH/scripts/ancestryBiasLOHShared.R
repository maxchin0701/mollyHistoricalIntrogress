# This script runs a binomial test comparing distribution of ancestries across
# shared LOH tracts. 
# Output: 
#     LOHAncBinomTest.tsv: tsv with test statistics

#### LIBRARY ####
library(coda)
library(broom)

#### READ IN PARAMS ####
args <- commandArgs(trailingOnly = TRUE)
inLOHAnc <- args[1]
outBinomSummary <- args[2]

#### LOAD IN DATA ####
LOHAncShared <- read.delim(paste0(inLOHAnc),
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
            file=outBinomSummary, 
            quote=FALSE, sep='\t',row.names=FALSE)
