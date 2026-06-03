#### CALCULATE AIC VALUES FOR ALL MODELS ####
#### LOAD IN DATA ####
#get model from snakemake
#take in command line arguments
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]

#load in bestLhood file
dat <- read.delim(paste0("data/",model,"/bestIter/",model,"/",model,".bestlhoods"),sep="\t")

#### CALCULATE AIC ####
k <- ncol(dat) - 2

AIC <- 2*k - 2 * (dat$MaxEstLhood/log10(exp(1)))

##### SAVE AIC ####
write.table(cbind(dat[1,(k+1):ncol(dat)],AIC),
            paste0("data/",model,"/bestIter/",model,".AIC"),
            row.names = F,col.names = T,sep = "\t",quote = F)



