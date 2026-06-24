#### GENERATE SCATTERPLOT FOR TRACT LENGTHS ####
#### PACKAGES ####
library(ggplot2)

#### LOAD IN DATA ####
LOHSharedLong <- read.delim("../../SNPBasedLOH/output/LOHRegionsSharedCombined/LOHRegionsSharedCombined.bed",
                          sep='\t',header=F)
chrs <- read.delim(paste0("../../SNPBasedLOH/data/chrIndex.tsv"),
                   sep="\t",row.names = NULL,header = T)

#cutoff GQ
parentalSpecies <- "LM"
cutoffGQ <- 10
cutoffDP <- 6

#set up dfs to store
LOHAnc <- as.data.frame(matrix(NA,nrow=0,ncol=4))

for(i in unique(chrs$CHR)){
  curLOHAnc <- read.delim(paste0("../../SNPBasedLOH/output/LOHAnc/",i,"LOHAnc_GQ",cutoffGQ,"_DP",cutoffDP,"_",parentalSpecies,".tsv"),sep="\t")
  LOHAnc <- rbind(LOHAnc,curLOHAnc)
}

#### PREPARE LOH SHARED LONG ANC DF ####
#df
LOHAncSharedLongIndex <- c()

#loop to get rows
for(i in 1:nrow(LOHAnc)){
  if(LOHAnc[i,2] %in% LOHSharedLong[,2] &&
     LOHAnc[i,3] %in% LOHSharedLong[,3]){
    LOHAncSharedLongIndex <- c(LOHAncSharedLongIndex,i)
  }
}

#subset
LOHAncSharedLong <- LOHAnc[LOHAncSharedLongIndex,]

#interval size
LOHAncSharedLong$intSize <- LOHAncSharedLong$end - LOHAncSharedLong$start

#### PLOT ####
cols <- colorRampPalette(c("#0000FF","#F37221"))(2)

#plot
ggplot()+
  geom_jitter(aes(y=LOHAncSharedLong$intSize,x=as.factor(LOHAncSharedLong$group),color=as.factor(LOHAncSharedLong$group)),
              alpha=0.5,width=0.2,height=0)+
  scale_y_continuous("Shared homozygosity tract size")+
  #scale_y_continuous("Log likelihood",limits=c(-1550000,-150000))+
  scale_x_discrete("Ancestry retained across tract")+
  scale_color_manual(values=cols,labels=c("plat" = "P. latipinna",
                                          "pmex" = "P. mexicana"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 7.5),
        axis.title = element_text(face = "bold",
                                  size = 7.5),
        legend.title = element_blank(),
        legend.position = c(0.8,0.9),
        legend.key = element_blank(),
        legend.box = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size=7.5),
        legend.key.height= unit(0.1,"inch"),
        legend.key.width = unit(0.1,"inch"),
        legend.spacing.y = unit(0.03,"inch"))






