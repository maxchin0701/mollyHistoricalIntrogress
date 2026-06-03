#### LIBRARY ####
library(ggplot2)

#### LOAD IN DATA ####
LOHRegions <- read.delim(paste0("../../SNPBasedLOH/output/LOHRegions/LOHRegionsAll.tsv"),
                         sep="\t")
LOHAnc <- read.delim(paste0("../../SNPBasedLOH/output/LOHAnc/LOHAncAll.tsv"),
                         sep="\t")

#get freq
LOHAnc <- cbind(LOHAnc,LOHRegions$group)

#new rownames
colnames(LOHAnc) <- c("chr","start","end","anc","freq")

#interval size
LOHAnc$intSize <- LOHAnc$end - LOHAnc$start

#remove mixed ancestr
LOHAnc <- LOHAnc[-which(LOHAnc$anc == "mix"),]

LOHAnc <- LOHAnc[which(LOHAnc$intSize > 1000),]

#set cols
cols <- colorRampPalette(c("#0000FF","#F37221"))(2)

#### PLOT ####

#plot
pdf(file=paste0("suppFigAncScatter.pdf"),
    width=6.49,height=5.15)
ggplot()+
  geom_jitter(aes(y=LOHAnc$intSize,x=as.factor(LOHAnc$freq),color=as.factor(LOHAnc$anc)),
              alpha=0.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))+
  scale_y_continuous("Shared homozygosity tract size (bp)")+
  #scale_y_continuous("Log likelihood",limits=c(-1550000,-150000))+
  scale_x_discrete("Frequency")+
  scale_color_manual(values=cols,labels=c("plat" = "P. latipinna",
                                          "pmex" = "P. mexicana"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.text.x = element_blank(),
        #axis.ticks.x=element_blank(),
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
dev.off()






