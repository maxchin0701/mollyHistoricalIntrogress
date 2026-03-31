#### GENERATE DENSITY PLOT FOR TRACT LENGTHS ####
#### PACKAGES ####
library(ggplot2)

#### LOAD IN DATA ####
LOHRegions <- read.delim(paste0("../../SNPBasedLOH/output/LOHRegions/LOHRegionsAll.tsv"),
                         sep="\t")
LOHRegionsLong <- LOHRegions[which(LOHRegions$intSize > 1000 & (LOHRegions$group >= 18 | LOHRegions$group <= 5)),]
chrs <- read.delim(paste0("../../SNPBasedLOH/data/chrIndex.tsv"),
                   sep="\t",row.names = NULL,header = T)

#### GET SHARED VS ISOLATED ####
LOHRegionsLong$freq <- apply(LOHRegionsLong,1,function(x) if(as.numeric(x["group"]) <= 5){
  return("Isolated")
  } else {
    return("Shared")
  })

#### PLOT ####
cols <- colorRampPalette(c("#0000FF","#F37221"))(2)

#plot
ggplot()+
  geom_jitter(aes(y=LOHRegionsLong$intSize,,x=as.factor(LOHRegionsLong$freq),
                  color=as.factor(LOHRegionsLong$freq)),
              alpha=0.5)+
  scale_y_continuous("Homozygosity tract size")+
  #scale_y_continuous("Log likelihood",limits=c(-1550000,-150000))+
  scale_x_discrete("Tract frequency")+
  scale_color_manual(values=cols)+
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
