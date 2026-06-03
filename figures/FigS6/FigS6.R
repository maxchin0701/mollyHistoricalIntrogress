#### LOAD PACKAGES ####
library(ggplot2)
library(coda)

#### LOAD IN DATA ####
#define models
# models <- c("noGeneFlow",
#          "modernGeneFlow",
#          "ancestralModernGeneFlowCons",
#          "ancestralModernGeneFlowShift")

models <- c("noGeneFlow",
            "modernGeneFlow",
            "ancestralGeneFlow",
            "ancestralModernGeneFlowCons",
            "ancestralModernGeneFlowShift")


modelsNeat <- c("No gene flow",
                "Modern gene flow",
                "Ancestral gene flow",
                "Modern + ancestral gene flow \n(constant rate)",
                "Modern + ancestral gene flow \n(rate shift)")

#load in lhoods
lhoods <- as.data.frame(matrix(NA,nrow=0,ncol=2))

for(i in models){
  curLhood <- read.delim(paste0("../../demMollyHistIntrogress/outputs/lhoods/",i,".lhoods"),
                         sep="\t",row.names = NULL,header = F)[1:100,2]
  lhoods <- rbind(lhoods,
                  cbind(curLhood,
                        rep(i,100)))
}
colnames(lhoods) <- c("lhood","model")
lhoods$lhood <- as.numeric(lhoods$lhood)
lhoods$model <- factor(lhoods$model,levels=unique(lhoods$model))
levels(lhoods$model) <- modelsNeat

#subset to null + alternative
lhoodsBest <- lhoods[which(lhoods$model=="Modern gene flow" |
                             lhoods$model=="Ancestral gene flow" |
                             lhoods$model=="Modern + ancestral gene flow \n(rate shift)" |
                             lhoods$model=="Modern + ancestral gene flow \n(constant rate)"),]
lhoodsBest$model <- factor(lhoodsBest$model,levels=unique(lhoodsBest$model))

#### GET HPD ####
#set up df
lhoodsBestHPD <- as.data.frame(matrix(data=NA,nrow=0,ncol=3))

#loop through models
for(i in 1:length(levels(lhoodsBest$model))){
  lhoodsBestHPD <- rbind(lhoodsBestHPD,
                         cbind(HPDinterval(as.mcmc(lhoodsBest$lhood[
                           which(lhoodsBest$model == levels(lhoodsBest$model)[i])]))[1],
                           -0.5e-4,
                           levels(lhoodsBest$model)[i]))
  lhoodsBestHPD <- rbind(lhoodsBestHPD,
                         cbind(HPDinterval(as.mcmc(lhoodsBest$lhood[
                           which(lhoodsBest$model == levels(lhoodsBest$model)[i])]))[2],
                           -0.5e-4,
                           levels(lhoodsBest$model)[i]))
}

colnames(lhoodsBestHPD) <- c("x","y","model")
lhoodsBestHPD$model <- factor(lhoodsBestHPD$model,
                              levels=unique(lhoodsBestHPD$model))
#### PLOT ####
cols <- colorRampPalette(c("#F37221","#0000FF"))(4)

ggplot()+
  geom_density(aes(x=lhoodsBest$lhood,fill=lhoodsBest$model),alpha=0.5)+
  geom_line(mapping=aes(x=as.numeric(lhoodsBestHPD$x),y=as.numeric(lhoodsBestHPD$y),
                        color=lhoodsBestHPD$model),
            show.legend = F,
            linewidth=1,alpha=0.5)+
  scale_y_continuous("Density")+
  xlab("Log10 likelihood")+
  xlim(-455000,-442000)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.5),
        axis.title = element_text(face = "bold",
                                  size = 7.5),
        legend.title = element_blank(),
        legend.position = c(0.7,0.9),
        legend.key = element_blank(),
        legend.box = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size=7.5),
        legend.key.height= unit(0.1,"inch"),
        legend.key.width = unit(0.1,"inch"),
        legend.spacing.y = unit(0.03,"inch"))
