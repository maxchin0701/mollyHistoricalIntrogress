#### LOAD PACKAGES ####
library(ggplot2)

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

#### PLOT ####

#get colors
cols <- colorRampPalette(c("#F37221","#0000FF"))(5)

ggplot()+
  geom_jitter(aes(y=lhoods$lhood,x=lhoods$model,color=lhoods$model),
               alpha=0.5,width=0.2,height=0)+
  scale_y_continuous("Log10 likelihood",limits=c(-670000,-430000))+
  #scale_y_continuous("Log likelihood",limits=c(-1550000,-150000))+
  scale_x_discrete("Model")+
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
        legend.position = c(0.8,0.5),
        legend.key = element_blank(),
        legend.box = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size=5),
        legend.key.height= unit(0.1,"inch"),
        legend.key.width = unit(0.1,"inch"),
        legend.spacing.y = unit(0.03,"inch"))
  
#save as 3.5 x 3.5 in