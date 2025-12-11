#### GENERATE VENN DIAGRAM COMPARING ANC AND MODERN LOH ####
#### PACKAGES ####
library(VennDiagram)

#### LOAD IN DATA ####
LOHSharedSNP<- read.delim("../output/LOHRegionsSharedCombined/LOHRegionsSharedCombined.bed",
           sep='\t',header=F)
AncUnsupported <- c("chr3_24298801_24309560",
                    "chr18_26967287_26981637",
                    "chr18_6101577_6102840",
                    "chr17_5363790_5372962")

#### CONVERT TO LIST ####
LOHSharedAll <- list(paste(LOHSharedSNP$V1,LOHSharedSNP$V2,LOHSharedSNP$V3,sep="_"))
LOHSharedAll[[2]] <- LOHSharedAll[[1]][-which(LOHSharedAll[[1]] %in% AncUnsupported)]

names(LOHSharedAll) <- c("SNP","AGR")

#### PLOT ####

venn.plot <- venn.diagram(x=LOHSharedAll,
             filename=NULL,
             disable.logging = T,scaled=T,ext.text=F,cex=c(0.75,3,5),
             fill=c("#F37221","#0000FF"),alpha=0.5,
             cat.col=c("#F37221","#7A3990"),cat.cex=c(3,3), 
             #cat.default.pos = "outer",
             cat.just = list(c(0.4,-0.5),c(0.4,-0.5))
)#,ext.line.lty = 0,ext.dist=-0.250)


grid.draw(venn.plot)
