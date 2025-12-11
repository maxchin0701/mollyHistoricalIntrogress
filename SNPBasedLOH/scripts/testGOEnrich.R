#### LIBRARIES ####
library("clusterProfiler")
library("AnnotationHub")

#### LOAD IN ORGDB ####
ah <- AnnotationHub()
# PForOrgDB <- query(ah,c("Poecilia formosa", "OrgDB"))
HSapOrgDB <- query(ah,c("Homo sapiens", "OrgDB"))[[1]]


#### LOAD IN DATA ####
data(geneList,package="DOSE")
genes <- names(geneList[abs(geneList)>2])

#### RUN GROUP GO #### 
outGroupGO <- groupGO(genes,
                      OrgDb = HSapOrgDB,
                      ont="BP",
                      level=3)

#### GENE OVERREPRESENTATION (ENRICHMENT) ####
outEnrichGO <- enrichGO(gene=genes,
                        OrgDb=HSapOrgDB,
                        universe=names(geneList),
                        ont="BP",
                        pvalueCutoff=0.0125,
                        pAdjustMethod="BH",
                        qvalueCutoff = 0.5,
                        readable=T)

#check results
outEnrichGO@result

#### REDO WITH GENE SYMBOLS
geneSymbols <- bitr(genes, fromType = "ENTREZID",
                toType = c("SYMBOL"),
                OrgDb = HSapOrgDB)

outEnrichGOSymbol <- enrichGO(gene=geneSymbols$SYMBOL,
                        OrgDb=HSapOrgDB,
                        universe=names(geneList),
                        ont="BP",
                        pvalueCutoff=0.0125,
                        pAdjustMethod="BH",
                        qvalueCutoff = 0.5,
                        keyType = "SYMBOL",
                        readable=T)










