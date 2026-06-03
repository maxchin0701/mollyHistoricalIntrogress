#### LIBRARIES ####
library("clusterProfiler")
library("AnnotationHub")
library("biomaRt")

#### LOAD IN ORGDB ####
ah <- AnnotationHub()
PForOrgDB <- query(ah,c("Poecilia formosa", "OrgDB"))[[1]]

AnnotationHub::query(ah,c("Poecilia","orgdb"))

#### LOAD IN DATA ####
genes <- read.delim(paste0("../output/LOHRegionsSharedCombined/LOHSharedGenes.tsv"),
                    sep="\t",row.names = NULL,header = F)
colnames(genes) <- c("gene","pepID","symbol")
#genes$symbol <- tolower(genes$symbol)

#### GET GENE ID FROM PEPTIDE ID ####
ensembl <- useEnsembl(biomart='ensembl',dataset="pformosa_gene_ensembl")

for(i in 1:nrow(genes)){
  genes$geneID[i] <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
        filters = "ensembl_peptide_id",
        values = genes$pepID[i],
        mart = ensembl)[1,2]
}

#### CONVERT DATA ####
genesEntrezID <- bitr(genes$geneID, fromType = "ENSEMBL",
                             toType = c("ENTREZID"),
                             OrgDb = PForOrgDB)

#### GENE OVERREPRESENTATION (ENRICHMENT) ####
outEnrichGOBP <- enrichGO(gene=as.numeric(genesEntrezID$ENTREZID),
                        OrgDb=PForOrgDB,
                        #universe=names(geneList),
                        keyType = "ENTREZID",
                        ont="BP",
                        pvalueCutoff=0.0125,
                        pAdjustMethod="BH",
                        qvalueCutoff = 0.05,
                        readable=T)

outEnrichGOCC <- enrichGO(gene=as.numeric(genesEntrezID$ENTREZID),
                        OrgDb=PForOrgDB,
                        #universe=names(geneList),
                        keyType = "ENTREZID",
                        ont="CC",
                        pvalueCutoff=0.0125,
                        pAdjustMethod="BH",
                        qvalueCutoff = 0.05,
                        readable=T)

outEnrichGOMF <- enrichGO(gene=as.numeric(genesEntrezID$ENTREZID),
                          OrgDb=PForOrgDB,
                          #universe=names(geneList),
                          keyType = "ENTREZID",
                          ont="MF",
                          pvalueCutoff=0.0125,
                          pAdjustMethod="BH",
                          qvalueCutoff = 0.05,
                          readable=T)



#check results
outDFBP <- outEnrichGOBP@result
outDFCC <- outEnrichGOCC@result
outDFMF <- outEnrichGOMF@result


