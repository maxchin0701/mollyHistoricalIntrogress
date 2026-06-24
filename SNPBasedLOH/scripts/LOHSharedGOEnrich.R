# This script performs GO enrichment analysis on genes intersecting empirical 
# shared tracts.
# output:
#   LOHGenesGOEnrichBP.tsv: tsv with Biological Process GO enrichment results 
#   LOHGenesGOEnrichCC.tsv: tsv with Cellular Component GO enrichment results 
#   LOHGenesGOEnrichMF.tsv: tsv with Molecular Function GO enrichment results 

#### LIBRARIES ####
library("clusterProfiler")
library("AnnotationHub")
library("biomaRt")

#### READ IN PARAMS ####
args <- commandArgs(trailingOnly = TRUE)
intersectGenes <- args[1]
outBP <- args[2]
outCC <- args[3]
outMF <- args[4]

#### LOAD IN ORGDB ####
ah <- AnnotationHub()
PForOrgDB <- query(ah,c("Poecilia formosa", "OrgDB"))[[1]]

AnnotationHub::query(ah,c("Poecilia","orgdb"))

#### LOAD IN DATA ####
genes <- read.delim(paste0(intersectGenes),
                    sep="\t",row.names = NULL,header = F)
colnames(genes) <- c("gene","pepID","symbol")
#genes$symbol <- tolower(genes$symbol)

#### GET GENE ID FROM PEPTIDE ID ####
ensembl <- useEnsembl(biomart='ensembl',dataset="pformosa_gene_ensembl")

for(i in 1:nrow(genes)){
  print(i)
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

#### SAVE ####
#BP
write.table(outDFBP, 
            file=outBP, 
            quote=FALSE, sep='\t',row.names=FALSE)
#CC
write.table(outDFCC, 
            file=outCC, 
            quote=FALSE, sep='\t',row.names=FALSE)
#MF
write.table(outDFMF, 
            file=outMF, 
            quote=FALSE, sep='\t',row.names=FALSE)


